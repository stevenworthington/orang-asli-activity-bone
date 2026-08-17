###############################################################################
# Static verification of the startup set and the scripts that consume it.
#
# Run with `just check`. Every item is an assertion; the script stops at the
# first failure and exits non-zero, so it can gate a long run.
#
# The point is to prove a change to the source is COMPLETE before spending an
# hour of sampling on it. Assertions here are the ones a reader of a diff
# cannot make for themselves: that no call resolves to an attached package
# without being namespaced, that the registry's 15 specs really carry the
# families and constants the specification names, that the prior scheme
# reproduces a hand computation, and that arguments documented as switches
# actually switch something.
#
# Fits nothing by default. `CHECK_FIT=1 just check` adds a smoke fit of one
# community-level spec, which is the only way to confirm that `re_formula`
# survives the trip through marginaleffects.
#
# Run: Rscript code/_checks/verify-startup.R
###############################################################################



FAILURES <- 0L

ok <- function(what) cat(sprintf("  ok    %s\n", what))

check <- function(what, condition, detail = NULL) {
  if (isTRUE(condition)) {
    ok(what)
  } else {
    FAILURES <<- FAILURES + 1L
    cat(sprintf("  FAIL  %s\n", what))
    if (!is.null(detail)) cat(paste0("          ", detail, collapse = "\n"), "\n")
  }
  invisible(isTRUE(condition))
}

section <- function(title) cat(sprintf("\n== %s ==\n", title))


# ---- The files this check covers ----
#
# The startup set, the nine reported analyses, the figure/table assembly, and
# the `_experiments/` scripts behind a reported artifact. `_targets.R` is
# excluded: it is a targets manifest written in that package's DSL, where
# `tar_target()` reads better bare.

STARTUP_FILES <- here::here("code", "_startup",
                            c("init.R", "functions.R", "packages.R", "options.R",
                              "data.R", "specifications.R", "pipeline-helpers.R"))
ANALYSIS_FILES   <- Sys.glob(here::here("code", "*.R"))
FINAL_FILES   <- Sys.glob(here::here("code", "_final", "*.R"))

# Two of these (effect-size-probability-fig.R, outcome-noise-anchors.R) are
# standalone base-R scripts that never source init.R, so nothing masks anything
# when they run. They are still held to the same rule, because the qualified
# form is correct under every search path and the unqualified form is only
# correct under theirs.
EXPERIMENT_FILES <- here::here("code", "_experiments",
  c("effect-size-probability-fig.R",        # Supplementary Figure 4
    "age-subset-amef.R",                    # Supplementary Figure 8
    "outcome-noise-anchors.R",              # measurement-resolution anchors
    "walking-experiment-rope.R",            # acute walking-experiment reverse-ROPE
    "industrialization-linear-summaries.R", # main-text industrialization summaries
    "_reference-scales.R",                  # shared SD / resolution-floor scales
    "pa-contrast-effects.R",                # main-text magnitude bounds
    "conf-slope-summaries.R"))              # confounder-DAG magnitude bounds

SWEPT_FILES <- c(STARTUP_FILES, ANALYSIS_FILES, FINAL_FILES, EXPERIMENT_FILES)


section("Files parse")

check("all 7 _startup/ files present", length(STARTUP_FILES) == 7 &&
        all(file.exists(STARTUP_FILES)),
      STARTUP_FILES[!file.exists(STARTUP_FILES)])
check("9 analysis scripts present", length(ANALYSIS_FILES) == 9,
      paste("found", length(ANALYSIS_FILES)))

parse_errs <- character(0)
for (f in SWEPT_FILES) {
  e <- tryCatch({ parse(f); NULL }, error = function(e) conditionMessage(e))
  if (!is.null(e)) parse_errs <- c(parse_errs, paste(basename(f), "--", e))
}
check(sprintf("all %d swept files parse", length(SWEPT_FILES)),
      length(parse_errs) == 0, parse_errs)


section("Session loads")

suppressMessages(source(here::here("code", "_startup", "init.R")))
check("source(init.R) completes", exists("model_templates") && exists("dat"))
check("dat has 1,007 rows", nrow(dat) == 1007, paste("got", nrow(dat)))
check("registry holds 15 specs", length(model_templates) == 15,
      paste("got", length(model_templates)))


section("No stale identifiers")

# None of these names should reach the interpreter. A surviving reference is a
# silent wrong answer, not an error: a renamed column resolves to NULL inside a
# formula, `spec$priors` to NULL where a prior object is expected, and a helper
# that no longer exists is looked up in whatever package happens to export a
# matching name.
STALE <- c("ad_steps_1k", "tibia_sos_1k", "steps_1k",
           "spec\\$priors", "HDInterval::hdi",
           "\\bscale_this\\b", "\\bcalc_probs\\b", "\\bapprox_flat\\b",
           "\\bapprox_linearity\\b", "\\bmake_age_cond_amef_panel\\b",
           "\\bpp_check_stats_grouped\\b",
           # the density-based interval primitive; this project uses the
           # sample-based `median_hdci` / `hdci` throughout
           "(?<![a-z_])median_hdi\\b(?!ci)")

# Matched against the DEPARSED code, so comments are excluded. Some of these
# names are legitimately DISCUSSED in prose -- data.R explains why a former
# column name is not recreated, and pipeline-helpers.R explains how `hdci` and
# `hdi` differ. Prose about a name is the documentation working; a live
# reference to one is the bug.
code_of <- function(f) unlist(lapply(parse(f, keep.source = FALSE), deparse))
code <- lapply(SWEPT_FILES, code_of)
names(code) <- SWEPT_FILES
for (pat in STALE) {
  hits <- unlist(lapply(names(code), function(f) {
    i <- grep(pat, code[[f]], perl = TRUE)
    if (length(i)) sprintf("%s: %s", basename(f), trimws(code[[f]][i]))
  }))
  check(sprintf("no live reference to /%s/", pat), length(hits) == 0, hits)
}


section("Namespacing is complete")

# Walk each swept file's AST and collect every call head that is a bare symbol.
# A head written as `pkg::fun` parses to a call to `::` and is skipped, so what
# remains is exactly the set of unqualified calls. Each is then resolved on the
# live search path: anything owned by an attached NON-base package is a miss.
#
# This is the check that makes "namespace everything" verifiable. A pattern
# search cannot do it -- it would have to know every symbol every attached
# package exports.
#
# Note what "base" means below: OWNERSHIP ON THE LIVE SEARCH PATH, not which
# package originally defined the name. That distinction is the point. On this
# path `sd` and `match` are owned by posterior, `intersect` and `setdiff` by
# lubridate, and `setequal` by dplyr, so those calls are required to be
# qualified even though they look like base R -- including the `sd()` that
# every autoscaled prior scale is computed with. The convention of leaving base
# and stats bare rests on their origin being obvious, and a masked name's is
# not.

BASE_PKGS <- c("base", "stats", "utils", "graphics", "grDevices", "methods",
               "datasets")

# Smooth constructors inside a model formula MUST stay unqualified: brms
# detects smooth terms by name when it parses the formula, and `mgcv::t2(...)`
# is not the name it looks for. Asserted below to occur only in the registry.
FORMULA_EXEMPT <- c("t2", "s")

collect <- function(e, heads, binds) {
  if (is.call(e)) {
    h <- e[[1]]
    if (is.name(h)) {
      nm <- as.character(h)
      if (nm %in% c("::", ":::")) return(invisible(NULL))   # already qualified
      assign(nm, TRUE, envir = heads)
      # A function can also be PASSED, not called: `apply(M, 2, sd)`,
      # `vapply(x, f, ...)`, `lapply(pkgs, library, ...)`. Those are symbols in
      # argument position, never call heads, so a head-only walk cannot see
      # them -- and `apply(..., sd)` is exactly the masked-name case above.
      # Bare symbol arguments to the apply family are collected as heads too.
      if (nm %in% c("apply", "lapply", "sapply", "vapply", "mapply", "Map",
                    "do.call", "Reduce", "Filter", "purrr::map")) {
        for (j in seq_along(e)[-1]) {
          arg_is_sym <- tryCatch(is.name(e[[j]]) &&
                                   !identical(e[[j]], quote(expr = )),
                                 error = function(err) FALSE)
          if (arg_is_sym) assign(as.character(e[[j]]), TRUE, envir = heads)
        }
      }
      # local bindings: assignment targets and function formals
      if (nm %in% c("<-", "=", "<<-") && length(e) >= 2 && is.name(e[[2]]))
        assign(as.character(e[[2]]), TRUE, envir = binds)
      if (nm == "function" && length(e) >= 2 && !is.null(names(e[[2]])))
        for (p in names(e[[2]])) assign(p, TRUE, envir = binds)
    }
    # Recurse without ever BINDING an element to a variable. A call can carry
    # R's empty-argument sentinel (`M[, cf]`, a formal with no default); binding
    # that to a name makes the name test as missing, and the next use of it
    # errors. Passing `e[[i]]` straight through never performs a symbol lookup,
    # so the sentinel stays inert -- and is skipped anyway.
    for (i in seq_along(e)) {
      skip <- tryCatch(identical(e[[i]], quote(expr = )), error = function(err) TRUE)
      if (!skip) collect(e[[i]], heads, binds)
    }
  }
  invisible(NULL)
}

heads <- new.env(parent = emptyenv())
binds <- new.env(parent = emptyenv())
for (f in SWEPT_FILES) for (e in parse(f)) collect(e, heads, binds)

owner_of <- function(nm) {
  w <- tryCatch(utils::find(nm, mode = "function"), error = function(e) character(0))
  if (!length(w)) return(NA_character_)
  sub("^package:", "", w[[1]])
}

is_operator <- function(nm) grepl("^(%.*%|[^A-Za-z.]|\\.\\.\\.)", nm)

unqualified <- sort(ls(heads))
local_names  <- ls(binds)
misses <- character(0)
unresolved <- character(0)
for (nm in unqualified) {
  if (is_operator(nm)) next
  if (nm %in% local_names) next          # defined in the swept source itself
  if (nm %in% FORMULA_EXEMPT) next
  own <- owner_of(nm)
  if (is.na(own)) { unresolved <- c(unresolved, nm); next }
  if (own %in% c(".GlobalEnv", BASE_PKGS)) next
  misses <- c(misses, sprintf("%s  (%s)", nm, own))
}

check(sprintf("no unqualified call into an attached package (%d call heads examined)",
              length(unqualified)),
      length(misses) == 0, misses)
check("every unqualified call head resolves to something",
      length(unresolved) == 0, unresolved)

# The mirror check: every QUALIFIED call must name something that exists.
#
# Namespacing tends to be applied in bulk, and a bulk rewrite can invent a
# function as easily as it can find one -- `ggplot2::plot()` is the shape of it,
# from a pattern that assumes a plotting call belongs to ggplot2 in a script
# that draws with base graphics. Nothing catches that until the script runs,
# which for a supplementary figure can be weeks. The check above cannot see it
# either: a qualified call is exactly what it skips.

qualified <- new.env(parent = emptyenv())
collect_qualified <- function(e) {
  if (is.call(e)) {
    if (is.name(e[[1]]) && as.character(e[[1]]) %in% c("::", ":::") && length(e) == 3)
      assign(paste0(as.character(e[[2]]), "::", as.character(e[[3]])), TRUE, envir = qualified)
    for (i in seq_along(e)) {
      skip <- tryCatch(identical(e[[i]], quote(expr = )), error = function(err) TRUE)
      if (!skip) collect_qualified(e[[i]])
    }
  }
  invisible(NULL)
}
for (f in SWEPT_FILES) for (e in parse(f)) collect_qualified(e)

# Resolved through the package's EXPORTS, not its namespace frame. A re-export
# -- brms::loo from loo, brms::pp_check from bayesplot, dplyr::any_of from
# tidyselect, ggplot2::unit from grid -- is a perfectly good `pkg::name` that
# does not appear in `pkg`'s own frame, so a frame test flags a handful of
# working calls and buries any real one among them.
dangling <- character(0)
for (q in ls(qualified)) {
  parts <- strsplit(q, "::", fixed = TRUE)[[1]]
  found <- tryCatch({ getExportedValue(parts[1], parts[2]); TRUE },
                    error = function(e) FALSE)
  if (!isTRUE(found)) dangling <- c(dangling, q)
}
check(sprintf("every pkg::name reference resolves (%d checked)", length(ls(qualified))),
      length(dangling) == 0, dangling)

# The formula exemption must not spread beyond the registry, the only file that
# writes a model formula. Anywhere else, a bare `s()` or `t2()` is an
# unqualified call into mgcv.
smooth_hits <- unlist(lapply(names(code), function(f) {
  i <- grep("(?<![A-Za-z._])(t2|s)\\s*\\(\\s*age_years", code[[f]], perl = TRUE)
  if (length(i)) basename(f)
}))
check("bare smooth constructors appear only in specifications.R",
      length(base::setdiff(unique(smooth_hits), "specifications.R")) == 0,
      base::setdiff(unique(smooth_hits), "specifications.R"))


section("Registry invariants")

PA_KEYS  <- setdiff(names(model_templates), URB_KEYS)
EXPECTED_FAMILY <- c(
  "sos-steps" = "student", "sos-enmo" = "student",
  "ctx-steps" = "gamma",   "ctx-enmo" = "gamma",
  "osteo-steps" = "lognormal", "osteo-enmo" = "lognormal",
  "sos-steps-conf" = "student", "sos-enmo-conf" = "student",
  "ctx-steps-conf" = "gamma",   "ctx-enmo-conf" = "gamma",
  "osteo-steps-conf" = "lognormal", "osteo-enmo-conf" = "lognormal",
  "sos-urb" = "gaussian", "steps-urb" = "gaussian", "enmo-urb" = "gaussian"
)

check("15 spec keys match the expected-family table",
      setequal(names(model_templates), names(EXPECTED_FAMILY)),
      setdiff(union(names(model_templates), names(EXPECTED_FAMILY)),
              intersect(names(model_templates), names(EXPECTED_FAMILY))))

bad <- character(0)
for (k in names(model_templates)) {
  sp <- model_templates[[k]]
  fam <- sp$bf$family$family
  if (!identical(fam, unname(EXPECTED_FAMILY[[k]])))
    bad <- c(bad, sprintf("%s: family %s, expected %s", k, fam, EXPECTED_FAMILY[[k]]))
  # EVERY variable the formula names, not just the outcome and exposure. A
  # covariate that `janitor::clean_names()` spelled differently would otherwise
  # surface only as a brms error part-way into a long run.
  miss <- base::setdiff(all.vars(sp$bf$formula), names(dat))
  if (length(miss))
    bad <- c(bad, sprintf("%s: formula names columns absent from dat: %s",
                          k, paste(miss, collapse = ", ")))
  if (is.null(sp$contrast_units) || is.null(sp$contrast_label))
    bad <- c(bad, sprintf("%s: missing contrast_units / contrast_label", k))
  if (is.null(sp$exposure) || is.null(EXPOSURE_AXIS_SCALE[[sp$exposure]]))
    bad <- c(bad, sprintf("%s: exposure '%s' unregistered in EXPOSURE_AXIS_SCALE", k, sp$exposure))
}
check("families, columns, contrasts and axis scales all resolve",
      length(bad) == 0, bad)

# The industrialization specs and only they carry the community-level settings.
urb_only <- character(0)
for (k in PA_KEYS) {
  sp <- model_templates[[k]]
  if (!identical(sp$prior_kappa, KAPPA_PA))
    urb_only <- c(urb_only, sprintf("%s: prior_kappa %s, expected %s", k, sp$prior_kappa, KAPPA_PA))
  for (fld in c("prior_sd_community", "grid_cluster", "grid_n", "re_formula"))
    if (!is.null(sp[[fld]]))
      urb_only <- c(urb_only, sprintf("%s: has %s, which belongs only to the urb specs", k, fld))
  if (!identical(sp$contrast_units, 1))
    urb_only <- c(urb_only, sprintf("%s: contrast_units %s, expected 1", k, sp$contrast_units))
}
for (k in URB_KEYS) {
  sp <- model_templates[[k]]
  if (!identical(sp$prior_kappa, KAPPA_URB))
    urb_only <- c(urb_only, sprintf("%s: prior_kappa %s, expected %s", k, sp$prior_kappa, KAPPA_URB))
  if (!identical(sp$prior_sd_community, SD_COMMUNITY_URB))
    urb_only <- c(urb_only, sprintf("%s: prior_sd_community %s, expected %s", k, sp$prior_sd_community, SD_COMMUNITY_URB))
  if (!identical(sp$grid_cluster, "village_id")) urb_only <- c(urb_only, sprintf("%s: grid_cluster", k))
  if (!identical(sp$grid_n, 41))                 urb_only <- c(urb_only, sprintf("%s: grid_n %s, expected 41", k, sp$grid_n))
  if (!identical(sp$contrast_units, 10))         urb_only <- c(urb_only, sprintf("%s: contrast_units %s, expected 10", k, sp$contrast_units))
  if (!(is.logical(sp$re_formula) && is.na(sp$re_formula)))
    urb_only <- c(urb_only, sprintf("%s: re_formula must be NA (marginal of community)", k))
}
check("kappa 0.25 / 0.5 split, and community-level fields confined to the urb specs",
      length(urb_only) == 0, urb_only)


section("Postestimation grid")

grid_bad <- character(0)
for (k in names(model_templates)) {
  sp <- model_templates[[k]]
  g  <- spec_grid(sp, prep_local_data(sp, get(sp$data)))
  want <- if (k %in% URB_KEYS) 41L else 51L
  if (length(g) != want)
    grid_bad <- c(grid_bad, sprintf("%s: %d grid points, expected %d", k, length(g), want))
}
check("51 grid points for every PA spec, 41 for every urb spec",
      length(grid_bad) == 0, grid_bad)

# The community grid must come from the DISTINCT community index values. Taking
# quantiles over individuals weights by community size, which shifts the grid's
# support and with it both ends of the reported interval -- by enough, on the
# activity outcomes, to change whether it excludes zero.
sp   <- model_templates[["steps-urb"]]
d    <- prep_local_data(sp, get(sp$data))
gc_  <- spec_grid(sp, d)
xv   <- as.numeric(tapply(d[[sp$exposure]], d[[sp$grid_cluster]], function(v) v[1]))
want <- quantile(xv, sp$grid_quantiles, na.rm = TRUE)
check("urb grid endpoints are quantiles of the distinct community values",
      isTRUE(all.equal(unname(range(gc_)), unname(want))),
      sprintf("grid [%g, %g] vs community quantiles [%g, %g]",
              min(gc_), max(gc_), want[1], want[2]))
check("individual-weighted quantiles would differ (so the distinction is live)",
      !isTRUE(all.equal(unname(want),
                        unname(quantile(d[[sp$exposure]], sp$grid_quantiles, na.rm = TRUE)))))


section("Prior construction")

# One spec per family, checked against a hand computation rather than against
# build_priors()'s own arithmetic.
prior_bad <- character(0)
FLAT_EXPECTED <- c("sos-steps" = 23L, "ctx-steps" = 8L, "osteo-steps" = 8L,
                   "sos-urb" = 0L, "steps-urb" = 0L, "enmo-urb" = 0L)
for (k in names(FLAT_EXPECTED)) {
  sp  <- model_templates[[k]]
  dl  <- prep_local_data(sp, get(sp$data))
  pri <- build_priors(sp, dl)
  n_flat <- attr(pri, "n_flat")
  if (!identical(as.integer(n_flat), unname(FLAT_EXPECTED[[k]])))
    prior_bad <- c(prior_bad, sprintf("%s: %d flat community indicators, expected %d",
                                      k, n_flat, FLAT_EXPECTED[[k]]))

  fam <- sp$bf$family$family
  cls <- pri$class
  if (identical(fam, "student") != ("nu" %in% cls))
    prior_bad <- c(prior_bad, sprintf("%s: nu prior present == %s, family %s", k, "nu" %in% cls, fam))
  if (identical(fam, "gamma") != ("shape" %in% cls))
    prior_bad <- c(prior_bad, sprintf("%s: shape prior present == %s, family %s", k, "shape" %in% cls, fam))
  if ((k %in% URB_KEYS) != ("sd" %in% cls))
    prior_bad <- c(prior_bad, sprintf("%s: community-SD prior present == %s", k, "sd" %in% cls))

  # community-SD prior: exponential(1 / (0.33 * sd_lp))
  if (k %in% URB_KEYS) {
    sd_lp <- attr(pri, "sd_lp")
    want  <- sprintf("exponential(%.8g)", 1 / (SD_COMMUNITY_URB * sd_lp))
    got   <- pri$prior[cls == "sd"][1]
    if (!identical(got, want))
      prior_bad <- c(prior_bad, sprintf("%s: sd prior %s, expected %s", k, got, want))
  }

  # one coefficient's normal(0, s), recomputed by hand from the design matrix
  b_rows <- which(cls == "b" & nzchar(pri$coef))
  if (length(b_rows)) {
    cf    <- pri$coef[b_rows[1]]
    sdat  <- brms::standata(sp$bf, data = dl)
    M     <- do.call(cbind, lapply(c("X", "Xs"), function(m)
               if (is.null(sdat[[m]])) NULL else as.matrix(sdat[[m]])))
    sx    <- stats::sd(M[, cf])
    sd_lp <- attr(pri, "sd_lp")
    want  <- sprintf("normal(0, %.8g)", sp$prior_kappa * sd_lp / sx)
    got   <- pri$prior[b_rows[1]]
    if (!identical(got, want))
      prior_bad <- c(prior_bad, sprintf("%s / %s: prior %s, hand-computed %s", k, cf, got, want))
  }
}
check("flat-indicator counts, family-specific classes, community SD and one autoscaled coefficient",
      length(prior_bad) == 0, prior_bad)


section("Interval primitives")

set.seed(11)
ng <- 21L; nd <- 2000L
gr <- seq(0, 2, length.out = ng)
E  <- matrix(rnorm(nd), nd, 1) %*% matrix(1, 1, ng)
Ym <- exp(0.5 * E + matrix(rep(gr, each = nd), nd, ng) * 0.3 + rnorm(nd * ng, 0, 0.15))
dd <- tibble::tibble(drawid = rep(seq_len(nd), times = ng),
                     xx     = rep(gr, each = nd),
                     draw   = as.vector(Ym))

lp <- linear_projection(dd, exposure = xx, level = 0.95)
check("linear_projection returns exactly one interval row", nrow(lp$beta_hpdi) == 1L,
      paste("got", nrow(lp$beta_hpdi), "rows"))
ref <- ggdist::median_hdci(lp$beta_draws, beta, .width = 0.95)
check("its bounds equal ggdist::median_hdci recomputed on the same draws",
      isTRUE(all.equal(c(lp$beta_hpdi$lo, lp$beta_hpdi$hi),
                       c(ref$.lower[1], ref$.upper[1]))),
      sprintf("[%.10g, %.10g] vs [%.10g, %.10g]",
              lp$beta_hpdi$lo, lp$beta_hpdi$hi, ref$.lower[1], ref$.upper[1]))

b_h <- simul_credible_bands(dd, exposure = xx, levels = 0.95, interval_type = "HPDI")
b_s <- simul_credible_bands(dd, exposure = xx, levels = 0.95, interval_type = "symmetric")
check("interval_type is a live conditional, not a label",
      !isTRUE(all.equal(c(b_h$lo, b_h$hi), c(b_s$lo, b_s$hi))))
check("both constructions carry their own name in the output",
      identical(unique(b_h$interval_type), "HPDI") &&
      identical(unique(b_s$interval_type), "symmetric"))

realised <- function(b) {
  lo <- b$lo[base::match(gr, b$xx)]; hi <- b$hi[base::match(gr, b$xx)]
  mean(rowSums(Ym >= rep(lo, each = nd) & Ym <= rep(hi, each = nd)) == ng)
}
check(sprintf("HPDI band achieves its stated simultaneous coverage (%.3f)", realised(b_h)),
      abs(realised(b_h) - 0.95) < 0.02)
check(sprintf("symmetric band achieves its stated simultaneous coverage (%.3f)", realised(b_s)),
      abs(realised(b_s) - 0.95) < 0.02)

# On a skewed posterior the HPD band must sit asymmetrically about the median;
# the symmetric one cannot, however it is calibrated. This is the substantive
# difference between the two, not just a different number.
medv <- apply(Ym, 2, median)
asym <- function(b) mean((b$hi - medv) / (medv - b$lo))
check(sprintf("HPDI band is asymmetric about the median (ratio %.2f vs %.2f symmetric)",
              asym(b_h), asym(b_s)),
      asym(b_h) > 1.5 * asym(b_s))

bm <- simul_credible_bands(dd, exposure = xx,
                           levels = c(0.05, 0.25, 0.50, 0.75, 0.95),
                           interval_type = "HPDI")
wid <- tapply(bm$hi - bm$lo, bm$level, mean)
check("nested levels give strictly nested bands", all(diff(wid) > 0),
      paste(names(wid), round(wid, 4), collapse = "; "))

old <- getOption("bone.band_type")
options(bone.band_type = "symmetric")
by_option <- unique(simul_credible_bands(dd, exposure = xx, levels = 0.95)$interval_type)
options(bone.band_type = old)
check("the default follows options(bone.band_type), so BAND_TYPE reaches every figure",
      identical(by_option, "symmetric"), by_option)


section("Contrast units reach the AMEF panel")

# The AMEF y-axis must be on the reported contrast. Built for a urb spec
# (contrast_units = 10) and a PA spec (1), and read back off the layer data.
band_for <- function(k) {
  sp <- model_templates[[k]]
  b  <- dd
  names(b)[names(b) == "xx"] <- sp$exposure
  list(spec = sp, bands = simul_credible_bands(b, exposure = !!rlang::sym(sp$exposure),
                                               levels = 0.95, function_type = "AMEF"),
       lp = linear_projection(b, exposure = !!rlang::sym(sp$exposure), level = 0.95))
}
cu_bad <- character(0)
for (k in c("steps-urb", "sos-steps")) {
  z  <- band_for(k)
  p  <- make_amef_panel(z$spec, z$bands, z$lp)
  yl <- p$labels$y
  drawn <- ggplot2::layer_data(p, 1)          # first nested ribbon (level 0.95)
  ratio <- max(drawn$ymax, na.rm = TRUE) / max(z$bands$hi, na.rm = TRUE)
  if (!isTRUE(all.equal(ratio, z$spec$contrast_units)))
    cu_bad <- c(cu_bad, sprintf("%s: panel/band ratio %.4f, contrast_units %g",
                                k, ratio, z$spec$contrast_units))
  if (!grepl(z$spec$contrast_label, yl, fixed = TRUE))
    cu_bad <- c(cu_bad, sprintf("%s: y-axis title '%s' omits '%s'",
                                k, gsub("\n", " / ", yl), z$spec$contrast_label))
}
check("AMEF bands scaled by contrast_units, and the label names the denominator",
      length(cu_bad) == 0, cu_bad)


section("Rug reflects the analytic sample")

# The rug must show the rows the fit used, not the whole cohort. On the
# biomarker specs the difference is roughly threefold.
sp   <- model_templates[["ctx-steps"]]
n_an <- nrow(prep_local_data(sp, get(sp$data)))
n_co <- sum(!is.na(get(sp$data)[[sp$exposure]]))
check(sprintf("analytic n (%d) is well below the cohort's exposure-complete n (%d)",
              n_an, n_co), n_an < 0.5 * n_co)
check("make_amef_panel's rug_data default is the analytic frame",
      identical(deparse(formals(make_amef_panel)$rug_data), "NULL") &&
        any(grepl("prep_local_data", deparse(body(make_amef_panel)), fixed = TRUE)))


# ---- Final figures: every PNG at least as new as its PDF ----
#
# `just figures-png` converts whatever PDFs are on disk; it has no way of
# knowing whether they are the current ones. Run it while a figure build is
# still writing and it exports the superseded image over the new one, with
# nothing raised. The PNG is the copy that gets inserted into a document, so a
# stale one is the version a reader ends up seeing.
#
# Pairs only. A PDF with no PNG is a skip, not a failure -- the PNG is a
# convenience export and not every figure has one. A missing outputs/ tree is a
# skip too, so a fresh clone with no figures built yet still passes.

section("Final figures")

fig_dir <- here::here("outputs", "figures", "final")
fig_pdfs <- if (dir.exists(fig_dir)) {
  list.files(fig_dir, pattern = "\\.pdf$", full.names = TRUE)
} else {
  character(0)
}

if (length(fig_pdfs) == 0L) {
  cat("  --    no built figures to check\n")
} else {
  stale <- vapply(fig_pdfs, function(p) {
    png <- sub("\\.pdf$", ".png", p)
    # A missing PNG is a skip; only one that exists and predates its PDF means
    # someone is looking at a superseded image.
    file.exists(png) && file.mtime(png) < file.mtime(p)
  }, logical(1))

  check(sprintf("all %d final PNG(s) are at least as new as their PDF",
                sum(file.exists(sub("\\.pdf$", ".png", fig_pdfs)))),
        !any(stale),
        if (any(stale)) {
          c(paste("stale:", paste(basename(fig_pdfs[stale]), collapse = ", ")),
            "run `just figures-png`, and only after every figure build has finished")
        })
}


# ---- Optional: one smoke fit end to end ----
#
# The only way to confirm that `re_formula = NA` survives the trip through
# marginaleffects into brms. Off by default because it compiles a Stan model.

if (nzchar(Sys.getenv("CHECK_FIT"))) {
  section("Smoke fit (CHECK_FIT set)")
  key <- "steps-urb"
  sp  <- model_templates[[key]]
  fit <- fit_one_spec(key, dat, mode = "smoke")
  check("smoke fit returns a brmsfit", inherits(fit, "brmsfit"))

  pd <- aerf_draws(fit, sp, dat)
  check("aerf_draws accepts re_formula = NA and returns the full grid",
        length(unique(pd[[sp$exposure]])) == 41L,
        paste("got", length(unique(pd[[sp$exposure]])), "grid points"))

  # re_formula must actually change something, or it is not being forwarded.
  sp_null <- sp; sp_null$re_formula <- NULL
  pd_null <- aerf_draws(fit, sp_null, dat)
  sd_na   <- stats::sd(pd$draw);  sd_null <- stats::sd(pd_null$draw)
  check(sprintf("re_formula = NA reaches brms and widens the posterior (sd %.4g vs %.4g)",
                sd_na, sd_null),
        sd_na > sd_null)

  sl <- amef_draws(fit, sp, dat)
  check("amef_draws returns draws on the same grid",
        setequal(unique(sl[[sp$exposure]]), unique(pd[[sp$exposure]])))
} else {
  cat("\n(skipping the smoke fit; set CHECK_FIT=1 to include it)\n")
}


# ---- Verdict ----

cat("\n")
if (FAILURES == 0L) {
  cat("PASS -- the startup set is consistent with the specification.\n")
} else {
  cat(sprintf("FAIL -- %d check(s) failed. Do not start the long runs.\n", FAILURES))
  quit(status = 1L)
}
