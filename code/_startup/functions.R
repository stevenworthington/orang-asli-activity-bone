###############################################################################
# Shared function definitions. Sourced by `_startup/init.R`, so every analysis
# script gets these for free. This file DEFINES and does not EXECUTE: sourcing
# it has no side effects beyond binding the functions below.
#
# What lives here, in the order it appears:
#   load_pkgs                  package-loading wrapper (no install)
#   theme_pub / theme_pub_leg  publication ggplot themes
#   skew                       Fisher-Pearson skewness, for the PP checks
#   prep_dat                   THE data prep: read the CSV, derive the analytic
#                              columns. Shared with `_targets.R`.
#   build_priors               THE prior scheme: autoscaled at fit time
#   calc_probs_direction       directional probabilities + median + HDI
#   compute_curvature_draws    finite-difference second derivative
#   linear_projection          per-draw best linear approximation of the AERF
#   nested_simul_ribbons       the five nested grey band layers
#   format_thousands           compact axis-tick formatter
#   x_scale_for                x-axis ticks in natural exposure units
#   make_aerf_panel            AERF panel (manuscript convention)
#   make_amef_panel            AMEF panel; puts the AMEF on the REPORTED
#                              contrast scale -- see its header
#   make_age_cond_amef_panel_with_bands
#   stack_subfig               AERF over AMEF, one subfigure
#   aerf_amef_plot             the 3-panel working diagnostic figure
#   pp_check_stats             six-statistic posterior-predictive check
#   simul_credible_bands       simultaneous credible bands, two constructions
#
# Package calls are namespaced explicitly. Base R and `stats` are NOT -- the
# convention earns its keep by disambiguating functions whose origin is
# unobvious, and `stats::median()` everywhere only makes the source harder to
# read. Operators and S3 methods (`+`, `|>`, `%in%`) stay bare, which is why
# `packages.R` still attaches the core set.
#
# WITH ONE IMPORTANT QUALIFICATION: a base or stats name that an ATTACHED
# package masks is not obvious, and gets qualified. On this search path `sd` and
# `match` resolve to posterior, `intersect` and `setdiff` to lubridate, and
# `setequal` to dplyr -- so an unqualified `sd(y)` reaches `posterior::sd`,
# including in `sd_lp`, the number every autoscaled prior is built from. Each
# masking method happens to dispatch back to the base implementation for the
# arguments used here, so no result depends on it; the point is that the source
# should say what it calls. The AST check in `code/_checks/verify-startup.R`
# enforces this, and will catch the next package that masks something.
#
# ONE DELIBERATE EXCEPTION: smooth constructors inside a model formula
# (`t2()`, `s()` in `specifications.R`) must stay unqualified. brms detects
# smooth terms by name when it parses the formula, and `mgcv::t2(...)` is not
# the name it looks for.
###############################################################################


# ---- Package loading wrapper ----
#
# Strictly LOADS packages (no install). If a package isn't in the project
# library, library() errors and the fix is `renv::restore()`.
load_pkgs <- function(...) {
  pkgs <- c(...)
  invisible(lapply(pkgs, library, character.only = TRUE))
}


# ---- ggplot themes ----

theme_pub <- function(base_size = 10) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.3, color = "grey90"),
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      plot.title       = ggplot2::element_text(size = base_size + 1, face = "bold"),
      axis.title       = ggplot2::element_text(size = base_size),
      legend.position  = "none"
    )
}

theme_pub_leg <- function(base_size = 10) {
  theme_pub(base_size) +
    ggplot2::theme(legend.position = "top")
}


# ---- Shape ----

# Fisher-Pearson skewness (NIST definition). Used as a PP-check statistic.
skew <- function(y) {
  y <- na.omit(y)
  n <- length(y)
  diff <- y - mean(y)
  (sqrt(n - 1) / (n - 2)) * n * (sum(diff^3) / (sum(diff^2)^1.5))
}


# ---- prep_dat: read + derive the analytic columns ----
#
# THE single definition of the project's data prep. `_startup/data.R` calls it
# to build `dat`; `_targets.R` calls it inside the `dat` target so the pipeline
# gets a file-format dependency on the CSV. One definition rather than two, so
# a scaling can only be defined in one place.
#
# Scaling rationale lives in `_startup/data.R`'s header.

prep_dat <- function(path) {
  readr::read_csv(path, show_col_types = FALSE) |>
    janitor::clean_names(case = "snake") |>
    dplyr::mutate(
      sex                   = factor(sex, levels = c("female", "male")),
      village_id            = factor(village_id),
      # One scaling per variable, used in every role it plays.
      tibia_sos_200         = tibia_sos / 200,                   # outcome (PA, urb)
      ad_steps_5k           = ad_tot_step_count_0_24hr / 5000,   # exposure (PA), outcome (urb)
      enmo_10               = ad_mean_enmo_mg_0_24hr / 10,       # exposure (PA), outcome (urb)
      osteocalcin_pg_ml_10k = osteocalcin_pg_ml / 10000,         # outcome (PA)
      # Body composition (z-scored), for the confounder-DAG variants. Fat mass
      # and lean body mass are one bundled DAG node but two regression
      # covariates.
      fat_mass_kg_z         = as.numeric(scale(fat_mass_kg)),
      fat_free_mass_kg_z    = as.numeric(scale(fat_free_mass_kg))
    )
}


# ---- build_priors: autoscaled weakly informative priors for one spec ----
#
# Priors are autoscaled to the ANALYTIC data, so they cannot be precomputed
# when the registry is defined -- hence a function called at fit time.
#
# The rule, rstanarm's autoscaling with student_t swapped for normal:
#
#   b          normal(0, kappa * sd_lp / sd(x_j))    per coefficient
#   Intercept  normal(mean_lp, 0.5 * sd_lp)
#   sds        exponential(2)                        smooth wiggliness
#   sigma      exponential(1 / sd_lp)                Gaussian / Student-t / lognormal
#   shape      exponential(0.1)                      Gamma only
#   nu         gamma(2, 0.1)                         Student-t only
#   sd         exponential(1 / (sd_community * sd_lp))   community RI, urb only
#
# `sd_lp` is the outcome SD on the LINEAR-PREDICTOR scale: log(y) under a log
# link, y otherwise. Autoscaling on the response scale under a log link would
# set every prior from a spread the linear predictor never sees, so the link is
# consulted here rather than assumed.
#
# `sds` is NOT autoscaled -- brms offers no autoscaling for it -- which is the
# whole reason data.R rescales the identity-link outcomes to sd ~ 1. It cannot
# do the same for the log-link outcomes: dividing y by a constant shifts the
# linear predictor without changing sd(log y). See data.R.
#
# COMMUNITY INDICATORS ARE LEFT FLAT in the PA specs, deliberately. The
# treatment-coded `village_id` indicators (23 on the SOS specs, 8 on the
# biomarker specs, which reach fewer communities) are nuisance parameters whose
# job is to absorb between-community variation completely, so the activity
# effect is identified strictly from within-community contrasts -- which is what
# community fixed effects are for. A finite prior shrinks them toward the
# reference category, letting between-community information back into the
# exposure estimate and making the answer depend on which community happens to
# be the reference. The industrialization specs have no indicators; they carry a
# community random intercept instead, whose SD prior IS scaled.
#
# Returns a brmsprior carrying `sd_lp` and `n_flat` attributes. `spec` must
# supply `outcome`, `bf`, `prior_kappa`, and (for the community-level specs)
# `prior_sd_community`.

build_priors <- function(spec, d) {
  fam    <- spec$bf$family$family
  is_log <- fam %in% c("lognormal", "gamma") ||
              identical(spec$bf$family$link, "log")
  y  <- d[[spec$outcome]]
  if (is_log && !all(y > 0, na.rm = TRUE))
    stop("non-positive outcome values under a log link: ", spec$outcome)
  lp    <- if (is_log) log(y) else y
  sd_lp <- stats::sd(lp, na.rm = TRUE)
  mu_lp <- mean(lp, na.rm = TRUE)

  # Per-coefficient SDs of the design matrix, so each prior is on that
  # coefficient's own scale. Xs carries the smooth null-space columns, where
  # the exposure's linear component -- the parameter the reported slope is
  # made of -- actually lives.
  sdat  <- brms::standata(spec$bf, data = d)
  colsd <- c()
  for (m in c("X", "Xs")) if (!is.null(sdat[[m]])) {
    M <- as.matrix(sdat[[m]])
    v <- apply(M, 2, stats::sd); names(v) <- colnames(M)
    colsd <- c(colsd, v)
  }
  colsd <- colsd[is.finite(colsd) & colsd > 0]

  gp    <- as.data.frame(brms::get_prior(spec$bf, data = d))
  coefs <- gp$coef[gp$class == "b" & gp$coef != ""]
  flat  <- grep("^village_id", coefs, value = TRUE)   # see header
  pri   <- NULL
  scaled <- base::setdiff(coefs, flat)
  for (cf in scaled) {
    # `[` not `[[`: a name absent from colsd makes `[[` error before any guard
    # can run, and the guard here needs to distinguish "no design column" from
    # "column with zero variance" -- both of which must be reported, not skipped.
    sx <- unname(colsd[cf])
    if (is.na(sx) || !is.finite(sx))
      stop("no usable design-matrix scale for coefficient '", cf, "' in ",
           spec$outcome, "'s model. build_priors() cannot autoscale it, and ",
           "leaving it flat would be a silent modelling change.")
    p <- brms::set_prior(sprintf("normal(0, %.8g)", spec$prior_kappa * sd_lp / sx),
                         class = "b", coef = cf)
    pri <- if (is.null(pri)) p else pri + p
  }

  # Every coefficient that should be scaled has been. The community indicators
  # are the only ones deliberately left without a prior; anything else slipping
  # through would be an unintended flat prior on a real covariate.
  got <- pri$coef[pri$class == "b" & nzchar(pri$coef)]
  if (!base::setequal(got, scaled))
    stop("build_priors(): coefficients without a prior: ",
         paste(base::setdiff(scaled, got), collapse = ", "))

  pri <- pri +
    brms::set_prior(sprintf("normal(%.8g, %.8g)", mu_lp, 0.5 * sd_lp), class = "Intercept") +
    brms::set_prior("exponential(2)", class = "sds")

  pri <- if (identical(fam, "gamma"))
    pri + brms::set_prior("exponential(0.1)", class = "shape")
  else
    pri + brms::set_prior(sprintf("exponential(%.8g)", 1 / sd_lp), class = "sigma")

  if (identical(fam, "student")) pri <- pri + brms::set_prior("gamma(2, 0.1)", class = "nu")

  if (!is.null(spec$prior_sd_community))
    pri <- pri + brms::set_prior(
      sprintf("exponential(%.8g)", 1 / (spec$prior_sd_community * sd_lp)), class = "sd")

  attr(pri, "sd_lp")   <- sd_lp
  attr(pri, "n_flat")  <- length(flat)
  pri
}


# ---- Posterior probability summaries ----

# Directional probabilities + median + 95% HPDI for a single draw vector.
calc_probs_direction <- function(draws) {
  hdi <- ggdist::hdci(draws, .width = 0.95)   # 1 x 2 matrix: lower, upper
  tibble::tibble(
    prob_positive = mean(draws > 0),
    prob_negative = mean(draws < 0),
    median        = median(draws),
    hdi_lower     = unname(hdi[1]),
    hdi_upper     = unname(hdi[2])
  )
}


# ---- Curvature draws (finite-difference second derivative) ----
#
# Kept in the helper set even though curvature is not part of the manuscript
# figure convention -- the reported figures show AERF + AMEF only, and the
# linear-projection slope carries the headline. The pipeline target
# `curvature_draws_<spec>` continues to produce these; the only quantity read
# from them downstream is a SIGN probability, which no positive rescaling can
# change.

compute_curvature_draws <- function(slope_draws, exposure) {
  exp_sym <- rlang::sym(exposure)
  slope_draws |>
    dplyr::arrange(drawid, !!exp_sym) |>
    dplyr::group_by(drawid) |>
    dplyr::mutate(
      x_lag     = dplyr::lag(!!exp_sym),
      x_lead    = dplyr::lead(!!exp_sym),
      draw_lag  = dplyr::lag(draw),
      draw_lead = dplyr::lead(draw),
      curvature = (draw_lead - draw_lag) / (x_lead - x_lag)
    ) |>
    dplyr::filter(!is.na(curvature)) |>
    dplyr::ungroup() |>
    dplyr::transmute(drawid, !!exp_sym := !!exp_sym, draw = curvature)
}


# ---- linear_projection: per-draw best linear approximation of the AERF ----
#
# For each posterior draw of the AERF, fit y = intercept + beta * x and
# return the slope summary plus a band on the fitted line over the exposure
# grid. The band comes from simul_credible_bands() at its default construction,
# so it reads as "joint coverage of the entire line at `level`" and matches the
# grey bands drawn beside it whichever construction is in force.
#
# UNITS. `beta` is per ONE unit of the modelled exposure column, always.
# Multiply by `spec$contrast_units` to reach the reported contrast. The
# scaling is deliberately NOT done here: the table scripts apply it at their
# own call sites, and doing it in both places would square it.
#
# Returns a list with:
#   $beta_draws : per-draw slope tibble (drawid, beta)
#   $beta_hpdi  : tibble(lo, mid, hi) with the slope's `level` HPDI
#   $band       : tibble(<exposure>, lo, mid (median), hi, level, ...) — the
#                 projection band for the AERF panel
#
# Used by make_aerf_panel() (band overlay) and make_amef_panel() (slope HPDI
# rendered as red-dashed horizontals + filled rect across the exposure range).

linear_projection <- function(pred_draws, exposure, level = 0.95,
                              value_col = "draw") {

  stopifnot(is.data.frame(pred_draws))

  exposure_sym <- rlang::ensym(exposure)
  exposure_nm  <- rlang::as_name(exposure_sym)

  req_cols <- c("drawid", exposure_nm, value_col)
  if (!all(req_cols %in% names(pred_draws))) {
    stop("pred_draws must contain drawid, exposure, and value columns")
  }
  if (!is.numeric(level) || length(level) != 1 || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("level must be a single number strictly between 0 and 1.")
  }

  x_grid <- pred_draws |>
    dplyr::distinct(.data[[exposure_nm]]) |>
    dplyr::arrange(.data[[exposure_nm]]) |>
    dplyr::pull(.data[[exposure_nm]])

  # Best linear projection of the AERF per posterior draw, computed as the
  # closed-form OLS coefficients (cov(x,y)/var(x)) -- a deterministic linear
  # functional of the Bayesian posterior, NOT a model fit. Numerically identical
  # to lm() coefficients; kept lm-free so no frequentist call appears anywhere on
  # the reported-inference path. (The power and calibration harnesses do fit
  # frequentist twins, deliberately and only as fast proxies.)
  lin_summ <- pred_draws |>
    dplyr::filter(!is.na(.data[[value_col]]), !is.na(.data[[exposure_nm]])) |>
    dplyr::group_by(drawid) |>
    dplyr::summarize(
      beta = {
        xc <- .data[[exposure_nm]] - mean(.data[[exposure_nm]])
        sum(xc * (.data[[value_col]] - mean(.data[[value_col]]))) / sum(xc^2)
      },
      intercept = mean(.data[[value_col]]) - beta * mean(.data[[exposure_nm]]),
      .groups   = "drop"
    )

  # draw-level projected line values on the exposure grid
  lin_proj_draws <- lin_summ |>
    dplyr::select(drawid, intercept, beta) |>
    tidyr::crossing(tibble::tibble(!!exposure_sym := x_grid)) |>
    dplyr::mutate(draw = intercept + beta * .data[[exposure_nm]])

  # Slope summary. `median_hdci` is SAMPLE-based (the narrowest continuous
  # interval spanning the required mass of draws) rather than density-based, so
  # it carries no kernel-bandwidth dependence and cannot return a split
  # interval. It is the primitive used for every scalar HPDI in the project.
  beta_summary <- lin_summ |>
    ggdist::median_hdci(beta, .width = level) |>
    dplyr::rename(lo = .lower, mid = beta, hi = .upper) |>
    dplyr::select(lo, mid, hi)

  # simul-band on the AERF projection. `interval_type` is left at its default so
  # the red envelope is built the same way as the grey bands beside it.
  band <- simul_credible_bands(
    draws_df      = lin_proj_draws,
    exposure      = !!exposure_sym,
    value_col     = "draw",
    levels        = level,
    function_type = "AERF"
  )

  list(
    level      = level,
    beta_draws = lin_summ |> dplyr::select(drawid, beta),
    beta_hpdi  = beta_summary,
    band       = band
  )
}


# ---- Nested simul-band ribbons (5/25/50/75/95% in greys) ----
#
# Helper used by make_aerf_panel / make_amef_panel. Returns five nested
# geom_ribbon layers in increasing alpha, one per level.
#
# It takes NO arguments: each layer selects its own rows with a `~ .x` lambda,
# which ggplot2 evaluates against the PLOT's data at render time. So the plot
# must be built on the long-format tibble `simul_credible_bands(levels =
# c(0.05, 0.25, 0.50, 0.75, 0.95))` returns, one row per exposure x level.

nested_simul_ribbons <- function() {
  list(
    ggplot2::geom_ribbon(data = ~ dplyr::filter(.x, level == 0.95),
                         ggplot2::aes(ymin = lo, ymax = hi),
                         fill = "grey20", alpha = 0.10),
    ggplot2::geom_ribbon(data = ~ dplyr::filter(.x, level == 0.75),
                         ggplot2::aes(ymin = lo, ymax = hi),
                         fill = "grey20", alpha = 0.12),
    ggplot2::geom_ribbon(data = ~ dplyr::filter(.x, level == 0.50),
                         ggplot2::aes(ymin = lo, ymax = hi),
                         fill = "grey20", alpha = 0.14),
    ggplot2::geom_ribbon(data = ~ dplyr::filter(.x, level == 0.25),
                         ggplot2::aes(ymin = lo, ymax = hi),
                         fill = "grey20", alpha = 0.15),
    ggplot2::geom_ribbon(data = ~ dplyr::filter(.x, level == 0.05),
                         ggplot2::aes(ymin = lo, ymax = hi),
                         fill = "grey20", alpha = 0.15)
  )
}


# ---- Panel + stack builders (manuscript convention) ----
#
# Each "cell" is a vertical stack of an AERF panel (top) over an AMEF panel
# (bottom), sharing the exposure x-axis. The AERF panel hides its x-axis
# (carried by the AMEF below). Both panels show nested grey simul-bands + a
# red-dashed linear-projection overlay. The AMEF panel additionally shows a
# zero line and a rug of the exposure values actually used by the fit.

# Compact axis-label formatter: 12500 -> "12.5k", 5000 -> "5k", 30 -> "30",
# 0.15 -> "0.15". Dynamic decimal precision so tick labels are
# distinguishable. Used as the default y-axis labeller on AERF / AMEF panels.
format_thousands <- function(x) {
  if (!is.numeric(x)) stop("Input must be numeric")
  if (length(x) == 0) return(character(0))

  big <- x[!is.na(x) & abs(x) >= 1000]
  digits <- 1L
  if (length(big) >= 2) {
    diffs_k <- diff(sort(unique(big))) / 1000
    diffs_k <- diffs_k[diffs_k > 0]
    if (length(diffs_k) > 0) {
      digits <- max(0L, as.integer(ceiling(-log10(min(diffs_k)))))
    }
  }

  sapply(x, function(y) {
    if (is.na(y)) return(NA_character_)
    if (abs(y) < 1000) return(as.character(y))
    paste0(formatC(round(y / 1000, digits), format = "fg"), "k")
  })
}

# X-axis tick formatting. Ticks are drawn in NATURAL units: counts of 1,000 or
# more get a "k" suffix (5k, 10k, ...), everything else prints as a number.
# The per-exposure multiplier comes from EXPOSURE_AXIS_SCALE in
# `_startup/specifications.R`, and an unregistered exposure STOPS rather than
# silently falling through to raw modelling units.
x_scale_for <- function(exposure_nm) {
  # Single-bracket lookup, deliberately: `[[` on a named vector ERRORS on a
  # missing name, which would make the guard below unreachable and turn an
  # unregistered exposure into "subscript out of bounds" instead of the message
  # this function exists to give. `[` yields NA, which the guard sees.
  f <- unname(EXPOSURE_AXIS_SCALE[exposure_nm])
  if (is.na(f))
    stop("no axis scale registered for exposure '", exposure_nm,
         "'. Add it to EXPOSURE_AXIS_SCALE in _startup/specifications.R.")
  ggplot2::scale_x_continuous(labels = function(x) {
    v <- x * f
    ifelse(is.na(v), NA_character_,
           ifelse(abs(v) >= 1000,
                  paste0(formatC(v / 1000, format = "fg", digits = 3), "k"),
                  formatC(v, format = "fg", digits = 3)))
  })
}

make_aerf_panel <- function(spec, simul_bands_aerf, lin_proj,
                            y_label = NULL, x_axis_label = NULL,
                            x_limits = NULL) {
  exp_sym <- rlang::sym(spec$exposure)
  bands   <- simul_bands_aerf |> dplyr::arrange(dplyr::desc(level))
  # Two-line y-axis label so the title fits inside the cell without clipping
  # at 2.6"-wide subfigures.
  if (is.null(y_label)) y_label <- paste0("Predicted\n", spec$outcome_label)

  p <- ggplot2::ggplot(bands, ggplot2::aes(x = !!exp_sym)) +
    nested_simul_ribbons() +
    # red-dashed linear-projection envelope on the AERF
    ggplot2::geom_ribbon(data = lin_proj$band,
                         ggplot2::aes(x = !!exp_sym, ymin = lo, ymax = hi),
                         fill = "red", alpha = 0.10, inherit.aes = FALSE) +
    ggplot2::geom_line(data = lin_proj$band,
                       ggplot2::aes(x = !!exp_sym, y = lo),
                       color = "red", linewidth = 0.5, linetype = "dashed",
                       inherit.aes = FALSE) +
    ggplot2::geom_line(data = lin_proj$band,
                       ggplot2::aes(x = !!exp_sym, y = hi),
                       color = "red", linewidth = 0.5, linetype = "dashed",
                       inherit.aes = FALSE) +
    x_scale_for(spec$exposure) +
    ggplot2::scale_y_continuous(labels = format_thousands) +
    ggplot2::labs(x = x_axis_label, y = y_label) +
    theme_pub() +
    ggplot2::theme(axis.title = ggplot2::element_text(face = "plain", size = 9))

  # Apply per-column x-limit (used by shared-axes figure assembly to align
  # cells in the same exposure column). coord_cartesian rather than
  # scale_x_continuous(limits=) because the latter clips data outside the
  # range, which would drop simul-band edges; coord_cartesian only changes
  # the view.
  if (!is.null(x_limits)) p <- p + ggplot2::coord_cartesian(xlim = x_limits)
  p
}

# AMEF panel.
#
# THE REPORTED CONTRAST IS APPLIED HERE, and only here. An AMEF is a rate and
# therefore always carries a denominator; upstream it is per ONE unit of the
# modelled exposure column, which is a modelling convenience rather than a
# reporting decision. This panel multiplies both the bands and the
# linear-projection slope HPDI by `spec$contrast_units` and names the
# denominator in the y-axis title, so a reader can read the interval quoted in
# the Results straight off the figure. For the PA specs `contrast_units` is 1
# (the data.R divisors equal the reported contrasts), so nothing moves; for the
# industrialization specs it is 10, so the figures and the tables quote the same
# denominator.
#
# Doing it here rather than in `amef_draws()` keeps every SAVED artifact -- the
# targets draws, `linear_projection()`'s beta, the two table scripts -- on one
# consistent per-unit scale, with a single presentation-layer conversion.
#
# `rug_data` is the data the rug is drawn from. It defaults to the analytic
# frame for this spec, i.e. the rows the fit actually used. Drawing it from the
# full cohort instead would overstate data support on the biomarker panels
# roughly threefold -- a rug of ~860 under an analytic n of 265 -- directly
# beneath the panel whose job is to show where the estimate is informed. A
# caller working from a subset (the age-subset supplement figure) passes that
# subset.

make_amef_panel <- function(spec, simul_bands_amef, lin_proj,
                            slope_draws = NULL, y_label = NULL,
                            x_axis_label = NULL,
                            x_limits = NULL,
                            show_x_axis = TRUE,
                            rug_data = NULL) {
  exp_sym <- rlang::sym(spec$exposure)
  cu      <- spec$contrast_units
  if (is.null(cu) || !is.finite(cu) || cu <= 0)
    stop("spec$contrast_units must be a positive number; got ",
         format(cu), " for exposure '", spec$exposure, "'.")

  # Everything plotted on the y-axis moves to the reported contrast together.
  bands <- simul_bands_amef |>
    dplyr::arrange(dplyr::desc(level)) |>
    dplyr::mutate(dplyr::across(dplyr::any_of(c("m", "sd", "med", "lo", "hi")),
                                function(v) v * cu))
  beta_lo <- lin_proj$beta_hpdi$lo * cu
  beta_hi <- lin_proj$beta_hpdi$hi * cu

  if (is.null(y_label))
    y_label <- paste0("Slope of predicted\n", spec$outcome_label,
                      "\n", spec$contrast_label)
  if (is.null(x_axis_label)) x_axis_label <- spec$exposure_label

  # Linear-projection slope band: bound to the simul-band grid (1st-99th
  # percentile of the exposure), matching the AERF panel's envelope extent.
  # annotate("segment") / annotate("rect") with explicit grid endpoints rather
  # than -Inf / +Inf, so the band doesn't extend beyond the data when
  # coord_cartesian(xlim = ...) widens the plot view (shared-axes layout).
  grid_lo <- min(bands[[spec$exposure]], na.rm = TRUE)
  grid_hi <- max(bands[[spec$exposure]], na.rm = TRUE)

  p <- ggplot2::ggplot(bands, ggplot2::aes(x = !!exp_sym)) +
    nested_simul_ribbons() +
    ggplot2::annotate("rect",
                      xmin = grid_lo, xmax = grid_hi,
                      ymin = beta_lo, ymax = beta_hi,
                      fill = "red", alpha = 0.10) +
    ggplot2::annotate("segment",
                      x = grid_lo, xend = grid_hi,
                      y = beta_lo, yend = beta_lo,
                      color = "red", linewidth = 0.5, linetype = "dashed") +
    ggplot2::annotate("segment",
                      x = grid_lo, xend = grid_hi,
                      y = beta_hi, yend = beta_hi,
                      color = "red", linewidth = 0.5, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
    x_scale_for(spec$exposure) +
    ggplot2::scale_y_continuous(labels = format_thousands) +
    ggplot2::labs(x = x_axis_label, y = y_label) +
    theme_pub() +
    ggplot2::theme(axis.title = ggplot2::element_text(face = "plain", size = 9))

  # Rug at the bottom of the AMEF panel, trimmed to the prediction grid.
  if (!is.null(slope_draws)) {
    rug_lo <- min(slope_draws[[spec$exposure]], na.rm = TRUE)
    rug_hi <- max(slope_draws[[spec$exposure]], na.rm = TRUE)
    if (is.null(rug_data)) rug_data <- prep_local_data(spec, get(spec$data))
    rug_df <- rug_data |>
      dplyr::filter(!is.na(.data[[spec$exposure]])) |>
      dplyr::filter(.data[[spec$exposure]] >= rug_lo,
                    .data[[spec$exposure]] <= rug_hi)
    p <- p + ggplot2::geom_rug(
      data = rug_df,
      ggplot2::aes(x = !!exp_sym),
      sides = "b", alpha = 0.3, length = ggplot2::unit(0.02, "npc"),
      inherit.aes = FALSE
    )
  }

  # Hide the x-axis text / title / ticks when this AMEF panel is not on the
  # bottom row of its column (shared-axes assembly).
  if (!show_x_axis) {
    p <- p + ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  }

  # Apply per-column x-limit (shared-axes assembly). coord_cartesian rather
  # than scale_x_continuous(limits=) so simul-band edges don't get clipped.
  if (!is.null(x_limits)) p <- p + ggplot2::coord_cartesian(xlim = x_limits)

  p
}


# Age-conditional AMEF panel with simultaneous credible bands per age slice.
# Renders posterior uncertainty for each age slice so apparent age-conditional
# divergences can be judged: wide bands in sparse data corners are tensor-smooth
# artifacts, tight bands that nevertheless differ between ages are signal. Bands
# only, with no posterior median drawn, matching the convention of every
# reported figure.
#
# Like make_amef_panel(), the y-axis is on the reported contrast scale.

make_age_cond_amef_panel_with_bands <- function(spec, amef_age_draws,
                                                y_label = NULL, x_axis_label = NULL,
                                                x_limits = NULL, show_x_axis = TRUE,
                                                show_legend = FALSE) {
  exp_sym <- rlang::sym(spec$exposure)
  cu      <- spec$contrast_units
  if (is.null(cu) || !is.finite(cu) || cu <= 0)
    stop("spec$contrast_units must be a positive number; got ",
         format(cu), " for exposure '", spec$exposure, "'.")

  # SIMULTANEOUS credible bands per age slice, not pointwise intervals: the
  # convention used by the main-text AMEF panels, calibrated for joint coverage
  # of the entire curve. Computed per age slice because each slice has its own
  # posterior of AMEF curves; pooling across ages would average over exactly the
  # conditioning this figure exists to expose.
  summary_df <- amef_age_draws |>
    dplyr::group_split(age_years) |>
    purrr::map_dfr(function(df) {
      bands <- simul_credible_bands(
        df, exposure = !!exp_sym,
        value_col = "draw", levels = 0.95,
        function_type = "AMEF"
      )
      bands$age_years <- df$age_years[1]
      bands
    }) |>
    dplyr::mutate(age_years = factor(age_years)) |>
    # Same conversion as make_amef_panel(), across the same columns: a
    # half-scaled bands object is a trap for whoever reads `med` next.
    dplyr::mutate(dplyr::across(dplyr::any_of(c("m", "sd", "med", "lo", "hi")),
                                function(v) v * cu))

  if (is.null(y_label))
    y_label <- paste0("Slope of predicted\n", spec$outcome_label,
                      "\n", spec$contrast_label)
  if (is.null(x_axis_label)) x_axis_label <- spec$exposure_label

  age_colors <- c("25" = "#0072B2",   # Okabe-Ito blue
                  "35" = "#009E73",   # bluish green
                  "50" = "#E69F00",   # orange
                  "65" = "#CC79A7")   # reddish purple

  p <- ggplot2::ggplot(summary_df,
                       ggplot2::aes(x = !!exp_sym,
                                    ymin = lo, ymax = hi,
                                    fill = age_years)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
    ggplot2::geom_ribbon(alpha = 0.22, color = NA) +
    ggplot2::scale_fill_manual(name  = "Age (years)", values = age_colors) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 1))) +
    x_scale_for(spec$exposure) +
    ggplot2::scale_y_continuous(labels = format_thousands) +
    ggplot2::labs(x = x_axis_label, y = y_label) +
    theme_pub() +
    ggplot2::theme(
      axis.title      = ggplot2::element_text(face = "plain", size = 9),
      legend.position = if (show_legend) "right" else "none"
    )

  if (!show_x_axis) {
    p <- p + ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  }

  if (!is.null(x_limits)) p <- p + ggplot2::coord_cartesian(xlim = x_limits)

  p
}


# Stack the AERF panel over the AMEF panel into a single subfigure. The
# AERF panel loses its x-axis (carried by the AMEF below); both panels are
# aligned vertically. Optional letter tag goes to the top-left of the AERF.

stack_subfig <- function(aerf, amef, tag = NULL) {
  aerf <- aerf + ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.text.x  = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  )
  if (!is.null(tag)) aerf <- aerf + ggplot2::labs(tag = tag)
  patchwork::wrap_plots(aerf, amef, ncol = 1, heights = c(1, 1)) &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "plain", size = 12))
}


# ---- aerf_amef_plot: the 3-panel working diagnostic figure ----
#
# Writes the AERF / AMEF / curvature composite that each of the nine analysis
# scripts saves to `outputs/figures/working/<stem>/aerf-amef.pdf`. It is the
# per-script diagnostic, and it shows the second derivative, which the
# manuscript figures deliberately omit.
#
# Axes stay in MODELLING units here -- per one unit of the exposure column --
# because this figure is for checking a fit, not for reporting an effect. The
# manuscript panels (make_amef_panel) are the ones that apply the reported
# contrast.

aerf_amef_plot <- function(spec, pred_draws, slope_draws,
                           curvature_draws = NULL, truth_df = NULL) {

  ex <- spec$exposure
  exp_sym <- rlang::sym(ex)

  if (is.null(curvature_draws)) {
    curvature_draws <- compute_curvature_draws(slope_draws, ex)
  }

  # X-axis ticks in natural units; see EXPOSURE_AXIS_SCALE.
  x_scale <- x_scale_for(ex)

  # Strip text combines panel role + y-axis quantity so each facet
  # carries both the name (AERF / AMEF / Curvature) and what's plotted
  # in one place. y-axis label is dropped (the strip carries it).
  panel_levels <- c(
    paste0("AERF: ",                       spec$outcome_label),
    paste0("AMEF: ",                       spec$amef_label),
    paste0("Curvature (slope of AMEF): ",  spec$curvature_label)
  )

  combined <- dplyr::bind_rows(
    pred_draws      |> dplyr::mutate(panel = panel_levels[1]),
    slope_draws     |> dplyr::mutate(panel = panel_levels[2]),
    curvature_draws |> dplyr::mutate(panel = panel_levels[3])
  ) |>
    dplyr::mutate(panel = factor(panel, levels = panel_levels))

  summary_df <- combined |>
    dplyr::group_by(panel, !!exp_sym) |>
    ggdist::median_hdci(draw, .width = c(0.5, 0.8, 0.95)) |>
    dplyr::ungroup()

  # Zero-reference line only on the AMEF and Curvature panels.
  zero_lines <- tibble::tibble(
    panel = factor(panel_levels[2:3], levels = panel_levels),
    yint  = 0
  )

  p <- ggplot2::ggplot(summary_df,
                       ggplot2::aes(x = !!exp_sym, y = draw,
                                    ymin = .lower, ymax = .upper)) +
    ggplot2::geom_hline(data = zero_lines,
                        ggplot2::aes(yintercept = yint),
                        linetype = "dashed", color = "grey50",
                        inherit.aes = FALSE) +
    ggdist::geom_lineribbon(alpha = 0.6) +
    ggplot2::scale_fill_brewer(palette = "Oranges") +
    x_scale +
    ggplot2::facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    ggplot2::labs(x = spec$exposure_label, y = NULL, fill = "HPDI") +
    theme_pub_leg() +
    ggplot2::theme(
      strip.text     = ggplot2::element_text(face = "bold", hjust = 0,
                                             size = ggplot2::rel(0.95)),
      panel.spacing  = ggplot2::unit(0.6, "lines")
    )

  if (!is.null(truth_df)) {
    truth_df <- truth_df |>
      dplyr::mutate(panel = factor(panel, levels = panel_levels))
    p <- p + ggplot2::geom_line(
      data = truth_df,
      ggplot2::aes(x = !!exp_sym, y = true),
      color = "black", linetype = "dashed", linewidth = 0.7,
      inherit.aes = FALSE
    )
  }
  p
}


# ---- Posterior predictive checks ----

# Six-panel PP check: min, mean, median, max, SD, skew.
pp_check_stats <- function(model, ndraws = 100) {
  p1 <- brms::pp_check(model, type = "stat", stat = "min",    ndraws = ndraws) + ggplot2::ggtitle("Min")    + theme_pub()
  p2 <- brms::pp_check(model, type = "stat", stat = "mean",   ndraws = ndraws) + ggplot2::ggtitle("Mean")   + theme_pub()
  p3 <- brms::pp_check(model, type = "stat", stat = "median", ndraws = ndraws) + ggplot2::ggtitle("Median") + theme_pub()
  p4 <- brms::pp_check(model, type = "stat", stat = "max",    ndraws = ndraws) + ggplot2::ggtitle("Max")    + theme_pub()
  p5 <- brms::pp_check(model, type = "stat", stat = "sd",     ndraws = ndraws) + ggplot2::ggtitle("SD")     + theme_pub()
  p6 <- brms::pp_check(model, type = "stat", stat = function(x) skew(x), ndraws = ndraws) + ggplot2::ggtitle("Skew") + theme_pub()
  (p1 + p2 + p3) / (p4 + p5 + p6)
}


# ---- simul_credible_bands: simultaneous credible bands (vector of levels) ----
#
# A band whose stated level is the probability that an ENTIRE posterior curve
# lies inside it, not the per-x probability. Two constructions, both calibrated
# by the same principle -- widen a one-parameter family until the fraction of
# whole curves contained reaches `level` -- but differing in the family:
#
#   "HPDI"       Nested pointwise highest-density intervals. At each grid point
#                the interval is the narrowest one spanning a common mass q of
#                that point's draws; q is searched upward from `level` until the
#                band contains whole curves at `level`. Asymmetric about the
#                median wherever the posterior is skewed, which is what makes it
#                an HPD region and what the manuscript's "simultaneous HPDI"
#                captions describe. This is the construction the reported
#                figures use.
#
#   "symmetric"  Standardized sup-norm: mean +/- c* x pointwise SD, with c* the
#                `level` quantile of the per-draw sup-norm statistic
#                max_x |y(x) - mean(x)| / sd(x). Symmetric about the MEAN by
#                construction, so it is not an HPD region however it is
#                calibrated. Kept so the two can be compared directly.
#
# The default comes from `getOption("bone.band_type")`, which `_startup/options.R`
# sets from the BAND_TYPE environment variable. Callers do not normally pass
# `interval_type` -- that is what lets one env var rebuild every simultaneous
# band the other way. `function_type` is metadata only; it is carried into the
# output for provenance and does not affect the computation.
#
# Returns a long-format tibble, one row per (exposure x level), for plotting
# nested bands on a single panel.

simul_credible_bands <- function(
    draws_df,
    exposure,
    value_col     = "draw",
    levels        = c(0.95),
    function_type = c("AERF", "AMEF"),
    interval_type = getOption("bone.band_type", "HPDI")
) {

  function_type <- match.arg(function_type)
  interval_type <- match.arg(interval_type, c("HPDI", "symmetric"))

  stopifnot(is.data.frame(draws_df))
  if (!("drawid" %in% names(draws_df))) {
    stop("draws_df must contain a column named 'drawid'.")
  }
  if (!(value_col %in% names(draws_df))) {
    stop(sprintf("draws_df must contain the value column '%s'.", value_col))
  }
  if (!is.numeric(levels) || any(!is.finite(levels)) || any(levels <= 0) || any(levels >= 1)) {
    stop("levels must be numeric values strictly between 0 and 1 (e.g., c(0.5, 0.8, 0.95)).")
  }

  exposure_sym <- rlang::ensym(exposure)
  exposure_nm  <- rlang::as_name(exposure_sym)
  if (!(exposure_nm %in% names(draws_df))) {
    stop(sprintf("draws_df must contain the exposure column '%s'.", exposure_nm))
  }

  levels <- sort(unique(levels))

  # Pivot to a draws x grid matrix once. Both constructions need whole curves,
  # not per-x marginals, so the long frame is never the right shape for either.
  keep <- !is.na(draws_df[[value_col]]) & !is.na(draws_df[[exposure_nm]])
  d      <- draws_df[keep, , drop = FALSE]
  x_grid <- sort(unique(d[[exposure_nm]]))
  ids    <- sort(unique(d[["drawid"]]))
  Y <- matrix(NA_real_, nrow = length(ids), ncol = length(x_grid))
  # Checked BEFORE the fill: a repeated (draw, exposure) pair overwrites its own
  # cell, so the matrix comes out complete and `anyNA` below sees nothing wrong.
  if (nrow(d) != length(ids) * length(x_grid)) {
    stop(sprintf(paste("draws contain %d rows for %d draws x %d grid points;",
                       "expected %d. Duplicated (draw, exposure) pairs would be",
                       "silently absorbed."),
                 nrow(d), length(ids), length(x_grid),
                 length(ids) * length(x_grid)))
  }
  Y[cbind(base::match(d[["drawid"]], ids), base::match(d[[exposure_nm]], x_grid))] <- d[[value_col]]
  if (anyNA(Y)) {
    stop(sprintf(paste("draws are not a complete draw x exposure grid:",
                       "%d of %d cells are missing."),
                 sum(is.na(Y)), length(Y)))
  }

  ng   <- length(x_grid)
  nlev <- length(levels)
  m    <- colMeans(Y)
  sdv  <- apply(Y, 2, stats::sd)
  med  <- apply(Y, 2, median)

  if (identical(interval_type, "symmetric")) {

    if (any(!is.finite(sdv)) || any(sdv <= 0)) {
      stop("Some exposure grid points have sd <= 0 across draws; ",
           "cannot form a standardized sup-norm band.")
    }
    z        <- sweep(sweep(Y, 2, m, "-"), 2, sdv, "/")
    t_sorted <- sort(apply(abs(z), 1, max))
    n_t      <- length(t_sorted)
    if (n_t == 0) stop("No finite sup-norm statistics found; cannot construct bands.")
    c_star <- t_sorted[pmax(1L, pmin(n_t, ceiling(levels * n_t)))]
    lo <- lapply(c_star, function(cs) m - cs * sdv)
    hi <- lapply(c_star, function(cs) m + cs * sdv)

  } else {

    Ys     <- apply(Y, 2, sort)     # one sort per grid point, reused by the search
    bounds <- lapply(levels, function(lv) .hpd_bounds(Ys, .simul_hpd_level(Y, Ys, lv)))
    lo <- lapply(bounds, function(b) b[1, ])
    hi <- lapply(bounds, function(b) b[2, ])

  }

  tibble::tibble(
    !!exposure_sym := rep(x_grid, times = nlev),
    m             = rep(m,   times = nlev),
    sd            = rep(sdv, times = nlev),
    med           = rep(med, times = nlev),
    lo            = unlist(lo, use.names = FALSE),
    hi            = unlist(hi, use.names = FALSE),
    level         = rep(levels, each = ng),
    function_type = function_type,
    interval_type = interval_type
  )
}

# Narrowest interval containing mass `q` of each column's draws.
#
# `Ys` must have each column sorted ascending. The interval spans
# floor(n * q) + 1 order statistics, and ties are broken toward the lower end,
# so the result is exact and deterministic. Returns a 2 x ncol(Ys) matrix.
#
# This is the same ESTIMATOR as `ggdist::hdci()` -- narrowest continuous
# interval carrying mass q -- computed differently, and the difference is worth
# knowing about. `hdci()` interpolates between order statistics (type-5
# quantiles) and locates the minimum with `optimize()`, a search over a
# piecewise-linear, non-convex width function; the version here takes the exact
# discrete minimum over order statistics. They agree to O(1/n): on 4,000
# lognormal draws the endpoints differ by ~0.3% of the interval width, which is
# far below the line width of any band drawn from them. The discrete form is
# used here because a band evaluates this tens of thousands of times per figure
# and must not inherit an optimizer's tolerance, or its local minima, into a
# published ribbon. Scalar HPDIs elsewhere use `ggdist::hdci()`.
.hpd_bounds <- function(Ys, q) {
  n <- nrow(Ys)
  m <- min(n, floor(n * q) + 1L)
  lo_i  <- seq_len(n - m + 1L)
  width <- Ys[lo_i + m - 1L, , drop = FALSE] - Ys[lo_i, , drop = FALSE]
  best  <- apply(width, 2L, which.min)
  cols  <- seq_len(ncol(Ys))
  rbind(Ys[cbind(best, cols)], Ys[cbind(best + m - 1L, cols)])
}

# Smallest pointwise HPD mass whose nested band contains whole curves at
# `level`. Containment is monotone non-decreasing in q, so bisection is exact
# to `tol`; the upper end of the final bracket is returned, which errs toward
# very slightly over-covering rather than under. The bracket is valid because
# q = level can only under-cover (simultaneous coverage never exceeds the
# pointwise level) and q = 1 spans every draw at every grid point.
.simul_hpd_level <- function(Y, Ys, level, tol = 1e-4) {
  n  <- nrow(Y)
  ng <- ncol(Y)
  # Named frac_inside, not `contains`: tidyselect exports a `contains()` that
  # the project attaches, and a local binding shadowing it would read as a
  # mistake even where it is not one.
  frac_inside <- function(q) {
    b <- .hpd_bounds(Ys, q)
    inside <- Y >= rep(b[1, ], each = n) & Y <= rep(b[2, ], each = n)
    mean(rowSums(inside) == ng)
  }
  lo <- level
  hi <- 1
  if (frac_inside(lo) >= level) return(lo)
  while (hi - lo > tol) {
    mid <- (lo + hi) / 2
    if (frac_inside(mid) >= level) hi <- mid else lo <- mid
  }
  hi
}
