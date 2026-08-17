###############################################################################
# Pipeline-callable functions. Used by `_targets.R` so each modeling step is
# a function call with explicit inputs/outputs (which `targets` can hash for
# content-level skip-if-up-to-date), rather than a script with implicit
# side effects.
#
# These functions assume the environment created by sourcing
# _startup/init.R (which sources functions.R, packages.R, options.R, data.R,
# specifications.R, this file) -- so `model_templates`, `SEED`, `WARMUP`,
# `ITER`, `THIN`, `CHAINS`, `BRMS_CONTROL`, and the helpers from functions.R
# are available.
#
# It also holds the READ side: `spec_draws()` and `spec_panel_data()`, the way
# every reported figure and table gets a fitted spec's artifacts back out of the
# cache. Producing and consuming the pipeline's outputs belong together -- that
# is what keeps `functions.R` free of any I/O.
#
# One deliberate exception: `_experiments/effect-size-probability-fig.R` reads
# `_targets/objects/` directly with qs2, because it is built to run without the
# startup set attached. It reads the same objects these helpers do, so the two
# cannot disagree.
#
# UNITS. Every draw frame produced here is in NATURAL outcome units (via
# `outcome_scale_factor`) and per ONE unit of the modelled exposure column. The
# reported contrast is applied downstream, at the presentation boundary -- see
# the `contrast_units` note in specifications.R.
###############################################################################


# ---- Drop rows with NA in any covariate the spec references ----

prep_local_data <- function(spec, dat_raw) {
  predictor_terms <- all.vars(spec$bf$formula)
  drop_cols <- unique(c(spec$outcome, spec$exposure, predictor_terms))
  drop_cols <- base::intersect(drop_cols, names(dat_raw))
  dat_raw |>
    tidyr::drop_na(tidyselect::all_of(drop_cols)) |>
    # droplevels() matters: dropping rows leaves the factor's full level set
    # attached, so the analytic frame can advertise 31 communities while the
    # model sees 24. marginaleffects then builds a prediction grid carrying
    # levels the fit never saw and errors ("New factor levels are not
    # allowed", surfaced as a generic "unable to compute predicted values").
    # Dropping unused levels is a no-op on the fit -- an unused level
    # contributes no design column -- and makes the fit, its stored data and
    # any prediction grid agree.
    droplevels()
}


# ---- The exposure grid a spec's AERF / AMEF is evaluated on ----
#
# One definition, used by aerf_draws(), amef_draws() and amef_at_age(), so the
# three cannot disagree about where the curve is evaluated.
#
# `spec$grid_cluster`, when set, names a grouping column the exposure is
# CONSTANT within -- the industrialization index is a community attribute, so
# every resident of a community shares one value. Quantiles must then be taken
# over the distinct community values, not over individuals: weighting by
# community size lets the largest communities determine the support of the
# reported curve. Getting this wrong is not cosmetic: on the community-level
# analyses it shifts both ends of the reported interval by enough to change
# whether it excludes zero.
#
# `spec$grid_n` is the number of points: 51 for the individual-level specs, 41
# for the community-level ones.

spec_grid <- function(spec, dat_local) {
  x <- dat_local[[spec$exposure]]
  if (!is.null(spec$grid_cluster)) {
    by_cluster <- split(x, dat_local[[spec$grid_cluster]], drop = TRUE)
    # Taking v[1] per cluster is only the community's index if the exposure
    # really is constant within the cluster. Asserted rather than assumed: if it
    # stopped holding, the grid would silently be built from arbitrary rows.
    varying <- names(by_cluster)[vapply(by_cluster,
                 function(v) length(unique(v[!is.na(v)])) > 1L, logical(1))]
    if (length(varying))
      stop("'", spec$exposure, "' is not constant within ", spec$grid_cluster,
           " for: ", paste(varying, collapse = ", "),
           ". spec_grid() assumes a cluster-level exposure.")
    x <- as.numeric(vapply(by_cluster, function(v) v[1], numeric(1)))
  }
  n <- if (is.null(spec$grid_n)) 51L else spec$grid_n
  seq(quantile(x, spec$grid_quantiles[1], na.rm = TRUE),
      quantile(x, spec$grid_quantiles[2], na.rm = TRUE),
      length.out = n)
}


# ---- Arguments a spec wants forwarded to postestimation ----
#
# Currently just `re_formula`. The three industrialization specs set NA, making
# predictions marginal of the community random intercept -- the across-community
# estimand the manuscript claims. Left unset (and so not forwarded) for the PA
# specs, which have no random effects at all. See specifications.R.

postestimation_args <- function(spec) {
  if (is.null(spec$re_formula)) list() else list(re_formula = spec$re_formula)
}


# ---- Fit a single spec by key ----
#
# spec_key: a key in `model_templates` (e.g. "sos-steps", "ctx-enmo", "sos-urb").
# dat_raw:  the project dataset (already loaded into `dat` by data.R).
# mode:     "full" = production sampling, "smoke" = fast 4-chain run for
#           plumbing tests.
#
# Returns:  a fitted brms object.

fit_one_spec <- function(spec_key, dat_raw, mode = "full") {
  spec <- model_templates[[spec_key]]
  if (is.null(spec)) stop("Unknown spec key: ", spec_key)

  if (identical(mode, "smoke")) {
    warmup_local <- SMOKE_WARMUP
    iter_local   <- SMOKE_ITER
    chains_local <- SMOKE_CHAINS
  } else {
    warmup_local <- WARMUP
    iter_local   <- ITER
    chains_local <- CHAINS
  }

  dat_local <- prep_local_data(spec, dat_raw)

  # Priors are autoscaled to the ANALYTIC data, so they are built here rather
  # than carried as a static field on the spec. See build_priors() in
  # _startup/functions.R for the scheme.
  priors_local <- build_priors(spec, dat_local)

  brms::brm(
    spec$bf,
    data      = dat_local,
    prior     = priors_local,
    warmup    = warmup_local,
    iter      = iter_local,
    thin      = THIN,
    chains    = chains_local,
    cores     = chains_local,
    seed      = SEED,
    backend   = "cmdstanr",
    control   = BRMS_CONTROL,
    save_pars = brms::save_pars(all = TRUE),
    refresh   = 0,
    silent    = 2
  )
}


# ---- AERF posterior draws on the spec's exposure grid ----
#
# Uses a counterfactual datagrid -- the full analytic dataset gets replicated
# for each exposure-grid value, so `avg_predictions(by = exposure)` truly
# averages over the cohort's observed covariate distribution. This matches the
# manuscript's Watson realized-causal-inference target of inference: the
# adjustment-set covariates distributed as observed in the cohort. The
# alternative, marginaleffects' `datagrid(model = fit)` default
# (`grid_type = "mean_or_mode"`), would give the AERF at typical covariate
# values rather than the population-average AERF.

aerf_draws <- function(fit, spec, dat_raw) {
  dat_local <- prep_local_data(spec, dat_raw)
  grid <- spec_grid(spec, dat_local)
  dg_args <- list(model = fit, grid_type = "counterfactual")
  dg_args[[spec$exposure]] <- grid
  new_grid <- do.call(marginaleffects::datagrid, dg_args)

  pred <- do.call(
    marginaleffects::avg_predictions,
    c(list(fit, newdata = new_grid, by = spec$exposure),
      postestimation_args(spec))
  ) |>
    marginaleffects::posterior_draws()
  pred$draw <- pred$draw * spec$outcome_scale_factor
  pred
}


# ---- AMEF (first-derivative) posterior draws on the same grid ----

amef_draws <- function(fit, spec, dat_raw) {
  dat_local <- prep_local_data(spec, dat_raw)
  grid <- spec_grid(spec, dat_local)
  dg_args <- list(model = fit, grid_type = "counterfactual")
  dg_args[[spec$exposure]] <- grid
  new_grid <- do.call(marginaleffects::datagrid, dg_args)

  slopes <- do.call(
    marginaleffects::avg_slopes,
    c(list(fit, variables = spec$exposure, newdata = new_grid,
           by = spec$exposure),
      postestimation_args(spec))
  ) |>
    marginaleffects::posterior_draws()
  slopes$draw <- slopes$draw * spec$outcome_scale_factor
  slopes
}


# ---- Curvature (second-derivative) draws via finite difference on AMEF ----

curvature_draws_from_amef <- function(slope_draws, spec) {
  compute_curvature_draws(slope_draws, spec$exposure)
}


# ---- Age-conditional AMEF posterior draws ----
#
# For each age value supplied, build a `datagrid()`-style newdata with age
# fixed at that value and the exposure varying over the spec's grid; query
# `avg_slopes` with `by = exposure` to recover the posterior of the
# per-exposure AMEF *conditional on age = age_value*. This is what the
# supplementary age-conditional figure is built from: whether the activity
# effect is larger before peak bone mass. Returns one long-format tibble with
# all age slices stacked, plus an `age_years` column.
#
# Note: this re-queries avg_slopes on the pre-existing fit; no refitting.
# 24 calls (6 specs x 4 ages) typically take ~30-60 s total.

amef_at_age <- function(fit, spec, age_value, dat_raw) {
  dat_local <- prep_local_data(spec, dat_raw)
  grid <- spec_grid(spec, dat_local)
  # Counterfactual datagrid: full analytic dataset replicated for each
  # (age_value x exposure) combination. With age_years fixed at a single value
  # and the exposure varying across the spec's grid, this gives n x length(grid)
  # rows. `avg_slopes(by = exposure)` then averages over the cohort's other
  # covariates within each exposure level -- the age-conditional analog of the
  # population-average AMEF in amef_draws().
  dg_args <- list(model = fit, grid_type = "counterfactual",
                  age_years = age_value)
  dg_args[[spec$exposure]] <- grid
  new_grid <- do.call(marginaleffects::datagrid, dg_args)

  slopes <- do.call(
    marginaleffects::avg_slopes,
    c(list(fit, variables = spec$exposure, newdata = new_grid,
           by = spec$exposure),
      postestimation_args(spec))
  ) |>
    marginaleffects::posterior_draws()
  slopes$draw <- slopes$draw * spec$outcome_scale_factor
  slopes$age_years <- age_value
  slopes
}


# ---- Reading the pipeline's output back ----
#
# THE single way to get a fitted spec's artifacts. Every reported figure and
# table reads through here, so no consumer can be looking at a different fit
# from its neighbour, and nothing depends on a hand-maintained copy of the
# cache that could go stale without erroring.
#
# `tar_map()` renames hyphens to dots when it expands target names, which is the
# only reason the key needs translating at all.

spec_draws <- function(key, what = c("pred", "slope", "curvature", "fit")) {
  what <- match.arg(what)
  if (!key %in% names(model_templates))
    stop("unknown spec key '", key, "'. Registry holds: ",
         paste(names(model_templates), collapse = ", "))
  prefix <- switch(what,
                   pred      = "pred_draws_",
                   slope     = "slope_draws_",
                   curvature = "curvature_draws_",
                   fit       = "fit_")
  target <- paste0(prefix, gsub("-", ".", key, fixed = TRUE))
  if (!isTRUE(targets::tar_exist_objects(target)))
    stop("target '", target, "' is not in the targets cache.\n",
         "  Build it with `just fit-all`, or `just fit-one ", target, "`.")
  targets::tar_read_raw(target)
}


# ---- Everything a manuscript panel needs for one spec ----
#
# The AERF and AMEF draws plus the derived bands and linear projection, in the
# shape `make_aerf_panel()` / `make_amef_panel()` consume. Shared by
# `_final/figures.R` (Figure 4 and the supplementary grids) and
# `_final/figure-industrialization.R` (Figure 3), so the two figures cannot
# drift apart in how their panels are built.
#
# THE AMEF COMES FROM `avg_slopes`, via the `slope_draws` target, for every
# spec -- never from a finite difference on the AERF. On the 41-point community
# grid a central difference spans 1.85 index units, 18.5% of the reported
# 10-unit contrast, wide enough to smooth away curvature the model actually
# estimated; and it would make two figures in the same paper answer the same
# question two different ways.
#
# `band_levels` defaults here rather than in each caller so the nested ribbon
# set is defined once. Named that rather than `levels` so it does not shadow
# base::levels() inside the function.

spec_panel_data <- function(key,
                            band_levels = c(0.05, 0.25, 0.50, 0.75, 0.95),
                            proj_level = 0.95) {
  # spec_draws() validates the key, and is called first so an unknown one gets
  # that message rather than an opaque failure inside rlang::sym(NULL).
  pred_draws  <- spec_draws(key, "pred")
  slope_draws <- spec_draws(key, "slope")
  spec        <- model_templates[[key]]
  exp_sym     <- rlang::sym(spec$exposure)
  list(
    spec             = spec,
    pred_draws       = pred_draws,
    slope_draws      = slope_draws,
    simul_bands_aerf = simul_credible_bands(pred_draws,  exposure = !!exp_sym,
                                            levels = band_levels, function_type = "AERF"),
    simul_bands_amef = simul_credible_bands(slope_draws, exposure = !!exp_sym,
                                            levels = band_levels, function_type = "AMEF"),
    lin_proj         = linear_projection(pred_draws, exposure = !!exp_sym,
                                         level = proj_level)
  )
}


# ---- Headline summary numbers from a single fit's draws ----
#
# The AMEF quantities are put on the reported contrast here, and the `contrast`
# column names the denominator, so `outputs/tables/spec-summary.csv` cannot be
# read against the wrong one. `p_concave_at_hi` is a sign probability and is
# invariant to the scaling; the curvature draws stay in modelling units.
#
# Interval primitive is `ggdist::hdci` -- sample-based, narrowest continuous
# interval -- the project-wide choice, and the one `linear_projection()` uses
# for the reported slopes. It is not interchangeable with `HDInterval::hdi`:
# `hdi` spans floor(n * q) + 1 order statistics and so carries a little more
# than q, while `hdci` interpolates between order statistics and carries
# exactly q. On 4,000 lognormal draws the endpoints differ by ~0.3% of the
# interval width. Nothing in the manuscript comes from this table -- the only
# column read downstream is `n_analytic` -- but the two primitives should not be
# swapped casually elsewhere.

summarize_one_fit <- function(spec_key, slope_draws, curvature_draws,
                              n_analytic = NA_integer_) {
  spec <- model_templates[[spec_key]]
  ex   <- spec$exposure
  cu   <- spec$contrast_units
  ex_lo <- min(slope_draws[[ex]])
  ex_hi <- max(slope_draws[[ex]])

  d_lo <- slope_draws |> dplyr::filter(.data[[ex]] == ex_lo) |>
                         dplyr::arrange(drawid) |> dplyr::pull(draw)
  d_hi <- slope_draws |> dplyr::filter(.data[[ex]] == ex_hi) |>
                         dplyr::arrange(drawid) |> dplyr::pull(draw)
  cu_hi <- curvature_draws |>
    dplyr::filter(.data[[ex]] == max(curvature_draws[[ex]])) |>
    dplyr::arrange(drawid) |> dplyr::pull(draw)

  hpdi_lo <- ggdist::hdci(d_lo * cu, .width = 0.95)
  hpdi_hi <- ggdist::hdci(d_hi * cu, .width = 0.95)

  data.frame(
    spec_key        = spec_key,
    n_analytic      = as.integer(n_analytic),
    contrast        = spec$contrast_label,
    p_declines      = round(mean(d_hi < d_lo), 3),
    p_concave_at_hi = round(mean(cu_hi < 0), 3),
    amef_lo_median  = round(unname(median(d_lo * cu)), 4),
    amef_lo_hpdi_lo = round(unname(hpdi_lo[1]), 4),
    amef_lo_hpdi_hi = round(unname(hpdi_lo[2]), 4),
    amef_hi_median  = round(unname(median(d_hi * cu)), 4),
    amef_hi_hpdi_lo = round(unname(hpdi_hi[1]), 4),
    amef_hi_hpdi_hi = round(unname(hpdi_hi[2]), 4)
  )
}
