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
###############################################################################


# ---- Drop rows with NA in any covariate the spec references ----

prep_local_data <- function(spec, dat_raw) {
  predictor_terms <- all.vars(spec$bf$formula)
  drop_cols <- unique(c(spec$outcome, spec$exposure, predictor_terms))
  drop_cols <- intersect(drop_cols, names(dat_raw))
  dat_raw |>
    tidyr::drop_na(tidyselect::all_of(drop_cols))
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
    warmup_local <- 500
    iter_local   <- 1500
    chains_local <- 4
  } else {
    warmup_local <- WARMUP
    iter_local   <- ITER
    chains_local <- CHAINS
  }

  dat_local <- prep_local_data(spec, dat_raw)

  brm(
    spec$bf,
    data      = dat_local,
    prior     = spec$priors,
    warmup    = warmup_local,
    iter      = iter_local,
    thin      = THIN,
    chains    = chains_local,
    cores     = chains_local,
    seed      = SEED,
    backend   = "cmdstanr",
    control   = BRMS_CONTROL,
    save_pars = save_pars(all = TRUE),
    refresh   = 0,
    silent    = 2
  )
}


# ---- AERF posterior draws on the spec's exposure grid ----
#
# Uses a counterfactual datagrid -- the full analytic dataset gets replicated
# for each exposure-grid value, so `avg_predictions(by = exposure)` truly
# averages over the cohort's observed covariate distribution. This matches the
# manuscript's Watson realized-causal-inference target of inference ("the
# adjustment-set covariates distributed as observed in the cohort") and the
# cross-language convention in `~/.config/agents/CODING.md` §
# Postestimation. The earlier `datagrid(model = fit)` default
# (`grid_type = "mean_or_mode"`) gave AERFs at typical covariate values
# rather than the population-average AERF -- corrected 2026-05-26.

aerf_draws <- function(fit, spec, dat_raw) {
  dat_local <- prep_local_data(spec, dat_raw)
  grid <- seq(
    quantile(dat_local[[spec$exposure]], spec$grid_quantiles[1], na.rm = TRUE),
    quantile(dat_local[[spec$exposure]], spec$grid_quantiles[2], na.rm = TRUE),
    length.out = 51
  )
  dg_args <- list(model = fit, grid_type = "counterfactual")
  dg_args[[spec$exposure]] <- grid
  new_grid <- do.call(marginaleffects::datagrid, dg_args)

  pred <- marginaleffects::avg_predictions(
    fit, newdata = new_grid, by = spec$exposure
  ) |>
    marginaleffects::posterior_draws()
  pred$draw <- pred$draw * spec$outcome_scale_factor
  pred
}


# ---- AMEF (first-derivative) posterior draws on the same grid ----

amef_draws <- function(fit, spec, dat_raw) {
  dat_local <- prep_local_data(spec, dat_raw)
  grid <- seq(
    quantile(dat_local[[spec$exposure]], spec$grid_quantiles[1], na.rm = TRUE),
    quantile(dat_local[[spec$exposure]], spec$grid_quantiles[2], na.rm = TRUE),
    length.out = 51
  )
  dg_args <- list(model = fit, grid_type = "counterfactual")
  dg_args[[spec$exposure]] <- grid
  new_grid <- do.call(marginaleffects::datagrid, dg_args)

  slopes <- marginaleffects::avg_slopes(
    fit, variables = spec$exposure,
    newdata = new_grid, by = spec$exposure
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
# fixed at that value and the exposure varying over the spec's 1st-99th-
# percentile grid; query `avg_slopes` with `by = exposure` to recover the
# posterior of the per-exposure AMEF *conditional on age = age_value*.
# Used by the age-conditional AMEF age-conditional reporting per Ian's 2026-05-25 ask
# ("are PA effects largest pre-peak-bone-mass?"). Returns one long-format
# tibble with all age slices stacked, plus an `age_years` column.
#
# Note: this re-queries avg_slopes on the pre-existing fit; no refitting.
# 24 calls (6 specs x 4 ages) typically take ~30-60 s total.

amef_at_age <- function(fit, spec, age_value, dat_raw) {
  dat_local <- prep_local_data(spec, dat_raw)
  grid <- seq(
    quantile(dat_local[[spec$exposure]], spec$grid_quantiles[1], na.rm = TRUE),
    quantile(dat_local[[spec$exposure]], spec$grid_quantiles[2], na.rm = TRUE),
    length.out = 51
  )
  # Counterfactual datagrid: full analytic dataset replicated for each
  # (age_value x exposure) combination. With age_years fixed at a single
  # value and exposure varying across 51 points, this gives n x 51 rows.
  # `avg_slopes(by = exposure)` then averages over the cohort's other
  # covariates within each exposure level -- the age-conditional analog of
  # the population-average AMEF in amef_draws(). Matches the convention in
  # ~/.config/agents/CODING.md § Postestimation.
  dg_args <- list(model = fit, grid_type = "counterfactual",
                  age_years = age_value)
  dg_args[[spec$exposure]] <- grid
  new_grid <- do.call(marginaleffects::datagrid, dg_args)

  slopes <- marginaleffects::avg_slopes(
    fit, variables = spec$exposure,
    newdata = new_grid, by = spec$exposure
  ) |>
    marginaleffects::posterior_draws()
  slopes$draw <- slopes$draw * spec$outcome_scale_factor
  slopes$age_years <- age_value
  slopes
}


# ---- Headline summary numbers from a single fit's draws ----

summarize_one_fit <- function(spec_key, pred_draws, slope_draws, curvature_draws,
                              n_analytic = NA_integer_) {
  spec <- model_templates[[spec_key]]
  ex   <- spec$exposure
  ex_lo <- min(slope_draws[[ex]])
  ex_hi <- max(slope_draws[[ex]])

  d_lo <- slope_draws |> dplyr::filter(.data[[ex]] == ex_lo) |>
                         dplyr::arrange(drawid) |> dplyr::pull(draw)
  d_hi <- slope_draws |> dplyr::filter(.data[[ex]] == ex_hi) |>
                         dplyr::arrange(drawid) |> dplyr::pull(draw)
  cu_hi <- curvature_draws |>
    dplyr::filter(.data[[ex]] == max(curvature_draws[[ex]])) |>
    dplyr::arrange(drawid) |> dplyr::pull(draw)

  hpdi_lo <- HDInterval::hdi(d_lo, credMass = 0.95)
  hpdi_hi <- HDInterval::hdi(d_hi, credMass = 0.95)

  data.frame(
    spec_key        = spec_key,
    n_analytic      = as.integer(n_analytic),
    p_declines      = round(mean(d_hi < d_lo), 3),
    p_concave_at_hi = round(mean(cu_hi < 0), 3),
    amef_lo_median  = round(unname(median(d_lo)), 4),
    amef_lo_hpdi_lo = round(unname(hpdi_lo[1]), 4),
    amef_lo_hpdi_hi = round(unname(hpdi_lo[2]), 4),
    amef_hi_median  = round(unname(median(d_hi)), 4),
    amef_hi_hpdi_lo = round(unname(hpdi_hi[1]), 4),
    amef_hi_hpdi_hi = round(unname(hpdi_hi[2]), 4)
  )
}
