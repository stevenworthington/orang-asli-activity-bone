###############################################################################
# Causal AERF/AMEF for tibial SOS (m/s, rescaled to /1000) given daily step
# count (rescaled to /1000). t2(age, steps) tensor smooth + DAG-implied
# adjustment {sex, pregnant_or_breastfeeding_n_y_0_1,
# functional_status_n_y_0_1}. Reads formula / priors / MCMC config from the
# registry in `code/_startup/specifications.R`.
###############################################################################


library(here)
source(here("code", "_startup", "init.R"))

# Smoke-test mode: `SAMPLING_MODE=smoke Rscript code/sos-steps.R` overrides
# the registry's full sampling settings with a fast 4-chain run.
SAMPLING_MODE <- Sys.getenv("SAMPLING_MODE", "full")
if (SAMPLING_MODE == "smoke") {
  WARMUP <- 500
  ITER   <- 1500
  CHAINS <- 4
}
cat("SAMPLING_MODE =", SAMPLING_MODE,
    " | WARMUP =", WARMUP, "ITER =", ITER, "CHAINS =", CHAINS, "\n")

set.seed(SEED)

SCRIPT_STEM <- "sos-steps"
spec        <- model_templates[[SCRIPT_STEM]]
dat_local   <- get(spec$data) |>
  tidyr::drop_na(tidyselect::all_of(c(spec$outcome, spec$exposure,
                                      "age_years", "sex",
                                      "pregnant_or_breastfeeding_n_y_0_1",
                                      "functional_status_n_y_0_1")))

cat("Fitting", SCRIPT_STEM, "on n =", nrow(dat_local), "rows\n")

out_dir <- here("outputs", "models", SCRIPT_STEM)
fig_dir <- here("outputs", "figures", "working", SCRIPT_STEM)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)


# ---- Fit ----

model_fit <- brm(
  spec$bf,
  data      = dat_local,
  prior     = spec$priors,
  warmup    = WARMUP, iter = ITER, thin = THIN,
  chains    = CHAINS, cores = CHAINS,
  seed      = SEED,
  backend   = "cmdstanr",
  control   = BRMS_CONTROL,
  save_pars = save_pars(all = TRUE),
  refresh   = 200
)

save(model_fit, file = file.path(out_dir, "model.Rdata"), compress = "gzip")
cat("\nSaved model.Rdata\n")


# ---- Convergence summary ----

cat("\n--- summary(model_fit) ---\n")
print(summary(model_fit))

cat("\n--- LOO ---\n")
tryCatch(print(loo(model_fit)),
         error = function(e) cat("LOO failed:", conditionMessage(e),
                                 "-- skipping; not load-bearing for AERF/AMEF\n"))


# ---- Posterior-predictive checks ----

ggsave(pp_check_stats(model_fit, ndraws = 100),
       file = file.path(fig_dir, "pp-check-stats.pdf"), height = 5, width = 7)


# ---- AERF on 51-point exposure grid (1st-99th percentile of observed) ----

grid <- seq(
  quantile(dat_local[[spec$exposure]], spec$grid_quantiles[1], na.rm = TRUE),
  quantile(dat_local[[spec$exposure]], spec$grid_quantiles[2], na.rm = TRUE),
  length.out = 51
)

# Build a 1-row "typical observation" via datagrid() (numeric covariates at
# their means, factors at their modes), then replicate over the exposure grid.
typical  <- marginaleffects::datagrid(model = model_fit)
new_grid <- typical[rep(1, length(grid)), , drop = FALSE]
new_grid[[spec$exposure]] <- grid

avg_predictions(
  model_fit,
  newdata = new_grid,
  by      = spec$exposure
) |>
  marginaleffects::posterior_draws() ->
pred_draws

# Backtransform to natural units so downstream consumers (figures,
# simul-bands, flatness / linearity probs) read natural-scale draws.
pred_draws$draw <- pred_draws$draw * spec$outcome_scale_factor

save(pred_draws, file = file.path(out_dir, "pred-draws.Rdata"), compress = "gzip")
cat("Saved pred-draws.Rdata (", nrow(pred_draws), "rows )\n")


# ---- AMEF (first derivative) on the same grid ----

avg_slopes(
  model_fit,
  variables = spec$exposure,
  newdata   = new_grid,
  by        = spec$exposure
) |>
  marginaleffects::posterior_draws() ->
slope_draws

slope_draws$draw <- slope_draws$draw * spec$outcome_scale_factor

save(slope_draws, file = file.path(out_dir, "slope-draws.Rdata"), compress = "gzip")
cat("Saved slope-draws.Rdata (", nrow(slope_draws), "rows )\n")


# ---- Curvature (d^2 AERF / dexposure^2) from finite difference on AMEF ----

curvature_draws <- compute_curvature_draws(slope_draws, spec$exposure)
save(curvature_draws, file = file.path(out_dir, "curvature-draws.Rdata"), compress = "gzip")
cat("Saved curvature-draws.Rdata (", nrow(curvature_draws), "rows )\n")


# ---- Working figures: AERF / AMEF / curvature stacked with HPDI ribbons ----

ggsave(aerf_amef_plot(spec, pred_draws, slope_draws, curvature_draws),
       file = file.path(fig_dir, "aerf-amef.pdf"),
       height = 8.4, width = 4)
cat("Saved aerf-amef.pdf\n")

cat("\nDone:", SCRIPT_STEM, "\n")
