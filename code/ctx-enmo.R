###############################################################################
# Causal AERF/AMEF for serum CTX-1 (ng/mL, unscaled, Gamma with log link)
# given mean daily ENMO (scaled /10, i.e. per 10 mg).
# t2(age, ENMO) tensor smooth + DAG-implied adjustment set
# {sex, pregnancy/lactation, smoking, alcohol, functional status,
#  community identity}.
#
# Community indicators are treated as FIXED effects and their priors are left
# FLAT, so the effect is identified strictly from within-community contrasts.
# Reads formula / family / prior settings / MCMC config from the registry in
# `code/_startup/specifications.R`; priors are built at fit time by
# build_priors(), with kappa 0.25 on the coefficient block.
###############################################################################


source(here::here("code", "_startup", "init.R"))

# Smoke-test mode: `SAMPLING_MODE=smoke Rscript code/ctx-enmo.R` overrides
# the registry's full sampling settings with a fast 4-chain run.
SAMPLING_MODE <- Sys.getenv("SAMPLING_MODE", "full")
if (SAMPLING_MODE == "smoke") {
  WARMUP <- SMOKE_WARMUP
  ITER   <- SMOKE_ITER
  CHAINS <- SMOKE_CHAINS
}
cat("SAMPLING_MODE =", SAMPLING_MODE,
    " | WARMUP =", WARMUP, "ITER =", ITER, "CHAINS =", CHAINS, "\n")

set.seed(SEED)

SCRIPT_STEM <- "ctx-enmo"
spec        <- model_templates[[SCRIPT_STEM]]
dat_local   <- prep_local_data(spec, get(spec$data))

cat("Fitting", SCRIPT_STEM, "on n =", nrow(dat_local), "rows\n")

out_dir <- here::here("outputs", "models", SCRIPT_STEM)
fig_dir <- here::here("outputs", "figures", "working", SCRIPT_STEM)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)


# ---- Fit ----

model_fit <- brms::brm(
  spec$bf,
  data      = dat_local,
  prior     = build_priors(spec, dat_local),
  warmup    = WARMUP, iter = ITER, thin = THIN,
  chains    = CHAINS, cores = CHAINS,
  seed      = SEED,
  backend   = "cmdstanr",
  control   = BRMS_CONTROL,
  save_pars = brms::save_pars(all = TRUE),
  refresh   = 200
)

save(model_fit, file = file.path(out_dir, "model.Rdata"), compress = "gzip")
cat("\nSaved model.Rdata\n")

# Remove the global dat_local before postestimation. marginaleffects/insight
# recovers a model's data by name, and can pick up a same-named global variable
# in preference to the model's own internal copy. prep_local_data() calls
# droplevels(), so the two now agree on village_id's level set and the recovery
# is harmless where it happens -- but relying on that is relying on a helper
# elsewhere never changing. This cannot arise in the targets pipeline, where the
# equivalent variable is local to fit_one_spec() and out of scope by this point.
rm(dat_local)


# ---- Convergence summary ----

cat("\n--- summary(model_fit) ---\n")
print(summary(model_fit))

cat("\n--- LOO ---\n")
tryCatch(print(brms::loo(model_fit)),
         error = function(e) cat("LOO failed:", conditionMessage(e),
                                 "-- skipping; not load-bearing for AERF/AMEF\n"))


# ---- Posterior-predictive checks ----

ggplot2::ggsave(pp_check_stats(model_fit, ndraws = 100),
       file = file.path(fig_dir, "pp-check-stats.pdf"), height = 5, width = 7)


# ---- AERF / AMEF / curvature via the shared pipeline-helpers.R functions ----
#
# Calls the same aerf_draws() / amef_draws() the targets pipeline calls, rather
# than repeating the postestimation logic inline. Two copies of that logic can
# disagree about the datagrid they build -- population-average over the
# cohort's observed covariates, versus a single typical observation -- and then
# produce different numbers from the same fit with neither path erroring.

pred_draws <- aerf_draws(model_fit, spec, get(spec$data))
save(pred_draws, file = file.path(out_dir, "pred-draws.Rdata"), compress = "gzip")
cat("Saved pred-draws.Rdata (", nrow(pred_draws), "rows )\n")

slope_draws <- amef_draws(model_fit, spec, get(spec$data))
save(slope_draws, file = file.path(out_dir, "slope-draws.Rdata"), compress = "gzip")
cat("Saved slope-draws.Rdata (", nrow(slope_draws), "rows )\n")

curvature_draws <- curvature_draws_from_amef(slope_draws, spec)
save(curvature_draws, file = file.path(out_dir, "curvature-draws.Rdata"), compress = "gzip")
cat("Saved curvature-draws.Rdata (", nrow(curvature_draws), "rows )\n")


# ---- Working figures: AERF / AMEF / curvature stacked with HPDI ribbons ----

ggplot2::ggsave(aerf_amef_plot(spec, pred_draws, slope_draws, curvature_draws),
       file = file.path(fig_dir, "aerf-amef.pdf"),
       height = 8.4, width = 4)
cat("Saved aerf-amef.pdf\n")

cat("\nDone:", SCRIPT_STEM, "\n")
