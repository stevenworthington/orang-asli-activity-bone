###############################################################################
# Causal AERF/AMEF for serum CTX-1 (ng/mL, lognormal) given mean daily
# ENMO (mg). t2(age, ENMO) tensor smooth + DAG-implied adjustment
# {sex, pregnancy/lactation, smoking, alcohol, functional status, village identity}. Reads
# formula / priors / MCMC config from the registry in
# `code/_startup/specifications.R`.
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

SCRIPT_STEM <- "ctx-enmo"
spec        <- model_templates[[SCRIPT_STEM]]
dat_local   <- prep_local_data(spec, get(spec$data))

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

# Remove the global dat_local before postestimation. village_id's factor
# levels reflect the full 31-village dataset even after drop_na() filters
# rows down to this spec's complete cases (drop_na doesn't call droplevels()),
# so dat_local's village_id carries more levels than the fitted model actually
# used. marginaleffects/insight's data recovery can pick up a same-named
# global variable instead of the model's own (clean) internal data, and then
# the stale, wider level set collides with the model's narrower one --
# "New factor levels are not allowed." This can't happen in the targets
# pipeline, where the equivalent variable is local to fit_one_spec() and goes
# out of scope before postestimation runs.
rm(dat_local)


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


# ---- AERF / AMEF / curvature via the shared pipeline-helpers.R functions ----
#
# Calls the same aerf_draws() / amef_draws() functions the targets pipeline
# uses, so the standalone-script and pipeline paths are guaranteed identical
# rather than merely pattern-identical.

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

ggsave(aerf_amef_plot(spec, pred_draws, slope_draws, curvature_draws),
       file = file.path(fig_dir, "aerf-amef.pdf"),
       height = 8.4, width = 4)
cat("Saved aerf-amef.pdf\n")

cat("\nDone:", SCRIPT_STEM, "\n")
