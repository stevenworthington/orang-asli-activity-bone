###############################################################################
# Calibration of the BAYESIAN one-stage hierarchical GAM for the community-level
# industrialization effect.
#
# WHY. Two estimators have already been calibrated: the Bayesian TWO-STAGE
# estimator (FP 0.040, coverage 0.960) and the FREQUENTIST one-stage twins
# (ml_t FP 0.049, ml_gam FP 0.049, coverage 0.954). The
# Bayesian ONE-STAGE model has never had its operating characteristics measured,
# so it cannot yet be described as independently validated -- only as agreeing
# with the two-stage estimator on the real data.
#
# This matters because the frequentist small-sample recommendations for few
# clusters (Satterthwaite / Kenward-Roger / between-within degrees of freedom;
# Leyrat et al. 2018, Donald & Lang 2007) have no literal Bayesian counterpart.
# The Bayesian analogue is to put a prior on the between-cluster SD and integrate
# over it, so the heavier tails arise by marginalization rather than correction.
# Whether that actually delivers nominal behaviour at 25 clusters is an empirical
# question, and this script answers it.
#
# DESIGN. Identical DGP, seed and truths to calibration-estimator-comparison.R:
# real SOS analytic sample, 25 villages, actual sizes 6-79, real age/sex/index,
# variance components from the fitted data.
#
#   MODEL   brm(y ~ s(age_years, k=5) + sex + s(industrial_index, k=4) +
#                   (1 | village_id))
#   ESTIMAND  across-gradient slope per index unit, by G-computation over the
#             cohort's age/sex distribution, marginal of village (re_formula = NA).
#             Same scalar the frequentist ml_gam arm reported, so the two are
#             directly comparable.
#   DECISION  95% HPDI of that slope excludes zero.
#
# THREE CONDITIONS.
#   null_default    FP rate, brms default (data-scaled) prior on the village SD
#   effect_default  power + coverage of the true slope, same prior
#   null_tight      FP rate under a DELIBERATELY BAD prior -- half-normal(0, 0.01)
#                   on the village SD, i.e. 10 m/s when the truth is ~70 m/s.
#                   This is the demonstration that in a Bayesian framework the
#                   prior on tau does the work a frequentist df correction does
#                   (Gelman 2006): squeeze tau toward zero and the model should
#                   revert toward the naive individual-level behaviour. If
#                   null_tight is badly anti-conservative while null_default is
#                   nominal, the prior is doing that work, and the choice is a
#                   substantive modelling decision rather than an incidental
#                   default.
#
# Outcome is modelled on the project's /1000 scale so the standing sds prior is
# on the right scale; slopes are converted back to m/s per index unit.
# Each prior setting compiles once; replicates reuse it via update(recompile = FALSE).
#
# Run: Rscript code/_experiments/calibration-onestage-bayes.R
# Env: ONESTAGE_N (default 120), ONESTAGE_N_TIGHT (default 60), ONESTAGE_SMOKE=1.
###############################################################################

suppressMessages({ library(here); library(lme4); library(brms) })
set.seed(2138)

out_dir <- here("outputs", "_experiments", "calibration")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

hpdi <- function(x, m = 0.95) { x <- sort(x); n <- length(x); k <- floor(m * n)
  i <- which.min(x[(k + 1):n] - x[1:(n - k)]); c(x[i], x[i + k]) }

SMOKE   <- nzchar(Sys.getenv("ONESTAGE_SMOKE"))
N_MAIN  <- if (SMOKE) 2 else as.integer(Sys.getenv("ONESTAGE_N", "120"))
N_TIGHT <- if (SMOKE) 2 else as.integer(Sys.getenv("ONESTAGE_N_TIGHT", "60"))
CFG     <- if (SMOKE) list(w = 200, i = 600, c = 2) else list(w = 400, i = 1250, c = 4)
SCALE   <- 1000


# ---- Real design ----

d <- read.csv(here("data", "processed", "Orang-Asli-pa-vs-bone-60126.csv"))
s <- d[is.finite(d$tibia_sos) & is.finite(d$age_years) & !is.na(d$sex) &
       is.finite(d$industrial_index) & !is.na(d$village_id), ]
s$village_id <- droplevels(factor(s$village_id))
s$age_c <- as.numeric(scale(s$age_years))
s$male  <- as.integer(s$sex == "male")
vid     <- as.integer(s$village_id)
index_v <- as.numeric(tapply(s$industrial_index, s$village_id, function(x) x[1]))
index_i <- index_v[vid]
nv      <- nlevels(s$village_id)
rng     <- range(index_v); span <- diff(rng)

m0 <- lmer(tibia_sos ~ age_c + male + (1 | village_id), data = s)
m1 <- lmer(tibia_sos ~ age_c + male + industrial_index + (1 | village_id), data = s)
mu <- fixef(m0)[["(Intercept)"]]; g_age <- fixef(m0)[["age_c"]]; d_sex <- fixef(m0)[["male"]]
tau_null <- attr(lme4::VarCorr(m0)$village_id, "stddev")[[1]]
tau_eff  <- attr(lme4::VarCorr(m1)$village_id, "stddev")[[1]]
sig      <- sigma(m0)
beta_eff <- fixef(m1)[["industrial_index"]]

sim_frame <- function(beta, tau) {
  u <- rnorm(nv, 0, tau)
  y <- mu + g_age * s$age_c + d_sex * s$male + u[vid] + beta * index_i +
       rnorm(nrow(s), 0, sig)
  data.frame(y = y / SCALE, age_years = s$age_years, sex = s$sex,
             village_id = s$village_id, industrial_index = s$industrial_index)
}

FORM <- bf(y ~ s(age_years, k = 5) + sex + s(industrial_index, k = 4) +
             (1 | village_id), family = gaussian())

# across-gradient slope per index unit, in m/s, by G-computation marginal of village
slope_draws <- function(fit, dd) {
  cohort <- dd[, c("age_years", "sex")]
  lo <- cohort; lo$industrial_index <- rng[1]
  hi <- cohort; hi$industrial_index <- rng[2]
  (rowMeans(posterior_epred(fit, newdata = hi, re_formula = NA)) -
   rowMeans(posterior_epred(fit, newdata = lo, re_formula = NA))) * SCALE / span
}

run_condition <- function(label, beta, tau, n_sim, extra_prior = NULL) {
  pri <- set_prior("student_t(3, 0, 2.5)", class = "sds")
  if (!is.null(extra_prior)) pri <- pri + extra_prior

  d0  <- sim_frame(beta, tau)
  fit <- brm(FORM, data = d0, prior = pri, warmup = CFG$w, iter = CFG$i,
             chains = CFG$c, cores = CFG$c, seed = 2138, backend = "cmdstanr",
             control = list(adapt_delta = 0.95, max_treedepth = 12),
             refresh = 0, silent = 2)

  res <- matrix(NA_real_, nrow = n_sim, ncol = 3,
                dimnames = list(NULL, c("est", "lo", "hi")))
  for (r in seq_len(n_sim)) {
    dd <- if (r == 1) d0 else sim_frame(beta, tau)
    f  <- if (r == 1) fit else update(fit, newdata = dd, recompile = FALSE,
                                      refresh = 0, silent = 2)
    b  <- slope_draws(f, dd)
    h  <- hpdi(b)
    res[r, ] <- c(mean(b), h[1], h[2])
    if (r %% 10 == 0 || r == n_sim) message(sprintf("  %s %d/%d", label, r, n_sim))
  }

  excl <- mean(res[, "lo"] > 0 | res[, "hi"] < 0, na.rm = TRUE)
  cov  <- mean(res[, "lo"] <= beta & res[, "hi"] >= beta, na.rm = TRUE)
  cat(sprintf("\n=== %s ===\n  true beta = %.3f m/s per index unit, village SD %.0f, %d villages, N_SIM = %d\n",
              label, beta, tau, nv, n_sim))
  cat(sprintf("  P(95%% HPDI excludes 0) = %.3f%s\n", excl,
              if (beta == 0) "   <- false-positive rate (target ~0.05)" else "   <- power"))
  cat(sprintf("  coverage of true beta   = %.3f\n", cov))
  cat(sprintf("  mean estimate           = %.3f\n", mean(res[, "est"], na.rm = TRUE)))
  cat(sprintf("  mean HPDI width         = %.3f\n",
              mean(res[, "hi"] - res[, "lo"], na.rm = TRUE)))
  data.frame(condition = label, true_beta = beta, tau = tau, n_sim = n_sim,
             excl0 = excl, coverage = cov,
             est_mean = mean(res[, "est"], na.rm = TRUE),
             width = mean(res[, "hi"] - res[, "lo"], na.rm = TRUE))
}

rows <- list(
  run_condition("null_default   (brms default prior on village SD)", 0, tau_null, N_MAIN),
  run_condition("effect_default (brms default prior on village SD)", beta_eff, tau_eff, N_MAIN),
  run_condition("null_tight     (half-normal(0, 10 m/s) on village SD)", 0, tau_null, N_TIGHT,
                extra_prior = set_prior("normal(0, 0.01)", class = "sd", lb = 0)))

tab <- do.call(rbind, rows)
write.csv(tab, file.path(out_dir, "onestage-bayes-calibration.csv"), row.names = FALSE)

cat("\nReference points:\n")
cat("  Bayesian two-stage:      FP 0.040, coverage 0.960\n")
cat("  frequentist ml_gam:      FP 0.049, coverage 0.954, power 0.409\n")
cat("  frequentist two_stage:   FP 0.047, coverage 0.954, power 0.419\n")
cat("  naive individual-level:  FP 0.418\n")
cat("\nsaved", file.path(out_dir, "onestage-bayes-calibration.csv"), "\n")
