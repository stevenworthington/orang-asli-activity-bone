###############################################################################
# Estimator comparison for the community-level industrialization -> SOS effect:
# single-level vs multilevel vs two-stage, calibrated on the REAL study design.
#
# The question is whether a separate two-stage estimator is needed at all, or
# whether an ordinary one-stage multilevel model -- the same shape as the rest
# of the paper's analyses -- recovers a cluster-constant exposure effect equally
# well. The existing calibration (calibration-cluster-honest.R) compares only
# the naive individual-level analysis against the two-stage estimator, so it
# cannot settle that. This script adds the multilevel arms.
#
# Five arms, all estimating the SAME scalar (m/s of tibial SOS per index unit):
#   single    lm(y ~ age + sex + index)                       z interval
#             ignores clustering entirely; the Moulton reference point.
#   ml_z      lmer(y ~ age + sex + index + (1 | village))     z interval
#             the standard contextual-effects model, large-sample interval.
#   ml_t      same fit                                        t interval, df = G-2
#             small-cluster df. lmerTest/pbkrtest are not in renv, so we use the
#             between-cluster df convention (G-2), which is also the df the
#             Hartung-Knapp two-stage arm uses -- so the two are compared fairly.
#   ml_gam    gam(y ~ s(age) + sex + s(index, k=4) + s(village, bs="re"))  z
#             the penalized-smooth version, and the configuration the Hodges &
#             Reich spatial-confounding argument concerns: a penalized index
#             smooth competes with the village random intercept for the same
#             between-village variance. Estimate = across-gradient contrast
#             / index range, so it is on the same per-unit scale as the rest.
#   two_stage village-FE adjusted means -> metafor REML meta-regression,
#             Hartung-Knapp t. The frequentist twin of the paper's Bayesian
#             se(sigma = TRUE) estimator, validated against it as a proxy.
#
# DGP, seed, truths and N_SIM are IDENTICAL to calibration-cluster-honest.R, so
# `single` and `two_stage` must reproduce that script's numbers (~0.418 and
# ~0.047 false-positive). That reproduction is the built-in check that the three
# new arms are trustworthy.
#
# Reported per arm: false-positive rate under the null; power, coverage of the
# true slope, and mean estimate (attenuation) under a real effect.
#
# Run: Rscript code/_experiments/calibration-estimator-comparison.R
# Smoke: CALIB_SMOKE=1 (N_SIM = 25).
###############################################################################

suppressMessages({ library(here); library(lme4); library(metafor); library(mgcv) })
set.seed(2138)

out_dir <- here("outputs", "_experiments", "calibration")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


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
idx_range <- diff(range(index_v))
nv <- nlevels(s$village_id)


# ---- DGP parameters from the real data ----

m0 <- lmer(tibia_sos ~ age_c + male + (1 | village_id), data = s)
m1 <- lmer(tibia_sos ~ age_c + male + industrial_index + (1 | village_id), data = s)
mu <- fixef(m0)[["(Intercept)"]]; g_age <- fixef(m0)[["age_c"]]; d_sex <- fixef(m0)[["male"]]
tau_null <- attr(lme4::VarCorr(m0)$village_id, "stddev")[[1]]
tau_eff  <- attr(lme4::VarCorr(m1)$village_id, "stddev")[[1]]
sig      <- sigma(m0)
beta_eff <- fixef(m1)[["industrial_index"]]

N_SIM <- if (nzchar(Sys.getenv("CALIB_SMOKE"))) 25 else 1000
ARMS  <- c("single", "ml_z", "ml_t", "ml_gam", "two_stage")

grid2 <- data.frame(industrial_index = range(index_v))   # endpoints for the GAM contrast


# ---- One replicate: fit all five arms, return estimate + 95% interval each ----

fit_arms <- function(y) {
  s$y <- y
  est <- lo <- hi <- setNames(rep(NA_real_, length(ARMS)), ARMS)

  # single-level, ignores clustering
  co <- summary(lm(y ~ age_c + male + industrial_index, data = s))$coefficients["industrial_index", ]
  est["single"] <- co[1]; lo["single"] <- co[1] - 1.96 * co[2]; hi["single"] <- co[1] + 1.96 * co[2]

  # one-stage multilevel, linear index at level 2
  fm <- tryCatch(lmer(y ~ age_c + male + industrial_index + (1 | village_id), data = s),
                 error = function(e) NULL, warning = function(w) NULL)
  if (!is.null(fm)) {
    b  <- fixef(fm)[["industrial_index"]]
    se <- sqrt(diag(as.matrix(vcov(fm)))[["industrial_index"]])
    tc <- qt(0.975, df = nv - 2)
    est["ml_z"] <- b; lo["ml_z"] <- b - 1.96 * se; hi["ml_z"] <- b + 1.96 * se
    est["ml_t"] <- b; lo["ml_t"] <- b - tc * se;   hi["ml_t"] <- b + tc * se
  }

  # one-stage multilevel with a penalized index smooth (the wash-out arm)
  gm <- tryCatch(gam(y ~ s(age_c, k = 5) + male + s(industrial_index, k = 4) +
                       s(village_id, bs = "re"), data = s, method = "REML"),
                 error = function(e) NULL)
  if (!is.null(gm)) {
    # across-gradient contrast from the lpmatrix, then put on a per-index-unit scale
    # the village random effect is constant across the two rows, so it cancels in the contrast
    Xp <- predict(gm, newdata = transform(grid2, age_c = 0, male = 0,
                                          village_id = factor(levels(s$village_id)[1],
                                                              levels = levels(s$village_id))),
                  type = "lpmatrix")
    cvec <- (Xp[2, ] - Xp[1, ]) / idx_range
    V  <- vcov(gm, unconditional = TRUE)
    b  <- as.numeric(cvec %*% coef(gm))
    se <- sqrt(as.numeric(t(cvec) %*% V %*% cvec))
    est["ml_gam"] <- b; lo["ml_gam"] <- b - 1.96 * se; hi["ml_gam"] <- b + 1.96 * se
  }

  # two-stage: village fixed-effect adjusted means -> Hartung-Knapp meta-regression
  f1 <- lm(y ~ 0 + village_id + age_c + male, data = s)
  c1 <- summary(f1)$coefficients
  vr <- grep("^village_id", rownames(c1))
  rr <- tryCatch(metafor::rma(yi = c1[vr, "Estimate"], sei = c1[vr, "Std. Error"],
                              mods = ~ index_v, method = "REML", test = "knha"),
                 error = function(e) NULL)
  if (!is.null(rr)) {
    est["two_stage"] <- rr$beta[2]; lo["two_stage"] <- rr$ci.lb[2]; hi["two_stage"] <- rr$ci.ub[2]
  }

  c(est = est, lo = lo, hi = hi)
}


# ---- Simulate under a known truth ----

simulate_truth <- function(beta, tau, label) {
  res <- t(replicate(N_SIM, {
    u <- rnorm(nv, 0, tau)
    y <- mu + g_age * s$age_c + d_sex * s$male + u[vid] + beta * index_i + rnorm(nrow(s), 0, sig)
    fit_arms(y)
  }))
  tab <- data.frame(
    arm      = ARMS,
    excl0    = vapply(ARMS, function(a) mean(res[, paste0("lo.", a)] > 0 |
                                             res[, paste0("hi.", a)] < 0, na.rm = TRUE), numeric(1)),
    coverage = vapply(ARMS, function(a) mean(res[, paste0("lo.", a)] <= beta &
                                             res[, paste0("hi.", a)] >= beta, na.rm = TRUE), numeric(1)),
    est_mean = vapply(ARMS, function(a) mean(res[, paste0("est.", a)], na.rm = TRUE), numeric(1)),
    width    = vapply(ARMS, function(a) mean(res[, paste0("hi.", a)] -
                                             res[, paste0("lo.", a)], na.rm = TRUE), numeric(1)),
    n_ok     = vapply(ARMS, function(a) sum(is.finite(res[, paste0("est.", a)])), integer(1)),
    row.names = NULL)

  cat(sprintf("\n=== %s ===\n  true beta = %.3f m/s per index unit (across-gradient %.0f m/s);  village SD %.0f, resid SD %.0f, %d villages, n = %d, N_SIM = %d\n",
              label, beta, beta * idx_range, tau, sig, nv, nrow(s), N_SIM))
  cat(sprintf("  %-10s %10s %10s %10s %10s %6s\n", "arm",
              if (beta == 0) "FP rate" else "power", "coverage", "mean est", "CI width", "n_ok"))
  for (i in seq_len(nrow(tab)))
    cat(sprintf("  %-10s %10.3f %10.3f %10.3f %10.2f %6d\n", tab$arm[i], tab$excl0[i],
                tab$coverage[i], tab$est_mean[i], tab$width[i], tab$n_ok[i]))
  list(raw = res, table = tab)
}

r_null <- simulate_truth(0,        tau_null, "NULL  (no industrialization effect; target FP ~ 0.05)")
r_eff  <- simulate_truth(beta_eff, tau_eff,  "EFFECT  (monotonic SOS decline of the estimated size; target coverage ~ 0.95)")


# ---- Persist ----

summary_tab <- rbind(cbind(truth = "null",   r_null$table),
                     cbind(truth = "effect", r_eff$table))
write.csv(summary_tab, file.path(out_dir, "estimator-comparison-summary.csv"), row.names = FALSE)
saveRDS(list(null = r_null, eff = r_eff,
             params = list(beta_eff = beta_eff, tau_null = tau_null, tau_eff = tau_eff,
                           sig = sig, idx_range = idx_range, nv = nv, n = nrow(s),
                           N_SIM = N_SIM, arms = ARMS)),
        file.path(out_dir, "estimator-comparison.rds"))

cat("\nReproduction check vs calibration-cluster-honest.R (same DGP + seed):\n")
cat(sprintf("  single    null FP = %.3f   (that script's naive arm: 0.418)\n",
            r_null$table$excl0[r_null$table$arm == "single"]))
cat(sprintf("  two_stage null FP = %.3f   (that script's cluster-honest arm: 0.047)\n",
            r_null$table$excl0[r_null$table$arm == "two_stage"]))
cat("\nsaved", file.path(out_dir, "estimator-comparison-summary.csv"), "\n")
