###############################################################################
# Calibration of the cluster-honest industrialization -> SOS inference.
#
# Validates the two-stage cluster-honest estimator (the frequentist twin of
# the village two-stage estimator: village fixed-effect adjusted means -> REML random-effects
# meta-regression on the index, small-sample Hartung-Knapp t) on data that
# matches the REAL design (25 villages, their actual sizes, real age/sex, real
# per-village index), under two KNOWN truths:
#   NULL   - no industrialization effect (all between-village SOS variation is
#            "other stuff"); checks the false-positive rate (should be ~5%).
#   EFFECT - a monotonic SOS decline of the estimated size; checks coverage of
#            the true slope (~95%) and power.
# The NAIVE individual-level regression (ignores clustering) is run alongside to
# show the Moulton over-rejection it suffers.
#
# DGP parameters are estimated from the real SOS data via lme4. The metafor REML
# meta-regression is the frequentist analog of the village two-stage estimator's Bayesian se(sigma=TRUE)
# meta-regression -- same two-stage structure, so its operating characteristics
# track the village two-stage estimator's. Base R + lme4 + metafor (no brms); fast.
#
# Run: script -q /dev/null Rscript code/_experiments/calibration-cluster-honest.R < /dev/null
###############################################################################

suppressMessages({ library(here); library(lme4); library(metafor) })
set.seed(2138)

out_dir <- here("outputs", "_experiments", "calibration")
if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE)
  try(system2("xattr", c("-w", "com.dropbox.ignored", "1", out_dir)), silent = TRUE) }

d <- read.csv(here("data", "processed", "Orang-Asli-pa-vs-bone-60126.csv"))
s <- d[is.finite(d$tibia_sos) & is.finite(d$age_years) & !is.na(d$sex) &
       is.finite(d$industrial_index) & !is.na(d$village_id), ]
s$village_id <- droplevels(factor(s$village_id))
s$age_c <- as.numeric(scale(s$age_years))
s$male  <- as.integer(s$sex == "male")
vid     <- as.integer(s$village_id)
index_v <- as.numeric(tapply(s$industrial_index, s$village_id, function(x) x[1]))  # per village, level order
index_i <- index_v[vid]                                                            # per individual
idx_range <- diff(range(index_v))
nv <- nlevels(s$village_id)

# ---- DGP parameters from the real data ----
m0 <- lmer(tibia_sos ~ age_c + male + (1 | village_id), data = s)
m1 <- lmer(tibia_sos ~ age_c + male + industrial_index + (1 | village_id), data = s)
mu <- fixef(m0)[["(Intercept)"]]; g_age <- fixef(m0)[["age_c"]]; d_sex <- fixef(m0)[["male"]]
tau_null <- attr(lme4::VarCorr(m0)$village_id, "stddev")[[1]]   # all between-village variation
tau_eff  <- attr(lme4::VarCorr(m1)$village_id, "stddev")[[1]]   # residual after the index
sig      <- sigma(m0)
beta_eff <- fixef(m1)[["industrial_index"]]

N_SIM <- 1000

run_methods <- function(y) {
  s$y <- y
  f1 <- lm(y ~ 0 + village_id + age_c + male, data = s)            # stage 1: adjusted village means
  co <- summary(f1)$coefficients
  vr <- grep("^village_id", rownames(co))
  rr <- tryCatch(metafor::rma(yi = co[vr, "Estimate"], sei = co[vr, "Std. Error"],
                              mods = ~ index_v, method = "REML", test = "knha"),
                 error = function(e) NULL)                          # stage 2: cluster-honest meta-reg
  if (is.null(rr)) { ch_excl <- NA; ch_lo <- NA; ch_hi <- NA }
  else { ch_lo <- rr$ci.lb[2]; ch_hi <- rr$ci.ub[2]; ch_excl <- (ch_lo > 0 | ch_hi < 0) }
  fn <- summary(lm(y ~ age_c + male + industrial_index, data = s))$coefficients["industrial_index", ]
  nv_lo <- fn[1] - 1.96 * fn[2]; nv_hi <- fn[1] + 1.96 * fn[2]      # naive individual-level
  c(ch_excl = as.numeric(ch_excl), ch_lo = as.numeric(ch_lo), ch_hi = as.numeric(ch_hi),
    nv_excl = as.numeric((nv_lo > 0) | (nv_hi < 0)))
}

simulate <- function(beta, tau, label) {
  res <- t(replicate(N_SIM, {
    u <- rnorm(nv, 0, tau)
    y <- mu + g_age * s$age_c + d_sex * s$male + u[vid] + beta * index_i + rnorm(nrow(s), 0, sig)
    run_methods(y)
  }))
  cat(sprintf("\n=== %s ===\n  beta = %.2f m/s per index unit  (across-gradient = %.0f m/s);  village SD = %.0f, resid SD = %.0f\n",
              label, beta, beta * idx_range, tau, sig))
  cat(sprintf("  CLUSTER-HONEST  P(95%% CI excludes 0) = %.3f%s\n",
              mean(res[, "ch_excl"], na.rm = TRUE),
              if (beta == 0) "   <- false-positive rate (target ~0.05)" else "   <- power"))
  if (beta != 0)
    cat(sprintf("  CLUSTER-HONEST  coverage of true beta = %.3f   (target ~0.95)\n",
                mean(res[, "ch_lo"] <= beta & res[, "ch_hi"] >= beta, na.rm = TRUE)))
  cat(sprintf("  NAIVE (indiv.)  P(95%% CI excludes 0) = %.3f%s\n",
              mean(res[, "nv_excl"], na.rm = TRUE),
              if (beta == 0) "   <- inflated => the Moulton over-rejection" else "   <- (overconfident)"))
  res
}

r_null <- simulate(0,        tau_null, "NULL  (no industrialization effect)")
r_eff  <- simulate(beta_eff, tau_eff,  "EFFECT  (monotonic SOS decline ~ estimated)")

saveRDS(list(null = r_null, eff = r_eff,
             params = list(beta_eff = beta_eff, tau_null = tau_null, tau_eff = tau_eff,
                           sig = sig, idx_range = idx_range, nv = nv, n = nrow(s), N_SIM = N_SIM)),
        file.path(out_dir, "sos-calibration.rds"))
cat("\nsaved", file.path(out_dir, "sos-calibration.rds"), "\n")
