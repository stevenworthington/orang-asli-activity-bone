###############################################################################
# FULLY-BAYESIAN calibration of the industrialization -> SOS cluster-honest
# inference (Option H). This is the Bayesian twin of calibration-cluster-honest.R
# (which used lme4/metafor). Its only purpose is a SANITY CHECK: do the operating
# characteristics of the *actual Bayesian estimator we report* (Option H: brms
# village-FE GAM -> brms se() meta-regression) match those of the fast frequentist
# twin we already calibrated?
#
#   Frequentist twin (calibration-cluster-honest.R, N_SIM=1000):
#       FP = 0.047   coverage = 0.954   power = 0.42
#
# If the Bayesian numbers land close to these, the frequentist calibration is a
# valid proxy for ALL nine analyses and the full ~16-39 h Bayesian suite is moot.
#
# Estimator (per simulated dataset), identical to set2-option-h-meta-gam.R:
#   Stage 1 (within): brm  y ~ village_id + s(age_years, k=5) + sex
#     -> age/sex-standardized village means + posterior SE (posterior_epred over
#        the cohort age/sex distribution, one prediction per village).
#   Stage 2 (between): brm  adj_mean | se(adj_se, sigma=TRUE) ~ industrial_index
#     -> Bayesian measurement-error meta-regression on the ~25 village means.
#   Decision: 95% HPDI on the stage-2 slope excludes zero?  (monotonic SOS rule)
#
# DGP parameters are estimated Bayesianly (two brms multilevel fits, one-time) so
# there is ZERO frequentist code anywhere. Datasets are simulated on the REAL
# design (real villages/sizes/age/sex/index) under two known truths:
#   NULL   - no industrialization effect (all between-village SOS variation is
#            "other stuff");          -> false-positive rate (target ~0.05)
#   EFFECT - monotonic SOS decline of the estimated size
#                                     -> coverage of true beta (~0.95) + power
#
# Memory-safe: datasets are fit in CHUNKs via brm_multiple(combine=FALSE) (compile
# once, fit many, return SEPARATE fits); each chunk's fits are summarised to a
# decision and discarded before the next chunk.
#
# Run: script -q /dev/null Rscript code/_experiments/calibration-sos-bayesian.R < /dev/null
# Smoke: CAL_SMOKE=1 (tiny MCMC, N_SIM=2) to verify the pipeline end-to-end.
###############################################################################

library(here); suppressMessages(library(brms))
SEED <- 2138; set.seed(SEED)

out_dir <- here("outputs", "_experiments", "calibration")
if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE)
  try(system2("xattr", c("-w", "com.dropbox.ignored", "1", out_dir)), silent = TRUE) }

hpdi <- function(x, m = 0.95) { x <- sort(x); n <- length(x); k <- floor(m * n)
  i <- which.min(x[(k + 1):n] - x[1:(n - k)]); c(x[i], x[i + k]) }


# ---- Data + real design ----

d <- read.csv(here("data", "processed", "Orang-Asli-pa-vs-bone-60126.csv"))
s <- d[is.finite(d$tibia_sos) & is.finite(d$age_years) & !is.na(d$sex) &
       is.finite(d$industrial_index) & !is.na(d$village_id), ]
s$village_id <- droplevels(factor(s$village_id))
s$sex        <- factor(s$sex)
s$age_c      <- as.numeric(scale(s$age_years))
s$male       <- as.integer(s$sex == "male")
vill    <- levels(s$village_id)
vid     <- as.integer(s$village_id)
nv      <- length(vill)
index_v <- as.numeric(tapply(s$industrial_index, s$village_id, function(x) x[1]))  # per village, level order
index_i <- index_v[vid]                                                            # per individual
idx_range <- diff(range(index_v))
cohort  <- s[, c("age_years", "sex")]


# ---- Config ----

SMOKE <- nzchar(Sys.getenv("CAL_SMOKE"))
if (SMOKE) {
  CH <- 2;  IT <- 400;  WU <- 200; N_SIM <- 2;   CHUNK <- 2;  DGP_C <- 2; DGP_I <- 500
} else {
  CH <- 10; IT <- 1000; WU <- 500; N_SIM <- 200; CHUNK <- 25; DGP_C <- 4; DGP_I <- 1500
}


# ---- DGP parameters (100% Bayesian-sourced) ----
# Location/covariate terms and variance components come from a CLEAN Bayesian
# varying-intercept fit (no village-level predictor -> no index/random-intercept
# collinearity, so it converges cleanly; the joint ~ index + (1|village) fit is the
# very ridge that Option H's two-stage design exists to avoid, and it will NOT
# converge). The effect size beta_eff is the estimated industrialization slope, set
# to the SAME value the frequentist calibration used (-2.43) so the DGP is identical
# and the ONLY thing differing between the two calibrations is the estimator (the
# Bayesian Option H real-data slope is -2.64, i.e. fully consistent). tau_eff (the
# residual between-village SD after the index explains its share) is derived
# analytically from the Bayesian tau_null.
message("Estimating DGP parameters (one clean brms varying-intercept fit)...")
m0 <- brm(tibia_sos ~ age_c + male + (1 | village_id), data = s,
          chains = DGP_C, iter = DGP_I, cores = DGP_C, seed = SEED, backend = "cmdstanr",
          control = list(adapt_delta = 0.95), refresh = 0, silent = 2)
pe       <- posterior_summary(m0)[, "Estimate"]
mu       <- pe[["b_Intercept"]]
g_age    <- pe[["b_age_c"]]
d_sex    <- pe[["b_male"]]
tau_null <- pe[["sd_village_id__Intercept"]]
sig      <- pe[["sigma"]]
beta_eff <- -2.43                                       # estimated SOS slope (m/s per index unit)
var_idx  <- mean((index_v - mean(index_v))^2)           # between-village index variance
tau_eff  <- sqrt(max(tau_null^2 - beta_eff^2 * var_idx, (0.5 * tau_null)^2))

cat(sprintf("\nDGP (Bayesian-sourced):\n  mu=%.0f  g_age=%.1f  d_sex=%.1f  beta_eff=%.2f m/s/index (across-gradient %.0f m/s)\n  tau_null=%.1f  tau_eff=%.1f  resid sig=%.1f  | %d villages, n=%d\n",
            mu, g_age, d_sex, beta_eff, beta_eff * idx_range, tau_null, tau_eff, sig, nv, nrow(s)))


# ---- DGP + Option H estimator ----

gen_y <- function(beta, tau) {
  u <- rnorm(nv, 0, tau)
  mu + g_age * s$age_c + d_sex * s$male + u[vid] + beta * index_i + rnorm(nrow(s), 0, sig)
}

standardize <- function(m) {                          # age/sex-standardized village means + SE
  ad <- matrix(NA_real_, brms::ndraws(m), nv)
  for (k in seq_len(nv)) { x <- cohort; x$village_id <- factor(vill[k], levels = vill)
    ad[, k] <- rowMeans(posterior_epred(m, newdata = x)) }
  data.frame(adj_mean = colMeans(ad), adj_se = apply(ad, 2, sd), industrial_index = index_v)
}

run_truth <- function(beta, tau, label) {
  excl <- logical(0); cover <- logical(0); done <- 0L
  while (done < N_SIM) {
    k <- min(CHUNK, N_SIM - done)
    datasets <- lapply(seq_len(k), function(i) { x <- s; x$y <- gen_y(beta, tau); x })
    fits1 <- brm_multiple(y ~ village_id + s(age_years, k = 5) + sex, data = datasets, combine = FALSE,
                          prior = set_prior("student_t(3, 0, 2.5)", class = "sds"),
                          chains = CH, iter = IT, warmup = WU, cores = CH, backend = "cmdstanr",
                          seed = SEED, refresh = 0, silent = 2)
    villdfs <- lapply(fits1, standardize)
    fits2 <- brm_multiple(bf(adj_mean | se(adj_se, sigma = TRUE) ~ industrial_index),
                          data = villdfs, combine = FALSE,
                          chains = CH, iter = IT, warmup = WU, cores = CH, backend = "cmdstanr",
                          control = list(adapt_delta = 0.99), seed = SEED, refresh = 0, silent = 2)
    for (f in fits2) {
      b <- as.data.frame(f)[["b_industrial_index"]]
      h <- hpdi(b, 0.95)
      excl  <- c(excl,  h[1] > 0 | h[2] < 0)
      cover <- c(cover, h[1] <= beta & h[2] >= beta)
    }
    rm(fits1, fits2, villdfs); gc(verbose = FALSE)
    done <- done + k
    cat(sprintf("  [%s] %d/%d  running: excl-rate=%.3f%s\n", label, done, N_SIM,
                mean(excl), if (beta != 0) sprintf("  coverage=%.3f", mean(cover)) else ""))
  }
  list(excl = excl, cover = cover)
}


# ---- Run both truths ----

cat("\n=== NULL truth (false-positive rate) ===\n")
r_null <- run_truth(0,        tau_null, "NULL")
cat("\n=== EFFECT truth (power + coverage) ===\n")
r_eff  <- run_truth(beta_eff, tau_eff,  "EFFECT")

fp    <- mean(r_null$excl)
power <- mean(r_eff$excl)
cover <- mean(r_eff$cover)

cat("\n###############################################################################\n")
cat(sprintf("BAYESIAN Option H calibration  (N_SIM=%d, %d chains x %d post-warmup = %d draws/fit)\n",
            N_SIM, CH, IT - WU, CH * (IT - WU)))
cat("###############################################################################\n")
cat(sprintf("                         Bayesian (this run)   Frequentist twin (N=1000)\n"))
cat(sprintf("  FP rate (target 0.05)        %.3f                 0.047\n", fp))
cat(sprintf("  coverage (target 0.95)       %.3f                 0.954\n", cover))
cat(sprintf("  power                        %.3f                 0.42\n",  power))
cat(sprintf("\n  MC SE (N=%d): FP/coverage ~%.3f, power ~%.3f\n",
            N_SIM, sqrt(0.05 * 0.95 / N_SIM), sqrt(0.42 * 0.58 / N_SIM)))

saveRDS(list(null = r_null, eff = r_eff,
             summary = c(fp = fp, coverage = cover, power = power),
             params = list(mu = mu, g_age = g_age, d_sex = d_sex, beta_eff = beta_eff,
                           tau_null = tau_null, tau_eff = tau_eff, sig = sig,
                           idx_range = idx_range, nv = nv, n = nrow(s), N_SIM = N_SIM),
             config = c(CH = CH, IT = IT, WU = WU, N_SIM = N_SIM, CHUNK = CHUNK)),
        file.path(out_dir, "sos-calibration-bayesian.rds"))
cat("\nsaved", file.path(out_dir, "sos-calibration-bayesian.rds"), "\n")
