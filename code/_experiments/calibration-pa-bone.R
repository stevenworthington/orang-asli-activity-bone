###############################################################################
# Operating-characteristic calibration for the PA -> bone specifications.
#
# WHY THIS EXISTS. The PA results are NULLS reported as magnitude bounds --
# "95% of the effect is below X" -- and that claim rests entirely on the
# interval having its nominal coverage. If the intervals under-cover, every
# reported U95 is overstated and the evidence-of-absence framing is stronger
# than the data support. This is not hypothetical: in the industrialization
# calibration a prior placed too tightly on a variance component measurably
# degraded coverage, which is what motivated measuring it here rather than
# assuming it.
#
# Two features of these specs make it worth measuring rather than assuming:
#   * ~24 community indicators carry FLAT priors, deliberately, so the exposure
#     effect is identified from within-community contrasts only. That is a lot
#     of unpenalized nuisance parameters, and their effect on coverage is not
#     obvious in either direction.
#   * all three likelihoods are non-Gaussian, and two of them are new here and
#     have never been calibrated: Student-t with an estimated nu, on the SOS arm
#     at n ~ 650, and Gamma(log) on CTX-1 at n ~ 265. Osteocalcin's lognormal is
#     unchanged and acts as the regression test.
#
# DESIGN. One spec per likelihood family, at the production settings:
#   sos-steps    Student-t, tibia_sos/200, n ~ 650
#   ctx-steps    Gamma(log), ctx1_ng_ml,   n ~ 265
#   osteo-steps  lognormal,  osteo/10000,  n ~ 265
#
# The DGP is fitted to the real data on each spec's own linear-predictor scale
# and carries the real design: real ages, sex, covariates and community
# identities. Community offsets are held FIXED across replicates, because the
# production model treats them as fixed effects -- redrawing them each replicate
# would calibrate a different estimator.
#
# ESTIMAND. The LINEAR-PROJECTION slope of the AERF by G-computation over a
# 21-point exposure grid, scaled to the reported contrast of 5,000 steps -- the
# quantity the manuscript reports. See the estimand note further down for why
# the endpoint contrast is not an adequate proxy for it, and why 21 points is
# enough (it reproduces the full 51-point projection to within 1.7%).
#
# Run: CAL_SPEC=sos-steps CAL_COND=null Rscript \
#        code/_experiments/calibration-pa-bone.R
# Env: CAL_SPEC, CAL_COND (null|effect), PA_CAL_N (200), PA_CAL_SMOKE=1
###############################################################################

library(here)
source(here("code", "_startup", "init.R"))
suppressMessages({ library(dplyr); library(brms) })
set.seed(2138)

CAL_SPEC <- Sys.getenv("CAL_SPEC", "sos-steps")
CAL_COND <- Sys.getenv("CAL_COND", "null")
stopifnot(CAL_SPEC %in% c("sos-steps", "ctx-steps", "osteo-steps"),
          CAL_COND %in% c("null", "effect"))
SMOKE  <- nzchar(Sys.getenv("PA_CAL_SMOKE"))
N_SIM  <- if (SMOKE) 2 else as.integer(Sys.getenv("PA_CAL_N", "200"))
CFG    <- if (SMOKE) list(w = 200, i = 600, c = 2) else list(w = 400, i = 1200, c = 4)
KAPPA  <- KAPPA_PA          # from the registry, not a second copy
# From the registry: since data.R scales steps by 5,000, one unit of
# ad_steps_5k IS the reported contrast, so this is 1.
CONTRAST_UNITS <- model_templates[[CAL_SPEC]]$contrast_units

out_dir <- here("outputs", "_experiments", "calibration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hpdi <- function(x, m = 0.95) { x <- sort(x); n <- length(x); k <- floor(m * n)
  i <- which.min(x[(k + 1):n] - x[1:(n - k)]); c(x[i], x[i + k]) }

SPECS <- list(
  "sos-steps"   = list(out = "tibia_sos_200", raw = "tibia_sos",         sf = 200,
                       fam = "student",   unit = "m/s"),
  "ctx-steps"   = list(out = "ctx1_ng_ml",    raw = "ctx1_ng_ml",        sf = 1,
                       fam = "gamma",     unit = "ng/mL"),
  "osteo-steps" = list(out = "osteocalcin_pg_ml_10k", raw = "osteocalcin_pg_ml", sf = 10000,
                       fam = "lognormal", unit = "pg/mL"))
S  <- SPECS[[CAL_SPEC]]
EX <- model_templates[[CAL_SPEC]]$exposure   # ad_steps_5k
COVS <- c("sex", "pregnant_or_breastfeeding_n_y_0_1", "smoking_binary_n_y_0_1",
          "alcohol_binary_n_y_0_1", "functional_status_n_y_0_1")


# ---- Real design ----------------------------------------------------------

# The scaled columns come from prep_dat(), the same derivation production
# uses, so the harness cannot calibrate a differently-scaled model.
keep <- c(S$out, EX, "age_years", "village_id", COVS)
s <- dat[stats::complete.cases(dat[, keep]), keep]
s$village_id <- droplevels(factor(s$village_id))
is_log <- S$fam %in% c("gamma", "lognormal")

# DGP on the scale the model is linear on, with community identity as a FIXED
# effect -- matching the production specification.
s$.lp <- if (is_log) log(s[[S$out]]) else s[[S$out]]
dgp_form <- as.formula(paste(".lp ~ age_years + I(age_years^2) +",
                             paste(COVS, collapse = " + "), "+ village_id +", EX))
m_dgp  <- lm(dgp_form, data = s)
beta_hat <- unname(coef(m_dgp)[[EX]])
sig_lp   <- sigma(m_dgp)
eta0     <- unname(fitted(m_dgp) - coef(m_dgp)[[EX]] * s[[EX]])   # everything but exposure
BETA     <- if (CAL_COND == "null") 0 else beta_hat

# Heavy-tail and shape parameters for the two non-Gaussian residual models,
# taken from the real fits rather than assumed.
NU    <- 6                                  # tibial SOS, matches nu ~ 6.0-6.3
SHAPE <- 1 / (sd(s[[S$out]]) / mean(s[[S$out]]))^2   # Gamma: shape = 1/CV^2

sim_frame <- function(beta) {
  eta <- eta0 + beta * s[[EX]]
  y <- switch(S$fam,
    student   = eta + sig_lp * rt(nrow(s), df = NU) / sqrt(NU / (NU - 2)),
    lognormal = exp(eta + rnorm(nrow(s), 0, sig_lp)),
    gamma     = rgamma(nrow(s), shape = SHAPE, rate = SHAPE / exp(eta)))
  d <- s; d[[S$out]] <- y; d
}

FAM <- switch(S$fam, student = student(), lognormal = lognormal(),
              gamma = Gamma(link = "log"))
FORM <- bf(as.formula(sprintf(
  "%s ~ t2(age_years, %s, k = c(5, 5)) + %s + village_id",
  S$out, EX, paste(COVS, collapse = " + "))), family = FAM)

grid_rng <- quantile(s[[EX]], c(0.01, 0.99), na.rm = TRUE)
span     <- diff(grid_rng)

# ESTIMAND. The LINEAR PROJECTION of the AERF over the exposure grid -- the
# quantity the paper reports -- and deliberately NOT an endpoint contrast (the
# AERF's difference between the two extreme exposure values, over the span).
# Those are not the same functional in practice: the fitted smooth wiggles even
# when the truth is linear, and the endpoint reads only the two points where
# the data are sparsest, so it is the noisier of the two. A calibration has to
# characterise the quantity actually reported, or it describes a different
# estimator. 21 grid points reproduce the full 51-point projection to within
# 1.7%, at a twenty-fifth of the cost.
#
# The true value under this DGP is linear in the exposure on the
# LINEAR-PREDICTOR scale, so under a log link the response-scale slope is not
# constant. It is therefore evaluated exactly as the estimate is, by
# G-computation on the DGP itself.
GRID <- seq(grid_rng[1], grid_rng[2], length.out = 21)
GXC  <- GRID - mean(GRID)
true_slope <- {
  bt <- if (is_log) function(z) exp(z) else function(z) z
  m  <- vapply(GRID, function(v) mean(bt(eta0 + BETA * v)), numeric(1))
  sum(GXC * (m - mean(m))) / sum(GXC^2) * CONTRAST_UNITS * S$sf
}

slope_draws <- function(fit, d) {
  M <- vapply(GRID, function(v) {
         nd <- d; nd[[EX]] <- v
         rowMeans(posterior_epred(fit, newdata = nd))
       }, numeric(brms::ndraws(fit)))
  as.numeric(M %*% GXC / sum(GXC^2)) * CONTRAST_UNITS * S$sf
}

# Production priors: kappa-autoscaled on interpretable coefficients, community
# indicators FLAT. Built once from the first replicate and held fixed, because
# brms writes prior values into the Stan code and re-deriving them per replicate
# would force 200 recompiles.
# Production build_priors(), not a second copy: a calibration harness has to
# exercise the production prior code, or it characterises something other than
# the model being reported. It takes a spec rather than a formula, so a minimal
# spec-alike is assembled; the fields it reads are bf, outcome, prior_kappa and
# (absent here) prior_sd_community.
cal_spec_for_priors <- list(bf = FORM, outcome = S$out, prior_kappa = KAPPA)

message(sprintf("=== %s | %s | n %d, %d communities | true slope %.5g %s per 5,000 steps ===",
                CAL_SPEC, CAL_COND, nrow(s), nlevels(s$village_id), true_slope, S$unit))

d0  <- sim_frame(BETA)
pri <- build_priors(cal_spec_for_priors, d0)
fit <- brm(FORM, data = d0, prior = pri, warmup = CFG$w, iter = CFG$i,
           chains = CFG$c, cores = CFG$c, seed = 2138, backend = "cmdstanr",
           control = list(adapt_delta = 0.95, max_treedepth = 12),
           refresh = 0, silent = 2)

res <- matrix(NA_real_, N_SIM, 4, dimnames = list(NULL, c("est", "lo", "hi", "u95")))
for (r in seq_len(N_SIM)) {
  dd <- if (r == 1) d0 else sim_frame(BETA)
  f  <- if (r == 1) fit else update(fit, newdata = dd, recompile = FALSE,
                                    refresh = 0, silent = 2)
  b  <- slope_draws(f, dd)
  h  <- hpdi(b)
  res[r, ] <- c(mean(b), h[1], h[2], unname(quantile(abs(b), 0.95)))
  if (r %% 10 == 0 || r == N_SIM) message(sprintf("  %s/%s %d/%d", CAL_SPEC, CAL_COND, r, N_SIM))
}

excl <- mean(res[, "lo"] > 0 | res[, "hi"] < 0, na.rm = TRUE)
cov  <- mean(res[, "lo"] <= true_slope & res[, "hi"] >= true_slope, na.rm = TRUE)
tab <- data.frame(
  spec = CAL_SPEC, condition = CAL_COND, family = S$fam, unit = S$unit,
  n = nrow(s), n_communities = nlevels(s$village_id), n_sim = N_SIM,
  true_slope = true_slope, excl0 = excl, coverage = cov,
  est_mean = mean(res[, "est"], na.rm = TRUE),
  width = mean(res[, "hi"] - res[, "lo"], na.rm = TRUE),
  u95_mean = mean(res[, "u95"], na.rm = TRUE),
  u95_covers = mean(res[, "u95"] >= abs(true_slope), na.rm = TRUE))

cat(sprintf("\n=== %s | %s ===\n", CAL_SPEC, CAL_COND))
cat(sprintf("  true slope        = %.5g %s per 5,000 steps\n", true_slope, S$unit))
cat(sprintf("  P(HPDI excludes 0)= %.3f%s\n", excl,
            if (CAL_COND == "null") "   <- false-positive rate (target ~0.05)" else "   <- power"))
cat(sprintf("  coverage          = %.3f   (target 0.95)\n", cov))
cat(sprintf("  mean estimate     = %.5g\n", tab$est_mean))
cat(sprintf("  mean HPDI width   = %.5g\n", tab$width))
cat(sprintf("  mean U95 bound    = %.5g\n", tab$u95_mean))
cat(sprintf("  U95 >= |true|     = %.3f   <- the magnitude-bound claim; should be ~1\n",
            tab$u95_covers))

write.csv(tab, file.path(out_dir, sprintf("pa-bone-calibration_%s_%s.csv", CAL_SPEC, CAL_COND)),
          row.names = FALSE)
cat("\nsaved", file.path(out_dir, sprintf("pa-bone-calibration_%s_%s.csv", CAL_SPEC, CAL_COND)), "\n")
