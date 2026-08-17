###############################################################################
# Calibration of the BAYESIAN one-stage hierarchical GAM for the community-level
# industrialization effect.
#
# WHY. This is the estimator the industrialization results are reported from,
# and its operating characteristics at ~25 clusters cannot be assumed. The
# frequentist small-sample recommendations for few clusters (Satterthwaite /
# Kenward-Roger / between-within degrees of freedom; Leyrat et al. 2018,
# Donald & Lang 2007) have no literal Bayesian counterpart.
# The Bayesian analogue is to put a prior on the between-cluster SD and integrate
# over it, so the heavier tails arise by marginalization rather than correction.
# Whether that actually delivers nominal behaviour at 25 clusters is an empirical
# question, and this script answers it.
#
# DESIGN. The real study design throughout: real SOS analytic sample, 25
# communities, actual sizes 6-79, real age/sex/index, variance components taken
# from the fitted data.
#
#   MODEL   brm(y ~ s(age_years, k=5) + sex + s(industrial_index, k=4) +
#                   (1 | village_id))
#   ESTIMAND  linear-projection slope of the AERF per index unit, by
#             G-computation over the cohort's age/sex distribution, marginal of
#             community (re_formula = NA). See the estimand note further down.
#   DECISION  95% HPDI of that slope excludes zero.
#
# TWO CONDITIONS, both at the (kappa, community-SD prior mean) cell requested.
#   null    FP rate                            -- target 0.05
#   effect  power + coverage of the true slope -- coverage target 0.95
#
# The outcome is modelled on the production scale for its own spec (SOS /200,
# steps /5,000, ENMO /10) so the standing sds prior is on the right scale;
# slopes are converted back to natural units per index unit.
# Each prior setting compiles once; replicates reuse it via update(recompile = FALSE).
#
# Run: Rscript code/_experiments/calibration-onestage-production.R
# Env: CAL_OUTCOME (sos|steps|enmo), CAL_KAPPA, CAL_SDMEAN, ONESTAGE_N (default
#      120), ONESTAGE_SMOKE=1. The defaults are the production cell on SOS, so
#      running the script bare calibrates what the registry reports.
###############################################################################

suppressMessages({ library(here); library(lme4); library(brms) })
set.seed(2138)

out_dir <- here("outputs", "_experiments", "calibration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hpdi <- function(x, m = 0.95) { x <- sort(x); n <- length(x); k <- floor(m * n)
  i <- which.min(x[(k + 1):n] - x[1:(n - k)]); c(x[i], x[i + k]) }

# The (kappa, community-SD prior mean) cell being calibrated. Both are swept
# rather than fixed: a community-SD prior sitting below the observed
# between-community SD shrinks the variance component and costs coverage, and a
# sweep is the only way to see where the acceptable region is.
# Defaults are the SELECTED PRODUCTION CELL, so running this script bare
# reproduces the calibration the registry quotes.
CAL_KAPPA  <- as.numeric(Sys.getenv("CAL_KAPPA", "0.5"))
CAL_SDMEAN <- as.numeric(Sys.getenv("CAL_SDMEAN", "0.33"))
# Which outcome to simulate from. All three are calibrated separately: the two
# activity outcomes sit on 29 communities against SOS's 25, with different
# residual scales and different between-community SDs, so SOS's operating
# characteristics do not transfer to them.
CAL_OUTCOME <- Sys.getenv("CAL_OUTCOME", "sos")
stopifnot(CAL_OUTCOME %in% c("sos", "steps", "enmo"))
OUT <- list(
  sos   = list(col = "tibia_sos",                scale = 200,  unit = "m/s"),
  steps = list(col = "ad_tot_step_count_0_24hr", scale = 5000, unit = "steps"),
  enmo  = list(col = "ad_mean_enmo_mg_0_24hr",   scale = 10,   unit = "mg")
)[[CAL_OUTCOME]]
SMOKE   <- nzchar(Sys.getenv("ONESTAGE_SMOKE"))
N_MAIN  <- if (SMOKE) 2 else as.integer(Sys.getenv("ONESTAGE_N", "120"))
CFG     <- if (SMOKE) list(w = 200, i = 600, c = 2) else list(w = 400, i = 1250, c = 4)
SCALE   <- OUT$scale   # production rescaling for this outcome


# ---- Real design ----

d <- read.csv(here("data", "processed", "orang-asli-pa-bone-data.csv"))
s <- d[is.finite(d[[OUT$col]]) & is.finite(d$age_years) & !is.na(d$sex) &
       is.finite(d$industrial_index) & !is.na(d$village_id), ]
s$village_id <- droplevels(factor(s$village_id))
s$age_c <- as.numeric(scale(s$age_years))
s$male  <- as.integer(s$sex == "male")
vid     <- as.integer(s$village_id)
index_v <- as.numeric(tapply(s$industrial_index, s$village_id, function(x) x[1]))
index_i <- index_v[vid]
nv      <- nlevels(s$village_id)
rng     <- range(index_v); span <- diff(rng)

s$.y <- s[[OUT$col]]
m0 <- lmer(.y ~ age_c + male + (1 | village_id), data = s)
m1 <- lmer(.y ~ age_c + male + industrial_index + (1 | village_id), data = s)
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

# ESTIMAND. The LINEAR PROJECTION of the AERF over the grid -- what production
# reports -- and deliberately NOT an endpoint contrast (the difference in the
# AERF between the two extreme index values, divided by the span). Those are not
# the same functional in practice: the fitted smooth wiggles even when the truth
# is linear, and the endpoint reads the two points where communities are
# sparsest, so it is the noisier of the two. On the production draws the two can
# differ by enough to disagree about whether the interval excludes zero. A
# calibration must characterise the quantity actually reported, so both
# harnesses compute the linear projection over a 21-point grid, which reproduces
# the full 41/51-point projection closely.
#
# Linear-projection slope per index unit, by G-computation marginal of community.
GRID <- seq(rng[1], rng[2], length.out = 21)
GXC  <- GRID - mean(GRID)
slope_draws <- function(fit, dd) {
  cohort <- dd[, c("age_years", "sex")]
  M <- vapply(GRID, function(v) {
         nd <- cohort; nd$industrial_index <- v
         rowMeans(posterior_epred(fit, newdata = nd, re_formula = NA))
       }, numeric(brms::ndraws(fit)))
  as.numeric(M %*% GXC / sum(GXC^2)) * SCALE
}

# The PRODUCTION prior set, built the same way production builds it. That is the
# point of the harness: a small prior mean on the community SD can sit below the
# between-community SD the data actually show. A prior that tight on a variance
# component can shrink community heterogeneity and under-cover, and calibration
# is the only thing that detects it -- which is why CAL_SDMEAN is swept rather
# than fixed.
#
# Priors are built ONCE from the first simulated dataset and then held fixed,
# because brms writes prior values into the Stan code and re-deriving them per
# replicate would force a recompile on all 120. Every replicate shares the same
# design and DGP, so sd_lp varies negligibly; production derives them from its
# own analytic data in exactly this way.
production_priors <- function(form, d, kappa = CAL_KAPPA, sd_mean = CAL_SDMEAN) {
  y <- d$y; sd_lp <- sd(y); mu_lp <- mean(y)
  sdat <- brms::standata(form, data = d); colsd <- c()
  for (m in c("X", "Xs")) if (!is.null(sdat[[m]])) {
    M <- as.matrix(sdat[[m]]); v <- apply(M, 2, sd); names(v) <- colnames(M)
    colsd <- c(colsd, v) }
  colsd <- colsd[is.finite(colsd) & colsd > 0]
  gp <- as.data.frame(brms::get_prior(form, data = d))
  pri <- NULL
  for (cf in gp$coef[gp$class == "b" & gp$coef != ""]) {
    sx <- colsd[[cf]]; if (is.null(sx) || !is.finite(sx)) next
    pp <- set_prior(sprintf("normal(0, %.8g)", kappa * sd_lp / sx), class = "b", coef = cf)
    pri <- if (is.null(pri)) pp else pri + pp }
  pri +
    set_prior(sprintf("normal(%.8g, %.8g)", mu_lp, 0.5 * sd_lp), class = "Intercept") +
    set_prior("exponential(2)", class = "sds") +
    set_prior(sprintf("exponential(%.8g)", 1 / (sd_mean * sd_lp)), class = "sd") +
    set_prior(sprintf("exponential(%.8g)", 1 / sd_lp), class = "sigma")
}

run_condition <- function(label, beta, tau, n_sim) {

  d0  <- sim_frame(beta, tau)
  pri <- production_priors(FORM, d0)
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
  cat(sprintf("\n=== %s | %s ===\n  true beta = %.4g %s per index unit, community SD %.4g, %d communities, N_SIM = %d\n",
              label, CAL_OUTCOME, beta, OUT$unit, tau, nv, n_sim))
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
  run_condition(sprintf("null   (kappa %g, sd mean %g x sd_lp)", CAL_KAPPA, CAL_SDMEAN),
                0, tau_null, N_MAIN),
  run_condition(sprintf("effect (kappa %g, sd mean %g x sd_lp)", CAL_KAPPA, CAL_SDMEAN),
                beta_eff, tau_eff, N_MAIN))

tab <- do.call(rbind, rows)
tab$kappa <- CAL_KAPPA; tab$sd_mean <- CAL_SDMEAN; tab$outcome <- CAL_OUTCOME
tab$n_communities <- nv; tab$n_obs <- nrow(s)
write.csv(tab, file.path(out_dir, sprintf("onestage-production-calibration_%s_k%g_sd%g.csv", CAL_OUTCOME, CAL_KAPPA, CAL_SDMEAN)), row.names = FALSE)

cat("\nReference points:\n")
cat("  nominal target:                       FP 0.050, coverage 0.950\n")
# Read from the file that generates it rather than restated. Hard-coded twin
# rates go stale silently as soon as the twin changes, and nothing here would
# catch it.
val_path <- here("outputs", "_experiments", "power-curves", "onestage-proxy-validation.csv")
if (file.exists(val_path)) {
  v <- utils::read.csv(val_path)
  cat(sprintf("  frequentist twin (SOS design):        FP %.3f, power %.3f\n",
              v$proxy[v$condition == "null"], v$proxy[v$condition == "effect"]))
} else {
  cat("  frequentist twin:                     run power-curves-industrialization-onestage.R\n")
}
cat("  naive individual-level, no clustering: far above nominal -- the reason\n")
cat("    inference cannot happen at the individual level. Set CAL_OUTCOME and\n")
cat("    drop the community term to measure it on this design.\n")
cat("\nsaved", file.path(out_dir, sprintf("onestage-production-calibration_%s_k%g_sd%g.csv", CAL_OUTCOME, CAL_KAPPA, CAL_SDMEAN)), "\n")
