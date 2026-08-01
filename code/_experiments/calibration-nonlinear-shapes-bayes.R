###############################################################################
# Bayesian confirmation of the non-linear-shape calibration.
#
# calibration-nonlinear-shapes.R runs the shape sweep with frequentist (mgcv)
# twins. The wash-out mechanism is specifically about PENALIZATION -- a penalized
# smooth losing a variance competition to a village random effect -- and mgcv's
# REML smoothing-parameter selection is not the same machinery as brms's prior on
# the smooth's wiggliness SD. So the key contrast is re-run here in brms, on the
# actual estimators the paper uses, to confirm the frequentist proxy is not
# hiding a prior-driven difference.
#
# Deliberately coarse: two truths (hump at matched signal, and a
# clearly-detectable hump), two arms, N_SIM small. This is a confirmation that
# the frequentist twins are faithful, not an independent calibration -- the
# frequentist sweep is the estimate of record.
#
#   ml_gam     brm(y ~ s(age) + sex + s(index, k=4) + (1 | village_id))
#              the one-stage model.
#   two_stage  stage 1  brm(y ~ village_id + s(age) + sex) -> G-computed
#                       age/sex-standardized village means + posterior SE
#              stage 2  brm(adj_mean | se(adj_se, sigma = TRUE) ~ s(index, k=4))
#              i.e. the reported estimator, exactly.
#
# Reported per arm: mean posterior swing against the true swing (attenuation),
# P(interior peak), shape RMSE, and the posterior mean of the smooth's wiggliness
# SD (sds_sindustrial_index_1) -- the Bayesian analogue of the frequentist edf,
# and the most direct read on whether the wiggly component is being shrunk away.
#
# Outcome is modelled on the project's /1000 scale so the standing sds prior is
# on the right scale; results are converted back to m/s for reporting.
# The Stan models are compiled once and reused via update(recompile = FALSE).
#
# Run: Rscript code/_experiments/calibration-nonlinear-shapes-bayes.R
# Smoke: SHAPEB_SMOKE=1 (N_SIM = 2, very short chains).
###############################################################################

suppressMessages({ library(here); library(lme4); library(brms) })
set.seed(2138)

out_dir <- here("outputs", "_experiments", "calibration")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

SMOKE <- nzchar(Sys.getenv("SHAPEB_SMOKE"))
N_SIM <- if (SMOKE) 2 else 25
CFG   <- if (SMOKE) list(w = 200, i = 600, c = 2) else list(w = 500, i = 1500, c = 4)
# Stage 2 uses the REPORTED estimator's production sampler settings rather than
# the lighter CFG above. With se() + sigma = TRUE on ~25 points the posterior has
# a weakly-identified tail that short chains wander into, producing smooth SDs in
# the 1e16 range and destroying the arm's mean. Diagnostics are recorded per
# replicate, non-converged fits are excluded, and MEDIANS are reported.
S2CFG <- if (SMOKE) list(w = 200, i = 600, c = 2) else list(w = 1000, i = 3000, c = 4)
SCALE <- 1000                                   # project outcome scaling


# ---- Real design (identical to calibration-nonlinear-shapes.R) ----

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
villages <- levels(s$village_id)

m0 <- lmer(tibia_sos ~ age_c + male + (1 | village_id), data = s)
m1 <- lmer(tibia_sos ~ age_c + male + industrial_index + (1 | village_id), data = s)
mu <- fixef(m0)[["(Intercept)"]]; g_age <- fixef(m0)[["age_c"]]; d_sex <- fixef(m0)[["male"]]
tau_eff <- attr(lme4::VarCorr(m1)$village_id, "stddev")[[1]]
sig     <- sigma(m0)
MOVE    <- abs(fixef(m1)[["industrial_index"]]) * span

f_hump <- function(x, S) { x0 <- mean(rng); -(S / (span / 2)^2) * (x - x0)^2 }
TRUTHS <- list(hump  = MOVE, hump2 = 2 * MOVE)

grid  <- seq(rng[1], rng[2], length.out = 41)
ctr   <- function(v) v - mean(v)
swing <- function(v) max(v) - min(v)

sim_frame <- function(S) {
  u <- rnorm(nv, 0, tau_eff)
  y <- mu + g_age * s$age_c + d_sex * s$male + u[vid] + f_hump(index_i, S) +
       rnorm(nrow(s), 0, sig)
  data.frame(y = y / SCALE, age_years = s$age_years, sex = s$sex,
             village_id = s$village_id, industrial_index = s$industrial_index)
}

converged <- function(fit) {
  rh <- suppressWarnings(max(brms::rhat(fit), na.rm = TRUE))
  dv <- sum(subset(brms::nuts_params(fit), Parameter == "divergent__")$Value)
  c(rhat = rh, divergent = dv, ok = as.numeric(is.finite(rh) && rh < 1.05))
}

# summaries of an AERF draws matrix (draws x grid), already in m/s
summarize_draws <- function(ep, f_true_grid) {
  fbar <- colMeans(ep)
  c(swing      = swing(fbar),                                  # posterior-mean curve
    swing_draw = mean(apply(ep, 1, swing)),                    # upward-biased, kept for reference
    int_peak   = mean(apply(ep, 1, which.max) %in% 2:(length(grid) - 1)),
    rmse       = sqrt(mean((ctr(fbar) - f_true_grid)^2)))
}


# ---- Compile both models once on a first dataset, then reuse ----

d0 <- sim_frame(MOVE)
prior_sds <- set_prior("student_t(3, 0, 2.5)", class = "sds")

fit_ml0 <- brm(bf(y ~ s(age_years, k = 5) + sex + s(industrial_index, k = 4) +
                    (1 | village_id), family = gaussian()),
               data = d0, prior = prior_sds, warmup = CFG$w, iter = CFG$i,
               chains = CFG$c, cores = CFG$c, seed = 2138, backend = "cmdstanr",
               control = list(adapt_delta = 0.95, max_treedepth = 12), refresh = 0, silent = 2)

fit_s10 <- brm(bf(y ~ village_id + s(age_years, k = 5) + sex, family = gaussian()),
               data = d0, prior = prior_sds, warmup = CFG$w, iter = CFG$i,
               chains = CFG$c, cores = CFG$c, seed = 2138, backend = "cmdstanr",
               control = list(adapt_delta = 0.95, max_treedepth = 12), refresh = 0, silent = 2)

vill_means <- function(fit, dd) {
  cohort <- dd[, c("age_years", "sex")]
  ad <- vapply(villages, function(v) {
    x <- cohort; x$village_id <- factor(v, levels = villages)
    rowMeans(posterior_epred(fit, newdata = x))
  }, numeric(brms::ndraws(fit))) * SCALE
  data.frame(adj_mean = colMeans(ad), adj_se = apply(ad, 2, sd),
             industrial_index = index_v)
}

# brms's default priors are computed once at compile time and carried through
# update(), so they would stay scaled to the FIRST replicate's village means. Set
# them explicitly at the scale brms's data-scaled default would choose anyway
# (village means are in m/s with SD ~70-100), so every replicate gets the same
# weakly-informative prior rather than a stale one.
prior_s2 <- set_prior("student_t(3, 0, 70)", class = "sigma") +
            set_prior("student_t(3, 0, 70)", class = "sds")

fit_s20 <- brm(bf(adj_mean | se(adj_se, sigma = TRUE) ~ s(industrial_index, k = 4)),
               data = vill_means(fit_s10, d0), family = gaussian(), prior = prior_s2,
               warmup = S2CFG$w, iter = S2CFG$i, chains = S2CFG$c, cores = S2CFG$c,
               seed = 2138, backend = "cmdstanr",
               control = list(adapt_delta = 0.99), refresh = 0, silent = 2)


# ---- Sweep ----

rows <- list()
for (tn in names(TRUTHS)) {
  S <- TRUTHS[[tn]]
  f_true_grid <- ctr(f_hump(grid, S))
  acc <- list(ml_gam = list(), two_stage = list())

  for (r in seq_len(N_SIM)) {
    dd <- sim_frame(S)
    cohort <- dd[, c("age_years", "sex")]

    # one-stage
    fml <- update(fit_ml0, newdata = dd, recompile = FALSE, refresh = 0, silent = 2)
    ep_ml <- vapply(grid, function(g) {
      x <- cohort; x$industrial_index <- g
      rowMeans(posterior_epred(fml, newdata = x, re_formula = NA)) * SCALE
    }, numeric(brms::ndraws(fml)))
    sds_ml <- mean(as.data.frame(fml)[["sds_sindustrial_index_1"]]) * SCALE

    # two-stage
    f1 <- update(fit_s10, newdata = dd, recompile = FALSE, refresh = 0, silent = 2)
    vm <- vill_means(f1, dd)
    f2 <- update(fit_s20, newdata = vm, recompile = FALSE, refresh = 0, silent = 2)
    ep_ts <- posterior_epred(f2, newdata = data.frame(industrial_index = grid,
                                                      adj_se = mean(vm$adj_se)))
    sds_ts <- mean(as.data.frame(f2)[["sds_sindustrial_index_1"]])

    acc$ml_gam[[r]]    <- c(summarize_draws(ep_ml, f_true_grid), sds = sds_ml, converged(fml))
    acc$two_stage[[r]] <- c(summarize_draws(ep_ts, f_true_grid), sds = sds_ts, converged(f2))
    message(sprintf("  %s rep %d/%d done", tn, r, N_SIM))
  }

  cat(sprintf("\n=== %s  (interior-peak hump, true swing %.0f m/s; %d villages, N_SIM = %d) ===\n",
              toupper(tn), swing(f_true_grid), nv, N_SIM))
  cat(sprintf("  %-11s %11s %11s %11s %11s %7s\n", "arm", "swing", "P(int pk)",
              "shape RMSE", "smooth SD", "n_ok"))
  for (a in names(acc)) {
    M  <- do.call(rbind, acc[[a]])
    ok <- M[, "ok"] == 1 & is.finite(M[, "swing"])
    m  <- apply(M[ok, , drop = FALSE], 2, median)      # median: robust to a stray bad fit
    cat(sprintf("  %-11s %11.1f %11.3f %11.1f %11.1f %7d\n",
                a, m[["swing"]], m[["int_peak"]], m[["rmse"]], m[["sds"]], sum(ok)))
    rows[[length(rows) + 1]] <- data.frame(truth = tn, arm = a,
                                           true_swing = swing(f_true_grid),
                                           swing = m[["swing"]], swing_draw = m[["swing_draw"]],
                                           int_peak = m[["int_peak"]],
                                           rmse = m[["rmse"]], smooth_sd = m[["sds"]],
                                           n_ok = sum(ok), n_sim = N_SIM)
    saveRDS(M, file.path(out_dir, sprintf("shape-bayes-perrep_%s_%s.rds", tn, a)))
  }
}

tab <- do.call(rbind, rows)
write.csv(tab, file.path(out_dir, "nonlinear-shape-calibration-bayes.csv"), row.names = FALSE)
cat("\nsaved", file.path(out_dir, "nonlinear-shape-calibration-bayes.csv"), "\n")
