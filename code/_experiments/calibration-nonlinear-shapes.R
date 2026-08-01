###############################################################################
# Does the village random intercept wash out a NON-LINEAR index effect?
#
# calibration-estimator-comparison.R found the one-stage multilevel GAM and the
# two-stage estimator indistinguishable, with no attenuation. But its DGP truth
# was LINEAR, and a linear trend lives almost entirely in the smooth's
# UNPENALIZED null space (brms/mgcv split s(x) into an unpenalized linear term
# plus a penalized wiggly part). The wash-out mechanism -- a penalized smooth
# losing a variance-competition against an unpenalized-cost random intercept --
# can only act on the wiggly part. So that run could not have detected a
# wash-out even if one existed. This script supplies the missing case.
#
# It also supplies a calibration check across simulated null, monotonic and
# hump truths.
#
# DESIGN. Real SOS analytic sample: 25 villages, actual sizes 6-79, real age /
# sex / per-village index, variance components from the fitted data. Four truths,
# all at the village level, with the MONOTONIC and HUMP truths scaled to the same
# curve movement (92 m/s) so shape is the only thing that differs:
#   null   f(x) = 0
#   mono   linear, across-gradient change 92 m/s  (the previous run's truth)
#   hump   interior-peak quadratic, peak-to-trough swing 92 m/s
#   hump2  same shape, swing 184 m/s (a clearly-detectable regime, to separate
#          "washed out" from "too small to see at 25 villages")
# The hump is the substantively motivated non-linear shape: the activity-vs-
# industrialization inverted-U seen in the observed data.
#
# ARMS. Only smooth-capable estimators; the linear arms cannot represent a hump.
#   single_gam  gam(y ~ s(age) + male + s(index, k=4))                  no clustering
#   ml_gam      + s(village_id, bs="re")                                the one-stage model
#   two_stage   lm village-FE means -> gam(vmean ~ s(index,k=4) + s(village,bs="re"),
#               weights = 1/se^2, scale = 1)
#               NOTE this is an exact frequentist analogue of the reported
#               `adj_mean | se(adj_se, sigma=TRUE) ~ s(index, k=4)`: fixing the
#               scale with 1/se^2 weights makes the sampling variance exactly
#               se_j^2 and the village random effect adds tau^2 on top. Writing
#               it out makes the structural point visible -- stage two ALSO
#               contains a penalized index smooth sitting beside a village-level
#               iid variance component, so if the wash-out mechanism were real it
#               would afflict the two-stage estimator too.
#
# METRICS. `edf` is the direct wash-out diagnostic: if the index smooth is being
# shrunk to flat, its effective degrees of freedom collapse toward 1.
#   p_sig     P(approximate p-value for s(index) < 0.05)  -> FP under null, power otherwise
#   edf       mean effective df of the index smooth
#   swing     mean estimated peak-to-trough of the fitted curve vs the truth (attenuation)
#   rmse      mean RMSE of the centred fitted curve against the centred truth (shape recovery)
#   int_peak  P(fitted curve's maximum is strictly interior)
#
# Frequentist twins throughout: the validated frequentist proxy is what these
# calibrations use. mgcv fits both the index smooth
# and the village random effect as penalized terms selected by REML, so the
# variance-competition the wash-out argument is about is represented faithfully.
# calibration-nonlinear-shapes-bayes.R re-checks the key contrast in brms.
#
# Run: Rscript code/_experiments/calibration-nonlinear-shapes.R
# Smoke: SHAPE_SMOKE=1 (N_SIM = 15).
###############################################################################

suppressMessages({ library(here); library(lme4); library(mgcv) })
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
nv      <- nlevels(s$village_id)
rng     <- range(index_v)
span    <- diff(rng)

m0 <- lmer(tibia_sos ~ age_c + male + (1 | village_id), data = s)
m1 <- lmer(tibia_sos ~ age_c + male + industrial_index + (1 | village_id), data = s)
mu <- fixef(m0)[["(Intercept)"]]; g_age <- fixef(m0)[["age_c"]]; d_sex <- fixef(m0)[["male"]]
tau_null <- attr(lme4::VarCorr(m0)$village_id, "stddev")[[1]]
tau_eff  <- attr(lme4::VarCorr(m1)$village_id, "stddev")[[1]]
sig      <- sigma(m0)
MOVE     <- abs(fixef(m1)[["industrial_index"]]) * span      # ~92 m/s of curve movement


# ---- Truths: village-level mean functions, centred so only shape matters ----

f_null <- function(x) rep(0, length(x))
f_mono <- function(x) -(MOVE / span) * (x - mean(rng))
f_hump <- function(x, S = MOVE) { x0 <- mean(rng); -(S / (span / 2)^2) * (x - x0)^2 }

TRUTHS <- list(
  null  = list(f = f_null,                              tau = tau_null, lab = "NULL  (flat)"),
  mono  = list(f = f_mono,                              tau = tau_eff,
               lab = sprintf("MONO  (linear, %.0f m/s across gradient)", MOVE)),
  hump  = list(f = function(x) f_hump(x, MOVE),         tau = tau_eff,
               lab = sprintf("HUMP  (interior peak, swing %.0f m/s)", MOVE)),
  hump2 = list(f = function(x) f_hump(x, 2 * MOVE),     tau = tau_eff,
               lab = sprintf("HUMP2 (interior peak, swing %.0f m/s)", 2 * MOVE)))

N_SIM <- if (nzchar(Sys.getenv("SHAPE_SMOKE"))) 15 else 500
ARMS  <- c("single_gam", "ml_gam", "two_stage")
grid  <- seq(rng[1], rng[2], length.out = 41)

ctr   <- function(v) v - mean(v)
swing <- function(v) max(v) - min(v)


# ---- One replicate ----

fit_arms <- function(y) {
  s$y <- y
  blank <- c(p = NA_real_, edf = NA_real_, swing = NA_real_, rmse = NA_real_,
             ipk = NA_real_, setNames(rep(NA_real_, length(grid)),
                                      paste0("fx", seq_along(grid))))
  out <- lapply(ARMS, function(a) blank)
  names(out) <- ARMS
  nd <- data.frame(industrial_index = grid, age_c = 0, male = 0,
                   village_id = factor(levels(s$village_id)[1], levels = levels(s$village_id)))

  summarize_gam <- function(gm, newdata, term = "s(industrial_index)") {
    st <- summary(gm)$s.table
    i  <- match(term, rownames(st))
    pt <- predict(gm, newdata = newdata, type = "terms", terms = term)
    fx <- as.numeric(pt[, term])
    c(p = st[i, "p-value"], edf = st[i, "edf"], swing = swing(fx),
      rmse = NA_real_, ipk = as.numeric(which.max(fx) %in% 2:(length(grid) - 1)),
      setNames(fx, paste0("fx", seq_along(fx))))
  }

  g1 <- tryCatch(gam(y ~ s(age_c, k = 5) + male + s(industrial_index, k = 4),
                     data = s, method = "REML"), error = function(e) NULL)
  g2 <- tryCatch(gam(y ~ s(age_c, k = 5) + male + s(industrial_index, k = 4) +
                       s(village_id, bs = "re"), data = s, method = "REML"),
                 error = function(e) NULL)

  # two-stage: village-FE adjusted means, then the se(sigma=TRUE) analogue
  f1 <- lm(y ~ 0 + village_id + age_c + male, data = s)
  c1 <- summary(f1)$coefficients
  vr <- grep("^village_id", rownames(c1))
  vd <- data.frame(vmean = c1[vr, "Estimate"], vse = c1[vr, "Std. Error"],
                   industrial_index = index_v,
                   village_id = factor(levels(s$village_id), levels = levels(s$village_id)))
  g3 <- tryCatch(gam(vmean ~ s(industrial_index, k = 4) + s(village_id, bs = "re"),
                     data = vd, weights = 1 / vd$vse^2, scale = 1, method = "REML"),
                 error = function(e) NULL)

  nd2 <- data.frame(industrial_index = grid,
                    village_id = factor(levels(s$village_id)[1], levels = levels(s$village_id)))
  for (nm in ARMS) {
    gm <- switch(nm, single_gam = g1, ml_gam = g2, two_stage = g3)
    if (is.null(gm)) next
    out[[nm]] <- summarize_gam(gm, if (nm == "two_stage") nd2 else nd)
  }
  unlist(out)
}


# ---- Sweep ----

rows <- list()
for (tn in names(TRUTHS)) {
  tr <- TRUTHS[[tn]]
  f_true_grid <- ctr(tr$f(grid))
  s_true <- swing(f_true_grid)

  res <- t(replicate(N_SIM, {
    u <- rnorm(nv, 0, tr$tau)
    y <- mu + g_age * s$age_c + d_sex * s$male + u[vid] + tr$f(index_i) + rnorm(nrow(s), 0, sig)
    fit_arms(y)
  }))

  cat(sprintf("\n=== %s ===\n  %d villages, n = %d, village SD %.0f, resid SD %.0f, N_SIM = %d, true swing %.0f m/s\n",
              tr$lab, nv, nrow(s), tr$tau, sig, N_SIM, s_true))
  cat(sprintf("  %-11s %8s %7s %11s %10s %10s\n", "arm",
              if (tn == "null") "FP" else "power", "edf", "mean swing", "shape RMSE", "P(int pk)"))
  for (a in ARMS) {
    fxcols <- grep(sprintf("^%s\\.fx", a), colnames(res))
    fx     <- res[, fxcols, drop = FALSE]
    rmse   <- mean(apply(fx, 1, function(r) sqrt(mean((ctr(r) - f_true_grid)^2))), na.rm = TRUE)
    r <- data.frame(truth = tn, arm = a, true_swing = s_true,
                    p_sig    = mean(res[, paste0(a, ".p")] < 0.05, na.rm = TRUE),
                    edf      = mean(res[, paste0(a, ".edf")], na.rm = TRUE),
                    swing    = mean(res[, paste0(a, ".swing")], na.rm = TRUE),
                    rmse     = rmse,
                    int_peak = mean(res[, paste0(a, ".ipk")], na.rm = TRUE))
    cat(sprintf("  %-11s %8.3f %7.2f %11.1f %10.1f %10.3f\n",
                a, r$p_sig, r$edf, r$swing, r$rmse, r$int_peak))
    rows[[length(rows) + 1]] <- r
  }
}

tab <- do.call(rbind, rows)
write.csv(tab, file.path(out_dir, "nonlinear-shape-calibration.csv"), row.names = FALSE)
cat("\nsaved", file.path(out_dir, "nonlinear-shape-calibration.csv"), "\n")
