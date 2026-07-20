###############################################################################
# A PRIORI power curves for the three community-level industrialization analyses,
# using the FREQUENTIST PROXY (validated against the Bayesian village two-stage estimator: matching
# FP 0.04/0.047 and coverage 0.96/0.954). A priori, NOT post-hoc: power is computed
# at pre-specified true effect sizes swept across a grid, independent of the
# observed estimates. N_SIM replicates per grid point.
#
# Estimator (two-stage cluster-honest, the frequentist twin of the village two-stage estimator):
#   stage 1  lm(Y ~ 0 + village + age_c + male)        -> adjusted village means + SE
#   stage 2  metafor::rma(means ~ index, REML, Hartung-Knapp)   [all three outcomes]
#   detect   95% slope CI excludes 0
#
# All three industrialization analyses (SOS, steps, ENMO) are summarized in the
# manuscript by the linear-projection slope of the AERF, so the power proxy detects
# every one through the second-stage meta-regression slope.
#
# DGP per outcome: mu + g_age*age_c + d_sex*male + u[village] + EFFECT + resid,
#   with (mu, g_age, d_sex, tau=between-village SD, sig=residual SD) estimated from
#   the real data (lme4), on the real design (real villages/sizes/age/sex/index),
#   and EFFECT = beta * index_i (grid in across-gradient units = beta*range).
#
# Outputs: outputs/_experiments/power-curves/industrialization-power.rds + .csv
#   (one row per analysis x grid point: power; plus landmark effect sizes).
#
# Run: script -q /dev/null Rscript code/_experiments/power-curves-industrialization.R < /dev/null
# Smoke: POW_SMOKE=1 (tiny N_SIM, few grid points).
###############################################################################

suppressMessages({ library(here); library(lme4); library(metafor) })
set.seed(2138)

out_dir <- here("outputs", "_experiments", "power-curves")
if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE)
  try(system2("xattr", c("-w", "com.dropbox.ignored", "1", out_dir)), silent = TRUE) }

SMOKE <- nzchar(Sys.getenv("POW_SMOKE"))
N_SIM <- if (SMOKE) 40 else 2000
NGRID <- if (SMOKE) 4  else 11

d <- read.csv(here("data", "processed", "Orang-Asli-pa-vs-bone-60126.csv"))

# analysis registry: outcome column, shape, label, units, grid max, resolution-floor landmark
analyses <- list(
  list(key = "sos-urb",   y = "tibia_sos",              lab = "Industrialization → tibial SOS",
       unit = "m/s (across gradient)",  grid_max = 200,  resolution = 32),   # qUS test-retest SEM
  list(key = "steps-urb", y = "ad_tot_step_count_0_24hr", lab = "Industrialization → daily steps",
       unit = "steps (across gradient)", grid_max = NA,   resolution = NA),   # grid_max = 1.5*SD (set below)
  list(key = "enmo-urb",  y = "ad_mean_enmo_mg_0_24hr",  lab = "Industrialization → mean ENMO",
       unit = "mg (across gradient)",    grid_max = NA,   resolution = NA)
)

# ---- per-analysis power curve ----

power_curve <- function(a) {
  s <- d[is.finite(d[[a$y]]) & is.finite(d$age_years) & !is.na(d$sex) &
         is.finite(d$industrial_index) & !is.na(d$village_id), ]
  s$village_id <- droplevels(factor(s$village_id))
  s$age_c <- as.numeric(scale(s$age_years))
  s$male  <- as.integer(s$sex == "male")
  s$Y     <- s[[a$y]]
  vid     <- as.integer(s$village_id)
  index_v <- as.numeric(tapply(s$industrial_index, s$village_id, function(x) x[1]))
  index_i <- index_v[vid]
  nv      <- nlevels(s$village_id)
  idx_rng <- diff(range(index_v))

  # DGP params from the real data
  m0    <- suppressMessages(lmer(Y ~ age_c + male + (1 | village_id), data = s))
  mu    <- fixef(m0)[["(Intercept)"]]; g_age <- fixef(m0)[["age_c"]]; d_sex <- fixef(m0)[["male"]]
  tau   <- attr(lme4::VarCorr(m0)$village_id, "stddev")[[1]]
  sig   <- sigma(m0)
  sd_y  <- sd(s$Y, na.rm = TRUE)

  # effect grid (in natural effect units) + landmarks
  grid_max <- if (is.na(a$grid_max)) 1.5 * sd_y else a$grid_max
  landmarks <- c(resolution = a$resolution, sd0.2 = 0.2 * sd_y, sd0.5 = 0.5 * sd_y)
  # include the exact landmark effect sizes as grid points (0.2-SD always; the
  # resolution floor where the outcome has one) so landmark power is a direct
  # simulated readout, not interpolated
  extra    <- c(landmarks[["sd0.2"]], if (!is.na(a$resolution)) a$resolution)
  eff_grid <- sort(unique(c(seq(0, grid_max, length.out = NGRID), extra)))

  detect_one <- function(eff) {
    beta <- eff / idx_rng                                      # across-gradient -> per-unit slope
    yvec <- mu + g_age * s$age_c + d_sex * s$male + rnorm(nv, 0, tau)[vid] + beta * index_i + rnorm(nrow(s), 0, sig)
    s$ysim <- yvec
    f1 <- lm(ysim ~ 0 + village_id + age_c + male, data = s)
    co <- summary(f1)$coefficients
    vr <- grep("^village_id", rownames(co))
    est <- co[vr, "Estimate"]; se <- co[vr, "Std. Error"]
    rr <- tryCatch(metafor::rma(yi = est, sei = se, mods = ~ index_v, method = "REML", test = "knha"),
                   error = function(e) NULL)
    if (is.null(rr)) return(NA)
    as.numeric(rr$ci.lb[2] > 0 | rr$ci.ub[2] < 0)              # slope CI excludes 0
  }

  pw <- vapply(eff_grid, function(eff) mean(replicate(N_SIM, detect_one(eff)), na.rm = TRUE), numeric(1))
  cat(sprintf("\n%s  (n=%d, %d villages, SD=%.3g, resid=%.3g, tau=%.3g)\n", a$lab, nrow(s), nv, sd_y, sig, tau))
  for (i in seq_along(eff_grid)) cat(sprintf("  effect=%8.2f  power=%.3f\n", eff_grid[i], pw[i]))
  data.frame(key = a$key, label = a$lab, unit = a$unit, effect = eff_grid, power = pw,
             resolution = landmarks[["resolution"]], sd0.2 = landmarks[["sd0.2"]], sd0.5 = landmarks[["sd0.5"]])
}

if (SMOKE) analyses <- analyses[1:2]
res <- do.call(rbind, lapply(analyses, power_curve))
saveRDS(res, file.path(out_dir, "industrialization-power.rds"))
write.csv(res, file.path(out_dir, "industrialization-power.csv"), row.names = FALSE)
cat("\nsaved", file.path(out_dir, "industrialization-power.csv"), "\n")
