###############################################################################
# A PRIORI power curves for the three community-level industrialization
# analyses, under the one-stage estimator the results are reported from.
#
# A priori, not post-hoc: power at pre-specified true effect sizes swept across
# a grid, independent of the observed estimates. N_SIM replicates per grid
# point.
#
# ESTIMATOR PROXY. The reported model is a Bayesian additive GAM with a
# community random intercept. As in the physical-activity -> bone power
# analysis, power is computed with a fast frequentist twin rather than the
# Bayesian model itself:
#
#   fit     gam(y ~ s(age_c, k=5) + male + s(industrial_index, k=4) +
#                    s(village_id, bs="re"), method = "REML")
#   detect  across-gradient LINEAR PROJECTION's 95% t interval on G-2 df
#           excludes 0
#
# This is the structural analogue of the twin already used for the
# physical-activity analyses -- flexible age, exposure entering additively, and
# community identity carried in the way the design requires (a random intercept
# here, because the industrialization index is a community-level attribute;
# fixed effects there, because activity varies within community).
#
# ESTIMAND. The best LINEAR PROJECTION of the exposure-response function over a
# 21-point index grid, averaged over the cohort's age/sex distribution via the
# lpmatrix, then multiplied by the gradient span so it reads in across-gradient
# units. This is the quantity the manuscript reports, and it is what the
# calibration harness measures.
#
# It is NOT an endpoint contrast. The two are different functionals in practice:
# the fitted smooth wiggles even when the truth is linear, and the endpoint reads
# the two index values where communities are sparsest, so it is the noisier of
# the two. A power analysis has to characterise the statistic actually reported.
#
# Both are linear functionals of the coefficient vector, so each is a single
# contrast vector `cv` against coef(gm), with se from cv' V cv. `cv` depends only
# on the covariate design and the spline bases -- neither of which varies across
# replicates, since only the simulated response changes -- so it is built once
# per analysis in setup() and reused. Verified against a per-replicate rebuild
# at script start.
#
# The proxy is validated against the reported estimator on every run: its null
# rejection rate and its power at the estimated effect are compared with the
# production one-stage additive model's directly measured values, on the same
# SOS design, at the same effect magnitude, and on this same projection estimand
# (BAYES_REF below, measured with calibration-onestage-production.R, whose
# defaults are the settings used here). The script warns if either drifts.
#
# DGP: additive community-level effect, matching power-curves-pa-bone.R. The
# estimand is the linear-projection slope, a population-averaged scalar, so an
# age x exposure interaction would redistribute that effect across ages without
# changing its cohort-average value, and leaving one out of the DGP does not
# bias the sweep.
#
# Outputs: outputs/_experiments/power-curves/industrialization-power-onestage.{rds,csv}
#          outputs/_experiments/power-curves/onestage-proxy-validation.csv
#
# Run:   Rscript code/_experiments/power-curves-industrialization-onestage.R
# Smoke: POW1S_SMOKE=1
###############################################################################

suppressMessages({ library(here); library(lme4); library(mgcv) })
set.seed(2138)

out_dir <- here("outputs", "_experiments", "power-curves")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

SMOKE <- nzchar(Sys.getenv("POW1S_SMOKE"))
N_SIM <- if (SMOKE) 20 else as.integer(Sys.getenv("POW1S_N", "2000"))
NGRID <- if (SMOKE) 3  else 11
N_VAL <- if (SMOKE) 20 else as.integer(Sys.getenv("POW1S_NVAL", "500"))

# Partial-run knobs, both empty by default so a bare run is the full sweep.
# POW1S_KEYS    comma-separated analysis keys to run, e.g. "steps-urb"
# POW1S_EFFECTS comma-separated effect sizes to evaluate INSTEAD of the grid
# A partial run writes to a suffixed filename, so re-measuring one cell cannot
# overwrite a full sweep's outputs. Adding a knob here rather than a second
# script keeps one implementation of the estimator: a separate one would be free
# to drift from the sweep it is meant to extend.
ONLY_KEYS <- Sys.getenv("POW1S_KEYS")
ONLY_EFFS <- Sys.getenv("POW1S_EFFECTS")
PARTIAL   <- nzchar(ONLY_KEYS) || nzchar(ONLY_EFFS)

# Reference operating characteristics for the REPORTED estimator -- the one-stage
# additive model at the production priors -- on the SOS design, at the same
# effect magnitude and on this same linear-projection estimand, measured over
# N = 300 replicates. The frequentist twin below is checked against them on every
# run: the twin is a stand-in for a Bayesian model, and these say how good a
# stand-in it is. Regenerate them with calibration-onestage-production.R, whose
# defaults are this cell.
BAYES_REF <- c(null = 0.063, effect = 0.380)

d <- read.csv(here("data", "processed", "orang-asli-pa-bone-data.csv"))

analyses <- list(
  list(key = "sos-urb",   y = "tibia_sos",                lab = "Industrialization → tibial SOS",
       unit = "m/s (across gradient)",   grid_max = 200, resolution = 32),
  # Daily steps has an instrument floor: the accelerometer's step-count error at
  # the sample median. Mean daily ENMO has none, hence NA there.
  list(key = "steps-urb", y = "ad_tot_step_count_0_24hr", lab = "Industrialization → daily steps",
       unit = "steps (across gradient)", grid_max = NA,  resolution = 1200),
  list(key = "enmo-urb",  y = "ad_mean_enmo_mg_0_24hr",   lab = "Industrialization → mean ENMO",
       unit = "mg (across gradient)",    grid_max = NA,  resolution = NA)
)


# ---- Design + DGP -------------------------------------------------------------

setup <- function(a) {
  s <- d[is.finite(d[[a$y]]) & is.finite(d$age_years) & !is.na(d$sex) &
         is.finite(d$industrial_index) & !is.na(d$village_id), ]
  s$village_id <- droplevels(factor(s$village_id))
  s$age_c <- as.numeric(scale(s$age_years)); s$male <- as.integer(s$sex == "male")
  s$Y <- s[[a$y]]
  vid <- as.integer(s$village_id)
  iv  <- as.numeric(tapply(s$industrial_index, s$village_id, function(x) x[1]))
  m0  <- suppressMessages(lmer(Y ~ age_c + male + (1 | village_id), data = s))
  idx_rng <- diff(range(iv))

  # The projection contrast vector. Built from a fit to the OBSERVED response
  # purely to obtain the spline bases: mgcv places knots from the covariates, so
  # the bases -- and hence `cv` -- are identical for every simulated response on
  # this design. Age, sex and the community random-effect columns are constant
  # across grid points and so drop out of the projection, which is what makes
  # this marginal of community, matching re_formula = NA in the reported model.
  gm0  <- gam(Y ~ s(age_c, k = 5) + male + s(industrial_index, k = 4) +
                  s(village_id, bs = "re"), data = s, method = "REML")
  grid <- seq(min(iv), max(iv), length.out = 21)
  gxc  <- grid - mean(grid)
  M    <- vapply(grid, function(v)
            colMeans(predict(gm0, newdata = transform(s, industrial_index = v),
                             type = "lpmatrix")),
            numeric(length(coef(gm0))))
  cv <- as.numeric(M %*% gxc / sum(gxc^2)) * idx_rng   # across-gradient units

  list(s = s, vid = vid, index_i = iv[vid], nv = nlevels(s$village_id),
       idx_rng = idx_rng, cv = cv,
       mu = fixef(m0)[[1]], g_age = fixef(m0)[[2]], d_sex = fixef(m0)[[3]],
       tau = attr(lme4::VarCorr(m0)$village_id, "stddev")[[1]],
       sig = sigma(m0), sd_y = sd(s$Y, na.rm = TRUE))
}

detect_one <- function(P, eff) {
  b <- eff / P$idx_rng
  s <- P$s
  s$ysim <- P$mu + P$g_age * s$age_c + P$d_sex * s$male +
            rnorm(P$nv, 0, P$tau)[P$vid] + b * P$index_i + rnorm(nrow(s), 0, P$sig)
  gm <- tryCatch(gam(ysim ~ s(age_c, k = 5) + male + s(industrial_index, k = 4) +
                       s(village_id, bs = "re"), data = s, method = "REML"),
                 error = function(e) NULL)
  if (is.null(gm)) return(NA_real_)
  est <- as.numeric(P$cv %*% coef(gm))
  se  <- sqrt(as.numeric(t(P$cv) %*% vcov(gm, unconditional = TRUE) %*% P$cv))
  as.numeric(abs(est / se) > qt(0.975, df = P$nv - 2))
}


# `cv` is precomputed from the observed-response fit. Confirm on one simulated
# replicate that a rebuilt basis gives the same contrast, so the reuse cannot
# silently drift if a future mgcv changes how knots are placed.
check_cv <- function(P) {
  s <- P$s
  s$ysim <- P$mu + P$g_age * s$age_c + P$d_sex * s$male +
            rnorm(P$nv, 0, P$tau)[P$vid] + rnorm(nrow(s), 0, P$sig)
  gm   <- gam(ysim ~ s(age_c, k = 5) + male + s(industrial_index, k = 4) +
                s(village_id, bs = "re"), data = s, method = "REML")
  grid <- seq(min(s$industrial_index), max(s$industrial_index), length.out = 21)
  gxc  <- grid - mean(grid)
  M    <- vapply(grid, function(v)
            colMeans(predict(gm, newdata = transform(s, industrial_index = v),
                             type = "lpmatrix")),
            numeric(length(coef(gm))))
  cv2 <- as.numeric(M %*% gxc / sum(gxc^2)) * P$idx_rng
  max(abs(cv2 - P$cv))
}


# ---- Validate the proxy against the reported estimator ------------------------
# Skipped on a partial run: it costs 2 x N_VAL fits on a design the partial run
# may not touch, and would rewrite an output that run is not recomputing.

if (!PARTIAL) {

cat("\n#### Proxy validation against the reported Bayesian estimator (SOS design) ####\n")
P_sos <- setup(analyses[[1]])
cat(sprintf("  contrast-vector reuse check: max |cv_rebuilt - cv| = %.3g\n", check_cv(P_sos)))
m1 <- suppressMessages(lmer(Y ~ age_c + male + industrial_index + (1 | village_id), data = P_sos$s))
eff_at_estimate <- fixef(m1)[["industrial_index"]] * P_sos$idx_rng

val <- do.call(rbind, lapply(list(list(lab = "null",   eff = 0),
                                  list(lab = "effect", eff = eff_at_estimate)), function(cd) {
  data.frame(condition = cd$lab, effect = cd$eff, n_sim = N_VAL,
             proxy = mean(replicate(N_VAL, detect_one(P_sos, cd$eff)), na.rm = TRUE),
             bayesian_reported = BAYES_REF[[cd$lab]])
}))
print(val, row.names = FALSE)
write.csv(val, file.path(out_dir, "onestage-proxy-validation.csv"), row.names = FALSE)

if (!SMOKE) {
  if (abs(val$proxy[1] - BAYES_REF[["null"]]) > 0.04)
    warning("proxy null rate has drifted from the reported estimator; re-validate")
  if (abs(val$proxy[2] - BAYES_REF[["effect"]]) > 0.10)
    warning("proxy power has drifted from the reported estimator; re-validate")
}

}   # end !PARTIAL


# ---- Sweep --------------------------------------------------------------------

power_curve <- function(a) {
  P <- setup(a)
  grid_max  <- if (is.na(a$grid_max)) 1.5 * P$sd_y else a$grid_max
  landmarks <- c(resolution = a$resolution, sd0.2 = 0.2 * P$sd_y, sd0.5 = 0.5 * P$sd_y)
  extra     <- c(landmarks[["sd0.2"]], if (!is.na(a$resolution)) a$resolution)
  eff_grid  <- sort(unique(c(seq(0, grid_max, length.out = NGRID), extra)))
  if (nzchar(ONLY_EFFS))
    eff_grid <- sort(as.numeric(strsplit(ONLY_EFFS, ",")[[1]]))

  cat(sprintf("\n%s  (n=%d, %d communities, SD=%.3g, resid=%.3g, tau=%.3g)\n",
              a$lab, nrow(P$s), P$nv, P$sd_y, P$sig, P$tau))
  pw <- vapply(eff_grid, function(eff) {
    p <- mean(replicate(N_SIM, detect_one(P, eff)), na.rm = TRUE)
    cat(sprintf("  effect=%9.2f  power=%.3f\n", eff, p)); p
  }, numeric(1))

  data.frame(key = a$key, label = a$lab, unit = a$unit, effect = eff_grid, power = pw,
             resolution = landmarks[["resolution"]], sd0.2 = landmarks[["sd0.2"]],
             sd0.5 = landmarks[["sd0.5"]], n = nrow(P$s), n_villages = P$nv)
}

if (SMOKE) analyses <- analyses[1:2]
if (nzchar(ONLY_KEYS)) {
  keep <- strsplit(ONLY_KEYS, ",")[[1]]
  analyses <- Filter(function(a) a$key %in% keep, analyses)
  if (!length(analyses)) stop("POW1S_KEYS matched no analysis key")
}
res <- do.call(rbind, lapply(analyses, power_curve))

stem <- if (PARTIAL) "industrialization-power-onestage-partial" else
                     "industrialization-power-onestage"
saveRDS(res, file.path(out_dir, paste0(stem, ".rds")))
write.csv(res, file.path(out_dir, paste0(stem, ".csv")), row.names = FALSE)
cat("\nWrote ", stem, ".{rds,csv} -> ", out_dir, "\n", sep = "")
