###############################################################################
# EXPERIMENT (Set 2, Option H): standardized village-level meta-GAM.
#   PRIMARY cluster-level anchor (per the 2026-06-15 second opinion).
#
# PREPPED, NOT YET RUN. Status: untested until executed.
#
# Two stages, inference at the level the exposure varies (the village):
#   Stage 1 (within): outcome ~ village_id + s(age_years, k=5) + sex
#     -> age/sex-standardized village means + posterior SE (one per village),
#        by predicting each village over the cohort's age/sex distribution.
#   Stage 2 (between): adj_mean | se(adj_se, sigma=TRUE) ~ s(industrial_index, k=4)
#     -> a Bayesian measurement-error meta-GAM on the ~24-31 village means. The
#        low-dimensional smooth can carry the SOS gradient AND the activity hump;
#        `se()` puts the stage-1 uncertainty in the likelihood (precision-weighting
#        emerges); sigma=TRUE estimates residual between-village SD.
# For SOS, also fit a linear stage-2 variant for comparison.
#
# Uncertainty is honest by construction (N = #villages). Note: the `se()` plug-in
# treats adj_se as known; full posterior propagation (carry stage-1 draws into
# stage 2) is the upgrade -- TODO if a reviewer presses.
#
# ISOLATED: writes only to outputs/_experiments/set2-option-h/. pty + sandbox off.
# Smoke: SET2_SMOKE=1 (1 spec, fast stage-1 MCMC).
###############################################################################

library(here)
source(here("code", "_startup", "init.R"))
suppressMessages({ library(brms); library(dplyr) })

make_out <- function(name) {
  p <- here("outputs", "_experiments", name)
  if (!dir.exists(p)) { dir.create(p, recursive = TRUE)
    try(system2("xattr", c("-w", "com.dropbox.ignored", "1", p)), silent = TRUE) }
  p
}
hpdi <- function(x, m = 0.95) { x <- sort(x); n <- length(x); k <- floor(m * n)
  i <- which.min(x[(k + 1):n] - x[1:(n - k)]); c(x[i], x[i + k]) }

SMOKE <- nzchar(Sys.getenv("SET2_SMOKE"))
specs <- c("sos-urb", "steps-urb", "enmo-urb")
if (SMOKE) specs <- "sos-urb"
cfg <- if (SMOKE) list(w = 500, i = 1500, c = 4) else list(w = WARMUP, i = ITER, c = CHAINS)
out_dir <- make_out("set2-option-h")

for (key in specs) {
  spec <- model_templates[[key]]; ex <- "industrial_index"; sf <- spec$outcome_scale_factor
  resp <- all.vars(spec$bf$formula)[1]
  d <- prep_local_data(spec, dat); d <- d[!is.na(d$village_id) & !is.na(d$industrial_index), ]
  d$village_id <- droplevels(factor(d$village_id))
  villages <- levels(d$village_id)

  # ---- Stage 1: village FE + age/sex; standardize each village over cohort age/sex ----
  s1_form <- bf(as.formula(sprintf("%s ~ village_id + s(age_years, k = 5) + sex", resp)),
                family = gaussian(link = "identity"))
  message(sprintf("Option H stage 1: %s  (n=%d, %d villages)", key, nrow(d), length(villages)))
  m1 <- brm(s1_form, data = d,
            prior = set_prior("student_t(3, 0, 2.5)", class = "sds"),
            warmup = cfg$w, iter = cfg$i, thin = THIN, chains = cfg$c, cores = cfg$c,
            seed = SEED, backend = "cmdstanr", control = BRMS_CONTROL, refresh = 0, silent = 2)

  cohort <- d[, c("age_years", "sex")]
  ndraws <- brms::ndraws(m1)
  adj_draws <- matrix(NA_real_, nrow = ndraws, ncol = length(villages))   # standardized means, scaled-outcome units
  for (k in seq_along(villages)) {
    x <- cohort; x$village_id <- factor(villages[k], levels = villages)
    adj_draws[, k] <- rowMeans(posterior_epred(m1, newdata = x))
  }
  adj_draws <- adj_draws * sf                                              # -> natural units
  villdf <- data.frame(
    village_id       = villages,
    industrial_index = vapply(villages, function(v) d$industrial_index[d$village_id == v][1], numeric(1)),
    adj_mean         = colMeans(adj_draws),
    adj_se           = apply(adj_draws, 2, sd))
  saveRDS(m1, file.path(out_dir, paste0("stage1-fit_", key, ".rds")))
  saveRDS(villdf, file.path(out_dir, paste0("village-means_", key, ".rds")))

  # ---- Stage 2: Bayesian measurement-error meta-GAM on the village means ----
  fit_smooth <- brm(bf(adj_mean | se(adj_se, sigma = TRUE) ~ s(industrial_index, k = 4)),
                    data = villdf, family = gaussian(),
                    warmup = 1000, iter = 3000, chains = 4, cores = 4,
                    seed = SEED, backend = "cmdstanr", control = list(adapt_delta = 0.99),
                    refresh = 0, silent = 2)

  grid <- seq(quantile(villdf$industrial_index, .01), quantile(villdf$industrial_index, .99), length.out = 41)
  nd   <- data.frame(industrial_index = grid, adj_se = mean(villdf$adj_se))   # adj_se dummy (epred ignores it)
  ep   <- brms::posterior_epred(fit_smooth, newdata = nd)                     # draws x grid (AERF)
  ng   <- length(grid)

  endpt  <- ep[, ng] - ep[, 1]                              # across-gradient contrast
  swing  <- apply(ep, 1, function(r) max(r) - min(r))       # peak-to-trough (hump magnitude)
  pk_loc <- grid[apply(ep, 1, which.max)]                   # argmax location
  interior <- mean(apply(ep, 1, which.max) %in% 2:(ng - 1)) # P(interior peak)
  he <- hpdi(endpt); hs <- hpdi(swing); hp <- hpdi(pk_loc)

  cat(sprintf("\n=== %s  (Option H, %s; %d villages) ===\n", key, spec$outcome_label, length(villages)))
  cat(sprintf("  smooth stage-2: across-gradient Δ HPDI [%.4g, %.4g]  P(decline)=%.3f | swing HPDI [%.4g, %.4g]\n",
      he[1], he[2], mean(endpt < 0), hs[1], hs[2]))
  cat(sprintf("                  peak-location HPDI [%.4g, %.4g] | P(interior peak)=%.2f\n", hp[1], hp[2], interior))

  res <- data.frame(spec = key, n_villages = length(villages),
                    endpt_lo = he[1], endpt_hi = he[2], p_decline = mean(endpt < 0),
                    swing_lo = hs[1], swing_hi = hs[2],
                    peak_lo = hp[1], peak_hi = hp[2], p_interior_peak = interior)

  if (identical(key, "sos-urb")) {                          # linear stage-2 variant for the monotonic SOS headline
    fit_lin <- brm(bf(adj_mean | se(adj_se, sigma = TRUE) ~ industrial_index),
                   data = villdf, family = gaussian(),
                   warmup = 1000, iter = 3000, chains = 4, cores = 4,
                   seed = SEED, backend = "cmdstanr", refresh = 0, silent = 2)
    b <- as.data.frame(fit_lin)[["b_industrial_index"]]
    span <- diff(range(villdf$industrial_index))
    hb <- hpdi(b); hc <- hpdi(b * span)
    cat(sprintf("  linear stage-2: slope HPDI [%.4g, %.4g] P(<0)=%.3f | across-gradient Δ HPDI [%.4g, %.4g]\n",
        hb[1], hb[2], mean(b < 0), hc[1], hc[2]))
    res$lin_slope_lo <- hb[1]; res$lin_slope_hi <- hb[2]; res$lin_p_neg <- mean(b < 0)
  }
  saveRDS(list(villdf = villdf, stage2_fit = fit_smooth, aerf_draws = ep, grid = grid, summary = res),
          file.path(out_dir, paste0("stage2_", key, ".rds")))
  write.csv(res, file.path(out_dir, paste0("summary_", key, ".csv")), row.names = FALSE)
}
cat("\nOption H done — per-spec summaries in", out_dir, "\n")
