###############################################################################
# Industrialization analysis: two-stage estimator with FULL propagation of the
# stage-one posterior, compared against the se() plug-in currently reported.
#
# WHY. The reported estimator (industrialization-village-two-stage.R) collapses the
# stage-one posterior of each village's standardized mean to (mean, SD) and feeds
# those into a Bayesian measurement-error meta-regression via
# `adj_mean | se(adj_se, sigma = TRUE)`. That plug-in makes two approximations:
#   (i)  each village mean's posterior is treated as Gaussian with KNOWN SD;
#   (ii) the village means are treated as INDEPENDENT measurements.
# (ii) is the substantive one: every village mean is standardized using the SAME
# posterior draw of the age smooth and the sex coefficient, so their posteriors
# are correlated by construction. `se()` cannot represent that.
#
# WHAT THIS SCRIPT DOES. Three estimators of the same across-gradient effect:
#
#   plugin   the reported estimator, read back from stage2_<key>.rds.
#
#   cut      "propagate the stage-one draws all the way through", implemented as
#            multiple imputation: draw M village-mean vectors from the stage-one
#            posterior, fit the stage-two smooth to each treating that vector as
#            exact, pool the draws (brm_multiple). This preserves the stage-one
#            correlation structure and any non-normality, because whole vectors
#            are carried forward rather than marginal summaries.
#            CAVEAT, stated up front: this is Plummer's *cut* posterior, not the
#            exact joint posterior. Stage two is never allowed to inform stage
#            one. Its known cost is that each imputation's residual SD absorbs
#            both genuine between-village variation AND stage-one measurement
#            noise, so `cut` is expected to be CONSERVATIVE (wider than exact).
#            It is reported because it is the literal reading of "propagate the
#            draws through", and because a conservative check that leaves the
#            conclusions intact is still informative.
#
#   joint    the exact fully-propagated model, for reference: one hierarchical
#            GAM on the individual-level data,
#              outcome ~ s(age) + sex + s(industrial_index, k = 4) + (1 | village_id)
#            Here the village effects have the index smooth as their prior mean
#            and tau as their SD, which IS the two-stage structure fit in one
#            step with no cut and no plug-in. Nothing is approximated.
#            NOTE: a penalized index smooth competing with a village random
#            intercept for the same between-village variance is the
#            configuration the Hodges & Reich spatial-confounding argument
#            concerns. Whether it bites at this design is settled empirically by
#            calibration-estimator-comparison.R; read the two together before
#            drawing a conclusion from `joint`.
#
# Also reports the mean off-diagonal correlation of the stage-one village-mean
# posterior, which quantifies exactly what approximation (ii) throws away.
#
# ISOLATED: writes only to outputs/_experiments/industrialization-village-two-stage-propagated/.
# Reuses the saved stage-one fits from industrialization-village-two-stage.R when present.
# Run: Rscript code/_experiments/industrialization-village-two-stage-propagated.R
# Smoke: PROP_SMOKE=1 (1 spec, M = 5, short chains).
###############################################################################

library(here)
source(here("code", "_startup", "init.R"))
suppressMessages({ library(brms); library(dplyr) })

make_out <- function(name) {
  p <- here("outputs", "_experiments", name)
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)
  p
}
hpdi <- function(x, m = 0.95) { x <- sort(x); n <- length(x); k <- floor(m * n)
  i <- which.min(x[(k + 1):n] - x[1:(n - k)]); c(x[i], x[i + k]) }

SMOKE   <- nzchar(Sys.getenv("PROP_SMOKE"))
specs   <- if (SMOKE) "sos-urb" else c("sos-urb", "steps-urb", "enmo-urb")
M_IMP   <- if (SMOKE) 5 else 50
imp_cfg <- if (SMOKE) list(w = 300, i = 800, c = 2) else list(w = 500, i = 1000, c = 2)
jnt_cfg <- if (SMOKE) list(w = 300, i = 800, c = 2) else list(w = 1000, i = 3000, c = 4)

out_dir <- make_out("industrialization-village-two-stage-propagated")
src_dir <- here("outputs", "_experiments", "industrialization-village-two-stage")

# summaries of an AERF draws matrix (draws x grid), on the natural outcome scale
summarize_aerf <- function(ep, grid, method) {
  ng    <- length(grid)
  endpt <- ep[, ng] - ep[, 1]
  swing <- apply(ep, 1, function(r) max(r) - min(r))
  he    <- hpdi(endpt); hs <- hpdi(swing)
  data.frame(method = method,
             endpt_lo = he[1], endpt_hi = he[2], p_decline = mean(endpt < 0),
             swing_lo = hs[1], swing_hi = hs[2],
             endpt_width = he[2] - he[1])
}

all_res <- list()

for (key in specs) {
  spec <- model_templates[[key]]
  sf   <- spec$outcome_scale_factor
  resp <- all.vars(spec$bf$formula)[1]
  d    <- prep_local_data(spec, dat)
  d    <- d[!is.na(d$village_id) & !is.na(d$industrial_index), ]
  d$village_id <- droplevels(factor(d$village_id))
  villages <- levels(d$village_id)


  # ---- Stage 1: reuse the saved fit when available ----

  s1_path <- file.path(src_dir, paste0("stage1-fit_", key, ".rds"))
  if (file.exists(s1_path)) {
    message(sprintf("Propagated two-stage, stage 1: %s  (reusing %s)", key, basename(s1_path)))
    m1 <- readRDS(s1_path)
  } else {
    message(sprintf("Propagated two-stage, stage 1: %s  (refitting; n=%d, %d villages)",
                    key, nrow(d), length(villages)))
    m1 <- brm(bf(as.formula(sprintf("%s ~ village_id + s(age_years, k = 5) + sex", resp)),
                 family = gaussian(link = "identity")),
              data = d, prior = set_prior("student_t(3, 0, 2.5)", class = "sds"),
              warmup = WARMUP, iter = ITER, thin = THIN, chains = CHAINS, cores = CHAINS,
              seed = SEED, backend = "cmdstanr", control = BRMS_CONTROL, refresh = 0, silent = 2)
  }

  # standardized village means, one column per village, one row per stage-1 draw
  cohort   <- d[, c("age_years", "sex")]
  ndraws   <- brms::ndraws(m1)
  adj_draws <- matrix(NA_real_, nrow = ndraws, ncol = length(villages))
  for (k in seq_along(villages)) {
    x <- cohort; x$village_id <- factor(villages[k], levels = villages)
    adj_draws[, k] <- rowMeans(posterior_epred(m1, newdata = x))
  }
  adj_draws <- adj_draws * sf                                    # natural units
  index_v <- vapply(villages, function(v) d$industrial_index[d$village_id == v][1], numeric(1))

  # what the plug-in's independence assumption discards
  cor_s1  <- cor(adj_draws)
  offdiag <- cor_s1[upper.tri(cor_s1)]
  message(sprintf("  stage-1 village-mean posterior correlation: mean %.3f, range [%.3f, %.3f]",
                  mean(offdiag), min(offdiag), max(offdiag)))

  grid <- seq(quantile(index_v, .01), quantile(index_v, .99), length.out = 41)
  res  <- list()


  # ---- cut: multiple imputation over stage-one draws ----

  imp_rows <- sample.int(ndraws, M_IMP)
  imp_data <- lapply(imp_rows, function(r)
    data.frame(adj_mean = adj_draws[r, ], industrial_index = index_v))

  fit_cut <- brm_multiple(bf(adj_mean ~ s(industrial_index, k = 4)),
                          data = imp_data, family = gaussian(), combine = TRUE,
                          warmup = imp_cfg$w, iter = imp_cfg$i,
                          chains = imp_cfg$c, cores = imp_cfg$c,
                          seed = SEED, backend = "cmdstanr",
                          control = list(adapt_delta = 0.99), refresh = 0, silent = 2)
  # Rhat from brm_multiple compares across imputations by design; per-imputation
  # convergence is what matters, so it is not used as a diagnostic here.
  ep_cut <- brms::posterior_epred(fit_cut, newdata = data.frame(industrial_index = grid))
  res$cut <- summarize_aerf(ep_cut, grid, "cut")


  # ---- joint: exact one-stage hierarchical GAM ----

  j_form <- bf(as.formula(sprintf(
    "%s ~ s(age_years, k = 5) + sex + s(industrial_index, k = 4) + (1 | village_id)", resp)),
    family = gaussian(link = "identity"))
  fit_joint <- brm(j_form, data = d,
                   prior = set_prior("student_t(3, 0, 2.5)", class = "sds"),
                   warmup = jnt_cfg$w, iter = jnt_cfg$i, chains = jnt_cfg$c, cores = jnt_cfg$c,
                   seed = SEED, backend = "cmdstanr", control = BRMS_CONTROL,
                   refresh = 0, silent = 2)
  # AERF by G-computation over the cohort's age/sex distribution, marginal of village
  ep_joint <- vapply(grid, function(g) {
    x <- cohort; x$industrial_index <- g
    rowMeans(posterior_epred(fit_joint, newdata = x, re_formula = NA)) * sf
  }, numeric(brms::ndraws(fit_joint)))
  res$joint <- summarize_aerf(ep_joint, grid, "joint")


  # ---- plugin: read back the reported estimator ----

  s2_path <- file.path(src_dir, paste0("stage2_", key, ".rds"))
  if (file.exists(s2_path)) {
    s2 <- readRDS(s2_path)
    res$plugin <- summarize_aerf(s2$aerf_draws, s2$grid, "plugin")
  } else {
    message("  no plug-in stage-2 file found; comparison will omit it")
  }

  tab <- do.call(rbind, res[intersect(c("plugin", "cut", "joint"), names(res))])
  tab <- cbind(spec = key, n_villages = length(villages),
               s1_cor_mean = mean(offdiag), tab, row.names = NULL)

  cat(sprintf("\n=== %s  (%s; %d villages) ===\n", key, spec$outcome_label, length(villages)))
  cat(sprintf("  %-8s %14s %14s %11s %12s\n",
              "method", "endpt lo", "endpt hi", "P(decline)", "CI width"))
  for (i in seq_len(nrow(tab)))
    cat(sprintf("  %-8s %14.4g %14.4g %11.3f %12.4g\n", tab$method[i],
                tab$endpt_lo[i], tab$endpt_hi[i], tab$p_decline[i], tab$endpt_width[i]))

  saveRDS(list(table = tab, adj_draws = adj_draws, index_v = index_v, grid = grid,
               ep_cut = ep_cut, ep_joint = ep_joint,
               fit_cut = fit_cut, fit_joint = fit_joint, s1_cor = cor_s1),
          file.path(out_dir, paste0("propagated_", key, ".rds")))
  all_res[[key]] <- tab
}

summary_tab <- do.call(rbind, all_res)
write.csv(summary_tab, file.path(out_dir, "propagation-comparison.csv"), row.names = FALSE)
cat("\nsaved", file.path(out_dir, "propagation-comparison.csv"), "\n")
