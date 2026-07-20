###############################################################################
# EXPERIMENT: village fixed effects for the six PA -> bone specs.
#
# Adds `+ village_id` (fixed-effect dummies) to each existing PA -> bone formula and
# re-fits, to estimate the WITHIN-village individual PA effect (absorbing all
# village-constant confounding: industrialization, diet/ecology, per-village
# operator, etc.). This is a *different estimand* from the pooled model, not a
# correction of it. Reports FE vs cached-pooled side-by-side on the same
# summaries (best-linear-projection slope HPDI, P(slope<0), per-person contrast
# HPDI + U95 vs the reconciled noise floor / outcome SD).
#
# ISOLATED: writes only to outputs/_experiments/pa-bone-village-fe/. Does NOT touch the
# registry, the targets store, or any canonical model/figure. Run with the pty
# wrapper + sandbox off (CODING.md). Smoke mode: FE_SMOKE=1 (1 spec, fast MCMC).
###############################################################################

library(here)
source(here("code", "_startup", "init.R"))
suppressMessages({ library(brms); library(dplyr) })

SMOKE <- nzchar(Sys.getenv("FE_SMOKE"))
pa_bone_specs  <- c("sos-steps", "sos-enmo", "ctx-steps", "ctx-enmo", "osteo-steps", "osteo-enmo")
if (SMOKE) pa_bone_specs <- "sos-steps"
cfg <- if (SMOKE) list(w = 500, i = 1500, c = 4) else list(w = WARMUP, i = ITER, c = CHAINS)

out_dir <- here("outputs", "_experiments", "pa-bone-village-fe")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hpdi <- function(x, m = 0.95) { x <- sort(x); n <- length(x); k <- floor(m * n)
  i <- which.min(x[(k + 1):n] - x[1:(n - k)]); c(x[i], x[i + k]) }

floors <- c(sos = 31.6,   ctx = 0.008, osteo = 990)        # reconciled measurement-noise floors
sds    <- c(sos = 182.15, ctx = 0.116, osteo = 19082.9)    # outcome SDs (natural units)
units  <- c(sos = "m/s",  ctx = "ng/mL", osteo = "pg/mL")
cwidth <- function(ex) if (identical(ex, "ad_steps_1k")) 5 else 10   # 5,000 steps / 10 mg ENMO

# best-linear-projection slope + per-person contrast from AERF draws (natural units)
summ <- function(pred, ex, w) {
  pred <- pred[is.finite(pred$draw), ]
  pred <- pred[order(pred$drawid, pred[[ex]]), ]
  x <- sort(unique(pred[[ex]])); ng <- length(x); nd <- nrow(pred) / ng
  M <- matrix(pred$draw, nrow = nd, ncol = ng, byrow = TRUE)
  xc <- x - mean(x); slope <- as.numeric((M %*% xc) / sum(xc^2)); eff <- slope * w
  list(sl = hpdi(slope), pneg = mean(slope < 0), ch = hpdi(eff),
       u95 = as.numeric(quantile(abs(eff), 0.95)))
}

rows <- list()
for (key in pa_bone_specs) {
  spec <- model_templates[[key]]; ex <- spec$exposure; w <- cwidth(ex)
  oc   <- sub("-.*", "", key)                                # sos / ctx / osteo
  fl   <- floors[[oc]]; sd <- sds[[oc]]; u <- units[[oc]]

  d <- prep_local_data(spec, dat)
  d <- d[!is.na(d$village_id), ]; d$village_id <- droplevels(factor(d$village_id))

  fe_bf <- bf(update(spec$bf$formula, . ~ . + village_id), family = spec$bf$family)
  message(sprintf("FE fit: %s  (n=%d, %d villages)", key, nrow(d), nlevels(d$village_id)))
  fit <- brm(fe_bf, data = d, prior = spec$priors,
             warmup = cfg$w, iter = cfg$i, thin = THIN, chains = cfg$c, cores = cfg$c,
             seed = SEED, backend = "cmdstanr", control = BRMS_CONTROL, refresh = 0, silent = 2)

  pred_fe <- aerf_draws(fit, spec, d)              # village-standardized within AERF, natural units
  saveRDS(pred_fe, file.path(out_dir, paste0("aerf_", key, ".rds")))
  s_fe <- summ(pred_fe, ex, w)

  pred_pl <- as.data.frame(qs2::qs_read(
    file.path("_targets", "objects", paste0("pred_draws_", gsub("-", ".", key)))))
  s_pl <- summ(pred_pl, ex, w)

  cat(sprintf("\n=== %s  (%s, contrast %s) ===\n", key, u,
              if (w == 5) "5,000 steps" else "10 mg ENMO"))
  cat(sprintf("  POOLED: slope [%9.4g, %9.4g]  P(<0)=%.2f | contrast [%9.4g, %9.4g] | U95 %8.4g  (%.2fx floor, %2.0f%% SD)\n",
      s_pl$sl[1], s_pl$sl[2], s_pl$pneg, s_pl$ch[1], s_pl$ch[2], s_pl$u95, s_pl$u95 / fl, 100 * s_pl$u95 / sd))
  cat(sprintf("  FE    : slope [%9.4g, %9.4g]  P(<0)=%.2f | contrast [%9.4g, %9.4g] | U95 %8.4g  (%.2fx floor, %2.0f%% SD)\n",
      s_fe$sl[1], s_fe$sl[2], s_fe$pneg, s_fe$ch[1], s_fe$ch[2], s_fe$u95, s_fe$u95 / fl, 100 * s_fe$u95 / sd))

  rows[[key]] <- data.frame(
    spec = key, unit = u,
    pooled_slope_lo = s_pl$sl[1], pooled_slope_hi = s_pl$sl[2], pooled_p_neg = s_pl$pneg, pooled_U95 = s_pl$u95,
    fe_slope_lo = s_fe$sl[1], fe_slope_hi = s_fe$sl[2], fe_p_neg = s_fe$pneg, fe_U95 = s_fe$u95,
    floor = fl, sd = sd)
}
write.csv(do.call(rbind, rows), file.path(out_dir, "fe-vs-pooled.csv"), row.names = FALSE)
cat("\nsaved", file.path(out_dir, "fe-vs-pooled.csv"), "\n")
