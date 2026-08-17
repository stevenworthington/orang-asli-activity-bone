###############################################################################
# Per-person contrast effects for the six mediator-DAG PA -> bone nulls. These
# are the magnitude bounds the main text quotes; `conf-slope-summaries.R` is the
# confounder-DAG counterpart and reports the same quantities the same way.
#
# The effect is the best-linear-projection slope of the AERF (per-draw OLS of
# the AERF on the exposure, matching the figure caption) times the reported
# contrast. Both come from shared code rather than being recomputed here:
# `linear_projection()` for the slope, `spec$contrast_units` for the contrast.
# The contrasts are 5,000 steps/day and 10 mg ENMO, and because data.R divides
# those two exposures by exactly those amounts, one modelled unit IS one
# reported contrast and `contrast_units` is 1.
#
# Reports the signed 95% HPDI, P(hypothesized direction), and
# U95 = 95th percentile of |effect| -- the "rule out effects larger than"
# bound, which is a plain quantile and has no lower end. U95 is then read
# against the two reference scales in `_reference-scales.R`: the analytic-sample SD and
# the measurement-resolution floor.
#
# Read-only; prints numbers, writes nothing.
###############################################################################


source(here::here("code", "_startup", "init.R"))
source(here::here("code", "_experiments", "_reference-scales.R"))

specs <- tibble::tribble(
  ~key,          ~hyp_pos,
  "sos-steps",   TRUE,
  "ctx-steps",   FALSE,
  "osteo-steps", TRUE,
  "sos-enmo",    TRUE,
  "ctx-enmo",    FALSE,
  "osteo-enmo",  TRUE
)

cat(sprintf("%-13s %-16s | proj slope/unit (95%% HPDI) | contrast 95%% HPDI | P(hyp) | U95 |effect| | vs floor | vs SD\n",
            "spec", "contrast"))

for (i in seq_len(nrow(specs))) {
  key     <- specs$key[i]
  hyp_pos <- specs$hyp_pos[i]
  spec    <- model_templates[[key]]
  ex      <- spec$exposure
  cu      <- spec$contrast_units
  oc      <- outcome_key(key)
  sd_out  <- analytic_sd(key)
  floor_o <- NOISE_FLOOR[[oc]]
  unit    <- UNIT_LABEL[[oc]]

  pred <- spec_draws(key, "pred")

  lp   <- linear_projection(pred, exposure = !!rlang::sym(ex), level = 0.95)
  beta <- lp$beta_draws$beta
  eff  <- beta * cu

  p_neg <- mean(beta < 0)
  p_hyp <- if (hyp_pos) 1 - p_neg else p_neg
  sl_h  <- lp$beta_hpdi
  ef_h  <- ggdist::hdci(eff, .width = 0.95)
  u95   <- as.numeric(stats::quantile(abs(eff), 0.95))

  cat(sprintf("%-13s %-16s | %.4g [%.4g, %.4g] | [%.4g, %.4g] | %.2f | %.4g %s | %.2fx | %.0f%%SD\n",
              key, spec$contrast_label,
              sl_h$mid, sl_h$lo, sl_h$hi,
              unname(ef_h[1]), unname(ef_h[2]),
              p_hyp, u95, unit, u95 / floor_o, 100 * u95 / sd_out))
}
