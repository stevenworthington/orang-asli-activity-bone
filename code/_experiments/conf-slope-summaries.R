###############################################################################
# Confounder-DAG linear-projection summaries for the six physical-activity ->
# bone `*-conf` specs, the sensitivity counterpart to the mediator-DAG numbers
# in `pa-contrast-effects.R`.
#
# Uses the SAME functionals as the main-text reporting, by calling the same
# code rather than reproducing it:
#   - `linear_projection()` for the per-draw OLS slope of the AERF on the
#     exposure, and its 95% HPDI -- the same function `_final/supp-slope-table.R`
#     uses, so the two cannot disagree about what a slope is
#   - per-person contrast effect = slope x `spec$contrast_units`, taken from the
#     registry rather than hard-coded, so a rescaling cannot silently change
#     what the reported number means
#   - U95 = 95th percentile of |effect|, the "rule out effects larger than"
#     bound, reported as % of the outcome SD and as a multiple of the
#     measurement-resolution floor (both from `_reference-scales.R`)
#
# U95 is a plain quantile, NOT an interval bound -- it answers "how large could
# the magnitude be" and has no lower end.
#
# Posterior probability of the HYPOTHESIZED direction:
#   - osteocalcin (formation) and tibial SOS (density): hypothesis positive,
#     P(hyp) = 1 - P(slope < 0)
#   - CTX-1 (resorption): hypothesis negative, P(hyp) = P(slope < 0)
#
# Read-only; prints numbers, writes nothing.
###############################################################################


source(here::here("code", "_startup", "init.R"))
source(here::here("code", "_experiments", "_reference-scales.R"))

specs <- tibble::tribble(
  ~key,               ~hyp_pos,
  "osteo-steps-conf", TRUE,
  "ctx-steps-conf",   FALSE,
  "sos-steps-conf",   TRUE,
  "osteo-enmo-conf",  TRUE,
  "ctx-enmo-conf",    FALSE,
  "sos-enmo-conf",    TRUE
)

for (i in seq_len(nrow(specs))) {
  key     <- specs$key[i]
  hyp_pos <- specs$hyp_pos[i]
  spec    <- model_templates[[key]]
  ex      <- spec$exposure
  cu      <- spec$contrast_units      # registry, not a hard-coded width
  oc      <- outcome_key(key)
  sd_out  <- analytic_sd(key)
  floor_o <- NOISE_FLOOR[[oc]]
  unit    <- UNIT_LABEL[[oc]]

  # The single read path out of the targets cache, so this script and the
  # figures it sits beside cannot be looking at different fits.
  pred <- spec_draws(key, "pred")

  lp   <- linear_projection(pred, exposure = !!rlang::sym(ex), level = 0.95)
  beta <- lp$beta_draws$beta
  eff  <- beta * cu

  p_neg  <- mean(beta < 0)                         # scale-invariant
  p_hyp  <- if (hyp_pos) 1 - p_neg else p_neg
  ch     <- ggdist::hdci(eff, .width = 0.95)
  u95    <- as.numeric(stats::quantile(abs(eff), 0.95))

  cat(sprintf("\n=== %-16s (%s, %s) ===\n", key, unit, spec$contrast_label))
  cat(sprintf("  P(hypothesized direction) = %.3f   [P(slope<0)=%.3f]\n", p_hyp, p_neg))
  cat(sprintf("  contrast effect 95%% HPDI  = [%.4g, %.4g] %s\n",
              unname(ch[1]), unname(ch[2]), unit))
  cat(sprintf("  U95 (95th pct |effect|)   = %.4g %s = %.1f%% of analytic SD = %.1fx resolution floor\n",
              u95, unit, 100 * u95 / sd_out, u95 / floor_o))
}
