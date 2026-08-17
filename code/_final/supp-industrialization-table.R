###############################################################################
# Supplementary slope table for the three industrialization analyses
# (Figure 3): linear-projection slope (95% HPDI, per 10 index units),
# P(decline), linearity threshold, and flatness threshold.
#
# HPDI, not SHPDI. Simultaneity is a property of a band over a whole curve; the
# linear-projection slope is one scalar per posterior draw, so its interval is an
# ordinary highest-posterior-density interval. The AERF and AMEF bands in the
# figures are the simultaneous ones.
#
# Reads the targets cache through spec_draws(), the same path Figure 3 and the
# PA -> bone slope table use, and takes the AMEF from the `slope_draws` target
# -- marginaleffects::avg_slopes -- exactly as supp-slope-table.R does. That
# matters for one column in particular: the flatness threshold appears in both
# supplementary tables under the same name, so it must be produced by the same
# derivative estimator in both, and every consumer must be reading the same
# fits.
#
# No point estimates are emitted (HPDI bounds and posterior probabilities only).
#
# The reported contrast comes from spec$contrast_units in the registry rather
# than a hard-coded multiplier, so it cannot drift from the contrast the fits
# were scaled for.
#
# This script PRINTS its numbers rather than writing a file -- they are
# transcribed into the supplement by hand.
###############################################################################


source(here::here("code", "_startup", "init.R"))

outcol <- c("steps-urb" = "ad_tot_step_count_0_24hr",
            "enmo-urb"  = "ad_mean_enmo_mg_0_24hr",
            "sos-urb"   = "tibia_sos")
lab    <- c("steps-urb" = "Average daily step count",
            "enmo-urb"  = "Mean daily ENMO",
            "sos-urb"   = "Tibial speed of sound")

for (key in c("steps-urb", "enmo-urb", "sos-urb")) {
  spec <- model_templates[[key]]; ex <- spec$exposure
  cu   <- spec$contrast_units          # 10 index units; from the registry, not hard-coded

  pred_draws  <- spec_draws(key, "pred")
  slope_draws <- spec_draws(key, "slope")
  nv <- nlevels(droplevels(prep_local_data(spec, get(spec$data))$village_id))

  # --- slope per 10 index units (HPDI bounds only) + P(decline), via linear_projection
  #     (the canonical method used for the main-text numbers) ---
  lp         <- linear_projection(pred_draws, exposure = !!rlang::sym(ex), level = 0.95)
  slope10_lo <- lp$beta_hpdi$lo * cu
  slope10_hi <- lp$beta_hpdi$hi * cu
  p_decline  <- mean(lp$beta_draws$beta < 0)

  # --- linearity: q95 max |deviation from best linear projection| / outcome range ---
  md <- pred_draws |> dplyr::filter(!is.na(draw)) |> dplyr::group_by(drawid) |>
    dplyr::summarize(m = { x <- .data[[ex]]; xc <- x - mean(x)
                   b <- sum(xc * (draw - mean(draw))) / sum(xc^2)
                   a <- mean(draw) - b * mean(x); max(abs(draw - (a + b * x))) },
              .groups = "drop") |> dplyr::pull(m)
  lin_pct <- unname(quantile(md, 0.95)) / diff(range(dat[[outcol[[key]]]], na.rm = TRUE)) * 100

  # --- flatness: q95 of per-draw max |AMEF|, expressed per 10 index units.
  #     Same expression as supp-slope-table.R, on the same kind of draws. ---
  max_abs <- slope_draws |> dplyr::filter(!is.na(draw)) |> dplyr::group_by(drawid) |>
    dplyr::summarize(ma = max(abs(draw)), .groups = "drop") |> dplyr::pull(ma)
  flat <- unname(quantile(max_abs, 0.95)) * cu

  cat(sprintf("\n%s  (%s)\n", lab[[key]], key))
  cat(sprintf("  n communities                      : %d\n", nv))
  cat(sprintf("  slope 95%% HPDI  %-19s: [%.4g, %.4g]\n", spec$contrast_label, slope10_lo, slope10_hi))
  cat(sprintf("  P(decline)                         : %.3f\n", p_decline))
  cat(sprintf("  linearity threshold (%% of range)   : %.1f%%\n", lin_pct))
  cat(sprintf("  flatness threshold %-16s: %.4g\n", spec$contrast_label, flat))
}
