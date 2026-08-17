###############################################################################
# Linear-projection slope + linearity threshold for the three industrialization
# analyses, from the AERF draws in the targets cache. Mirrors the PA -> bone
# reporting exactly -- linear_projection + q95 max-deviation-from-linear as a
# percentage of the outcome's range -- so both analysis sets are summarized on
# one framework and their shape diagnostics are directly comparable.
# Read-only; prints numbers, writes nothing.
###############################################################################
source(here::here("code", "_startup", "init.R"))
outcol <- c("sos-urb" = "tibia_sos", "enmo-urb" = "ad_mean_enmo_mg_0_24hr",
            "steps-urb" = "ad_tot_step_count_0_24hr")

for (key in c("steps-urb", "enmo-urb", "sos-urb")) {
  spec <- model_templates[[key]]; ex <- spec$exposure
  pred <- spec_draws(key, "pred")
  grid <- sort(unique(pred[[ex]]))

  lp <- linear_projection(pred, exposure = !!rlang::sym(ex), level = 0.95)
  slo <- lp$beta_hpdi$lo; shi <- lp$beta_hpdi$hi; p_neg <- mean(lp$beta_draws$beta < 0)
  rng <- diff(range(grid))                                   # index span on the prediction grid

  md <- pred |> dplyr::filter(!is.na(draw)) |> dplyr::group_by(drawid) |>
    dplyr::summarize(m = { x <- .data[[ex]]; xc <- x - mean(x)
                   b <- sum(xc * (draw - mean(draw))) / sum(xc^2)
                   a <- mean(draw) - b * mean(x); max(abs(draw - (a + b * x))) },
              .groups = "drop") |> dplyr::pull(m)
  orange  <- diff(range(dat[[outcol[[key]]]], na.rm = TRUE))
  lin_pct <- unname(quantile(md, 0.95)) / orange * 100

  cat(sprintf("\n%s\n  linear-projection slope/index-unit:  [%.4g, %.4g]   P(slope<0)=%.3f\n", key, slo, shi, p_neg))
  cat(sprintf("  implied across-gradient (slope x span):[%.4g, %.4g]\n", slo * rng, shi * rng))
  cat(sprintf("  linearity threshold: %.1f%% of outcome range  (smaller = more linear)\n", lin_pct))
}
