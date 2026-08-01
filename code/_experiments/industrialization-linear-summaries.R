###############################################################################
# Linear-projection slope + linearity threshold for the three industrialization
# analyses, from the village-level two-stage AERF draws. Mirrors the PA->bone reporting
# (linear_projection + q95 max-deviation-from-linear, as % of the outcome's range)
# so the industrialization Results can use the same framework rather than the
# post-hoc interior-peak test. Read-only; prints numbers, edits nothing.
###############################################################################
library(here); source(here("code", "_startup", "init.R")); suppressMessages(library(dplyr))
village_dir <- here("outputs", "_experiments", "industrialization-village-two-stage")
outcol <- c("sos-urb" = "tibia_sos", "enmo-urb" = "ad_mean_enmo_mg_0_24hr",
            "steps-urb" = "ad_tot_step_count_0_24hr")

for (key in c("steps-urb", "enmo-urb", "sos-urb")) {
  spec <- model_templates[[key]]; ex <- spec$exposure
  stage2 <- readRDS(file.path(village_dir, paste0("stage2_", key, ".rds")))
  ad <- stage2$aerf_draws; grid <- stage2$grid; nd <- nrow(ad); ng <- ncol(ad)
  pred <- tibble::tibble(drawid = rep(seq_len(nd), times = ng), draw = as.vector(ad))
  pred[[ex]] <- rep(grid, each = nd)

  lp <- linear_projection(pred, exposure = !!rlang::sym(ex), level = 0.95)
  slo <- lp$beta_hpdi$lo; shi <- lp$beta_hpdi$hi; p_neg <- mean(lp$beta_draws$beta < 0)
  rng <- diff(range(grid))                                   # index span on the prediction grid

  md <- pred |> filter(!is.na(draw)) |> group_by(drawid) |>
    summarize(m = { x <- .data[[ex]]; xc <- x - mean(x)
                   b <- sum(xc * (draw - mean(draw))) / sum(xc^2)
                   a <- mean(draw) - b * mean(x); max(abs(draw - (a + b * x))) },
              .groups = "drop") |> pull(m)
  orange  <- diff(range(dat[[outcol[[key]]]], na.rm = TRUE))
  lin_pct <- unname(quantile(md, 0.95)) / orange * 100

  cat(sprintf("\n%s\n  linear-projection slope/index-unit:  [%.4g, %.4g]   P(slope<0)=%.3f\n", key, slo, shi, p_neg))
  cat(sprintf("  implied across-gradient (slope x span):[%.4g, %.4g]\n", slo * rng, shi * rng))
  cat(sprintf("  linearity threshold: %.1f%% of outcome range  (smaller = more linear)\n", lin_pct))
}
