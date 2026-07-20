###############################################################################
# Supplementary Table 1 numbers: linear-projection slope (95% SHPDI, per 10 index
# units), P(decline), linearity threshold, and flatness threshold for the three
# industrialization analyses (Figure 3), from the Option H two-stage AERF draws
# (outputs/_experiments/set2-option-h/stage2_*.rds).
#
# Flatness is computed from the central finite-difference AMEF, matching
# figure-3-option-h.R, so it is the industrialization analogue of the flatness
# threshold in Supplementary Table 2. No point estimates are emitted (HPDI bounds
# and posterior probabilities only).
###############################################################################


library(here)
source(here("code", "_startup", "init.R"))
suppressMessages(library(dplyr))

oh_dir <- here("outputs", "_experiments", "set2-option-h")
outcol <- c("steps-urb" = "ad_tot_step_count_0_24hr",
            "enmo-urb"  = "ad_mean_enmo_mg_0_24hr",
            "sos-urb"   = "tibia_sos")
lab    <- c("steps-urb" = "Average daily step count",
            "enmo-urb"  = "Mean daily ENMO",
            "sos-urb"   = "Tibial speed of sound")

# central finite-difference AMEF (matches figure-3-option-h.R::fd_deriv)
fd_deriv <- function(ad, grid) {
  ng <- ncol(ad); d <- matrix(NA_real_, nrow(ad), ng)
  d[, 1]  <- (ad[, 2]  - ad[, 1])      / (grid[2]  - grid[1])
  d[, ng] <- (ad[, ng] - ad[, ng - 1]) / (grid[ng] - grid[ng - 1])
  for (g in 2:(ng - 1)) d[, g] <- (ad[, g + 1] - ad[, g - 1]) / (grid[g + 1] - grid[g - 1])
  d
}

for (key in c("steps-urb", "enmo-urb", "sos-urb")) {
  spec <- model_templates[[key]]; ex <- spec$exposure
  oh   <- readRDS(file.path(oh_dir, paste0("stage2_", key, ".rds")))
  ad   <- oh$aerf_draws; grid <- oh$grid; nd <- nrow(ad)
  nv   <- oh$summary$n_villages

  pred <- tibble::tibble(drawid = rep(seq_len(nd), times = ncol(ad)), draw = as.vector(ad))
  pred[[ex]] <- rep(grid, each = nd)

  # --- slope per 10 index units (HPDI bounds only) + P(decline), via linear_projection
  #     (the canonical method used for the main-text numbers) ---
  lp         <- linear_projection(pred, exposure = !!rlang::sym(ex), level = 0.95)
  slope10_lo <- lp$beta_hpdi$lo * 10
  slope10_hi <- lp$beta_hpdi$hi * 10
  p_decline  <- mean(lp$beta_draws$beta < 0)

  # --- linearity: q95 max |deviation from best linear projection| / outcome range ---
  md <- pred |> filter(!is.na(draw)) |> group_by(drawid) |>
    summarize(m = { x <- .data[[ex]]; xc <- x - mean(x)
                   b <- sum(xc * (draw - mean(draw))) / sum(xc^2)
                   a <- mean(draw) - b * mean(x); max(abs(draw - (a + b * x))) },
              .groups = "drop") |> pull(m)
  lin_pct <- unname(quantile(md, 0.95)) / diff(range(dat[[outcol[[key]]]], na.rm = TRUE)) * 100

  # --- flatness: q95 of per-draw max |AMEF|, expressed per 10 index units ---
  amef <- fd_deriv(ad, grid)
  flat <- unname(quantile(apply(abs(amef), 1, max), 0.95)) * 10

  cat(sprintf("\n%s  (%s)\n", lab[[key]], key))
  cat(sprintf("  n villages                         : %d\n", nv))
  cat(sprintf("  slope 95%% SHPDI / 10 index units   : [%.4g, %.4g]\n", slope10_lo, slope10_hi))
  cat(sprintf("  P(decline)                         : %.3f\n", p_decline))
  cat(sprintf("  linearity threshold (%% of range)   : %.1f%%\n", lin_pct))
  cat(sprintf("  flatness threshold / 10 index units: %.4g\n", flat))
}
