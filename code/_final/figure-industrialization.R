###############################################################################
# Manuscript Figure 3 (industrialization -> {SOS, ENMO, steps}) rendered from the
# VILLAGE-LEVEL two-stage estimator (Supplementary Material 9), NOT the
# individual-level GAM. Reuses the same panel machinery as Figure 4 (make_aerf_panel
# / make_amef_panel / stack_subfig) so the styling matches exactly; only the source
# draws differ -- here the second-stage AERF posterior (one curve per ~25-29 village
# means, smoothed over the industrialization index) instead of the n~860 individual
# GAM. The wide bands and weak/uncertain shapes are the honest village-level picture.
#
# Reads outputs/_experiments/industrialization-village-two-stage/stage2_<key>.rds (aerf_draws + grid),
# builds pred_draws / slope_draws in the pipeline's long format, computes the simul
# bands + linear projection, and writes outputs/figures/final/fig-3-urb.pdf.
#
# Run: script -q /dev/null Rscript code/_final/figure-industrialization.R < /dev/null
###############################################################################

library(here)
source(here("code", "_startup", "init.R"))
suppressMessages(library(dplyr))
set.seed(SEED)

village_dir       <- here("outputs", "_experiments", "industrialization-village-two-stage")
final_dir    <- here("outputs", "figures", "final")
BANDS_LEVELS <- c(0.05, 0.25, 0.50, 0.75, 0.95)

# central finite-difference derivative of an AERF draw-matrix over the grid
fd_deriv <- function(ad, grid) {
  ng <- ncol(ad); d <- matrix(NA_real_, nrow(ad), ng)
  d[, 1]  <- (ad[, 2]  - ad[, 1])    / (grid[2]  - grid[1])
  d[, ng] <- (ad[, ng] - ad[, ng - 1]) / (grid[ng] - grid[ng - 1])
  for (g in 2:(ng - 1)) d[, g] <- (ad[, g + 1] - ad[, g - 1]) / (grid[g + 1] - grid[g - 1])
  d
}

build_oh_summary <- function(key) {
  spec <- model_templates[[key]]
  ex   <- spec$exposure                                   # "industrial_index"
  stage2   <- readRDS(file.path(village_dir, paste0("stage2_", key, ".rds")))
  ad   <- stage2$aerf_draws                                   # draws x grid, natural units
  grid <- stage2$grid
  nd   <- nrow(ad); ng <- ncol(ad)

  long <- function(mat) {
    df <- tibble::tibble(drawid = rep(seq_len(nd), times = ng),
                         draw   = as.vector(mat))
    df[[ex]] <- rep(grid, each = nd)
    df
  }
  pred_draws  <- long(ad)
  slope_draws <- long(fd_deriv(ad, grid))

  list(
    spec             = spec,
    pred_draws       = pred_draws,
    slope_draws      = slope_draws,
    simul_bands_aerf = simul_credible_bands(pred_draws,  exposure = !!rlang::sym(ex),
                                            levels = BANDS_LEVELS, function_type = "AERF", interval_type = "HPDI"),
    simul_bands_amef = simul_credible_bands(slope_draws, exposure = !!rlang::sym(ex),
                                            levels = BANDS_LEVELS, function_type = "AMEF", interval_type = "HPDI"),
    lin_proj         = linear_projection(pred_draws, exposure = !!rlang::sym(ex), level = 0.95)
  )
}

keys <- c("steps-urb", "enmo-urb", "sos-urb")          # citation order: steps (3A, strongest), ENMO (3B), SOS (3C)
tags <- c("A", "B", "C")
S    <- setNames(lapply(keys, build_oh_summary), keys)
cat("Built village two-stage summaries for", paste(keys, collapse = ", "), "\n")

panels <- Map(function(key, tag) {
  s    <- S[[key]]
  aerf <- make_aerf_panel(s$spec, s$simul_bands_aerf, s$lin_proj)
  amef <- make_amef_panel(s$spec, s$simul_bands_amef, s$lin_proj, slope_draws = s$slope_draws)
  stack_subfig(aerf, amef, tag = tag)
}, keys, tags)

fig <- patchwork::wrap_plots(panels, ncol = 3)
ggsave(fig, file = file.path(final_dir, "fig-3-urb.pdf"), height = 2.93, width = 8.1)
cat("Saved fig-3-urb.pdf (village-level two-stage)\n")
