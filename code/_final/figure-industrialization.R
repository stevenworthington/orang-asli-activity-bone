###############################################################################
# Manuscript Figure 3: industrialization -> {daily steps, ENMO, tibial SOS}.
#
# THIS FILE IS THE SOLE OWNER OF fig-3-urb.pdf. Panel order is A = daily steps,
# B = ENMO, C = tibial SOS, which is what the manuscript caption states.
#
# Rendered from the primary estimator: the one-stage hierarchical model with
# ADDITIVE smooths s(age, k=5) + s(index, k=4), sex, and a community random
# intercept (1 | village_id), Gaussian identity link. The AERF is taken marginal
# of community (re_formula = NA on these specs), so it is the across-community
# effect.
#
# IDENTICAL MACHINERY TO FIGURE 4 in every respect that could make the two
# disagree: the same `spec_panel_data()` reads the draws and builds the bands,
# the same panel builders draw them, and both read the same targets cache. Only
# the specs differ. In particular the AMEF comes from `avg_slopes` here as it
# does there, never from a finite difference of this figure's own AERF -- see
# `spec_panel_data()` in `_startup/pipeline-helpers.R` for why that distinction
# matters on the 41-point community grid.
#
# The bands are wide because the industrialization index varies only BETWEEN the
# ~25 communities, so that is the number of units carrying the effect however
# many individuals are measured.
#
# The AMEF y-axis is on the reported contrast (per 10 index units), applied by
# make_amef_panel() from spec$contrast_units, so the figure and the
# industrialization slope table quote the same denominator.
#
# Run: Rscript code/_final/figure-industrialization.R
###############################################################################

source(here::here("code", "_startup", "init.R"))
set.seed(SEED)

# Created here rather than assumed: this script is runnable on its own, so it
# cannot rely on figures.R having been run first to make the directory.
final_dir <- here::here("outputs", "figures", "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

# Panel order matches the manuscript caption: (A) daily step count,
# (B) mean daily ENMO, (C) tibial speed of sound.
keys <- c("steps-urb", "enmo-urb", "sos-urb")
tags <- c("A", "B", "C")
S    <- setNames(lapply(keys, spec_panel_data), keys)
cat("Built industrialization panel data for", paste(keys, collapse = ", "), "\n")

panels <- Map(function(key, tag) {
  s    <- S[[key]]
  aerf <- make_aerf_panel(s$spec, s$simul_bands_aerf, s$lin_proj)
  amef <- make_amef_panel(s$spec, s$simul_bands_amef, s$lin_proj,
                          slope_draws = s$slope_draws)
  stack_subfig(aerf, amef, tag = tag)
}, keys, tags)

fig <- patchwork::wrap_plots(panels, ncol = 3)
ggplot2::ggsave(fig, file = file.path(final_dir, "fig-3-urb.pdf"),
                height = 2.93, width = 8.1)
cat("Saved fig-3-urb.pdf\n")
