###############################################################################
# Manuscript figure assembly. Reads AERF / AMEF posterior draws from the
# targets cache, computes simul-bands + linear-projection per spec, and
# assembles publication-ready PDFs: each cell is an AERF panel stacked over an
# AMEF panel sharing the exposure x-axis, with nested grey simul-bands
# (5/25/50/75/95%) + red-dashed linear-projection overlay.
#
# Outputs:
#   outputs/figures/final/fig-4-pa-bone.pdf      PA -> bone = manuscript Figure 4: 3 outcomes x 2 exposures
#                                                (rows: osteo / CTX / SOS;
#                                                 cols: steps / ENMO).
#                                                Column-major letters A-F.
#   outputs/figures/final/supp-fig-6-pa-bone-conf.pdf   confounder-DAG counterpart of Figure 4.
#   outputs/figures/final/supp-fig-7-age-conditional-bands.pdf
#   outputs/tables/spec-summary.csv              Headline numbers per spec (all 15).
#
# Manuscript Figure 3 (industrialization) is NOT built here; see the section
# below that says where it lives and why it has a single owner.
#
# Population-average summaries only -- the PA -> bone t2(age, exposure) tensor
# smooth is marginalized over the cohort's observed covariate distribution via
# a counterfactual datagrid (see aerf_draws() / amef_draws() in
# code/_startup/pipeline-helpers.R), not evaluated at typical covariate values.
# Age-conditional AMEFs are produced further down this file; the young/old
# age-subset refits are a separate analysis in
# code/_experiments/age-subset-amef.R.
#
# Band construction is NOT pinned here: simul_credible_bands() takes it from the
# BAND_TYPE environment variable via _startup/options.R, so
# `BAND_TYPE=symmetric just <recipe>` rebuilds every panel with the symmetric
# construction for comparison. Each band's provenance travels with it in the
# `interval_type` column.
###############################################################################


source(here::here("code", "_startup", "init.R"))

set.seed(SEED)


# ---- Output dirs ----

final_dir  <- here::here("outputs", "figures", "final")
tables_dir <- here::here("outputs", "tables")
for (d in c(final_dir, tables_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}


# ---- Per-spec panel data: draws, simul-bands AERF + AMEF, linear projection ----

# `spec_panel_data()` (in _startup/pipeline-helpers.R) reads the targets cache
# and derives the bands. It is shared with figure-industrialization.R so the two
# manuscript figures cannot drift apart in how their panels are built.
#
# Derived here rather than as targets because it is cheap from the stored draws
# and only the figure consumers need it. If something else ever wants the
# simul-bands / lin-proj objects, promote them to per-spec tar_map() targets.

spec_keys <- names(model_templates)
summaries <- setNames(lapply(spec_keys, spec_panel_data), spec_keys)
cat("Built panel data (draws + simul-bands + linear projection) for",
    length(spec_keys), "specs\n")


# ---- PA -> bone manuscript figure (3 outcomes x 2 exposures) ----

# Visual layout, column-major letters:
#
#   Col 1: steps          Col 2: ENMO
#   ┌─────────────┐       ┌─────────────┐
#   │ A osteo     │       │ D osteo     │
#   ├─────────────┤       ├─────────────┤
#   │ B CTX       │       │ E CTX       │
#   ├─────────────┤       ├─────────────┤
#   │ C SOS       │       │ F SOS       │
#   └─────────────┘       └─────────────┘
#
# patchwork::wrap_plots(panels, ncol = 2) fills row-major, so the panel
# list order is: [A, D, B, E, C, F].
#
# Shared axes: the x-axis is unified per column (steps for col 1, ENMO for
# col 2) by pre-computing the union range across the 3 specs in each column and
# passing it as `x_limits` to every cell. Only the bottom row (C, F) shows
# x-axis text + title; rows 1+2 hide both via the AMEF panel's
# `show_x_axis = FALSE`. Y-axis stays per-row because outcomes differ.

# Compute per-exposure-column union x-range across the 3 specs.
union_x_range <- function(keys) {
  rng <- lapply(keys, function(k) {
    bands <- summaries[[k]]$simul_bands_aerf
    range(bands[[summaries[[k]]$spec$exposure]], na.rm = TRUE)
  })
  c(min(vapply(rng, `[`, numeric(1), 1)),
    max(vapply(rng, `[`, numeric(1), 2)))
}

x_range_steps <- union_x_range(c("sos-steps", "ctx-steps", "osteo-steps"))
x_range_enmo  <- union_x_range(c("sos-enmo",  "ctx-enmo",  "osteo-enmo"))

set_1_cells <- list(
  list(key = "osteo-steps", tag = "A", col = 1, row = 1),
  list(key = "osteo-enmo",  tag = "D", col = 2, row = 1),
  list(key = "ctx-steps",   tag = "B", col = 1, row = 2),
  list(key = "ctx-enmo",    tag = "E", col = 2, row = 2),
  list(key = "sos-steps",   tag = "C", col = 1, row = 3),
  list(key = "sos-enmo",    tag = "F", col = 2, row = 3)
)

set_1_panels <- lapply(set_1_cells, function(cell) {
  s         <- summaries[[cell$key]]
  x_limits  <- if (cell$col == 1) x_range_steps else x_range_enmo
  is_bottom <- cell$row == 3
  aerf <- make_aerf_panel(s$spec, s$simul_bands_aerf, s$lin_proj,
                          x_limits = x_limits)
  amef <- make_amef_panel(s$spec, s$simul_bands_amef, s$lin_proj,
                          slope_draws = s$slope_draws,
                          x_limits = x_limits,
                          show_x_axis = is_bottom)
  stack_subfig(aerf, amef, tag = cell$tag)
})

fig_1 <- patchwork::wrap_plots(set_1_panels, ncol = 2)
ggplot2::ggsave(fig_1, file = file.path(final_dir, "fig-4-pa-bone.pdf"),
       height = 8.8, width = 5.4)
cat("Saved fig-4-pa-bone.pdf\n")


# ---- Industrialization manuscript figure: NOT BUILT HERE ----
#
# Figure 3 is written by `code/_final/figure-industrialization.R`, which is its
# single owner. Two generators writing one artifact is a defect in its own
# right: whichever ran last would win, and any disagreement between them --
# panel order, band construction, how the AMEF is computed -- would surface as
# a figure that silently changed depending on the order the scripts were run.
# The manuscript caption fixes the panel order as "(A) average daily step
# count, (B) mean daily ENMO, and (C) tibial speed of sound", and that file is
# where it is implemented.


# ---- Supplementary figure: PA -> bone under the CONFOUNDER DAG (6 specs) ----

# Mirrors Fig 4's 3x2 grid layout, but uses the *-conf spec keys. The
# adjustment set adds fat_mass_kg_z + fat_free_mass_kg_z, operationalizing the
# bundled "Fat mass & lean body mass" DAG node as two separate regression
# covariates. Confounder-DAG fits use the same prior structure as the
# mediator-DAG fits plus the two body-composition coefficients. Per-column
# x-range unified the same way Fig 4 does; only the bottom row carries the
# x-axis text.

x_range_steps_conf <- union_x_range(c("sos-steps-conf", "ctx-steps-conf", "osteo-steps-conf"))
x_range_enmo_conf  <- union_x_range(c("sos-enmo-conf",  "ctx-enmo-conf",  "osteo-enmo-conf"))

set_1_conf_cells <- list(
  list(key = "osteo-steps-conf", tag = "A", col = 1, row = 1),
  list(key = "osteo-enmo-conf",  tag = "D", col = 2, row = 1),
  list(key = "ctx-steps-conf",   tag = "B", col = 1, row = 2),
  list(key = "ctx-enmo-conf",    tag = "E", col = 2, row = 2),
  list(key = "sos-steps-conf",   tag = "C", col = 1, row = 3),
  list(key = "sos-enmo-conf",    tag = "F", col = 2, row = 3)
)

set_1_conf_panels <- lapply(set_1_conf_cells, function(cell) {
  s         <- summaries[[cell$key]]
  x_limits  <- if (cell$col == 1) x_range_steps_conf else x_range_enmo_conf
  is_bottom <- cell$row == 3
  aerf <- make_aerf_panel(s$spec, s$simul_bands_aerf, s$lin_proj,
                          x_limits = x_limits)
  amef <- make_amef_panel(s$spec, s$simul_bands_amef, s$lin_proj,
                          slope_draws = s$slope_draws,
                          x_limits = x_limits,
                          show_x_axis = is_bottom)
  stack_subfig(aerf, amef, tag = cell$tag)
})

fig_supp_conf <- patchwork::wrap_plots(set_1_conf_panels, ncol = 2)
ggplot2::ggsave(fig_supp_conf,
       file   = file.path(final_dir, "supp-fig-6-pa-bone-conf.pdf"),
       height = 8.8, width = 5.4)
cat("Saved supp-fig-6-pa-bone-conf.pdf\n")


# ---- Age-conditional AMEF supplementary figure ----
#
# Mirrors the 3x2 Fig 4 grid (mediator-DAG PA -> bone), but each cell shows the
# AMEF at four ages (25 / 35 / 50 / 65) rather than marginalized over the
# cohort's age distribution. The question it answers is whether the activity
# effect on bone is larger before peak bone mass (~age 35), and whether the
# older-adult share of the cohort could be masking one. The
# t2(age_years, exposure) tensor smooth already encodes an age-varying exposure
# effect; this figure slices it.

age_slices <- c(25, 35, 50, 65)

age_cond_amef_for_key <- function(spec_key) {
  spec <- model_templates[[spec_key]]
  fit  <- spec_draws(spec_key, "fit")
  do.call(rbind, lapply(age_slices, function(a) {
    amef_at_age(fit, spec, age_value = a, dat_raw = get(spec$data))
  }))
}

cat("Computing age-conditional AMEFs (6 specs x 4 ages = 24 avg_slopes calls)...\n")
age_cond_amef_data <- setNames(
  lapply(vapply(set_1_cells, function(c) c$key, character(1)), age_cond_amef_for_key),
  vapply(set_1_cells, function(c) c$key, character(1))
)

# Age-conditional AMEFs with 95% simultaneous HPDI bands and no posterior
# median drawn, matching the convention of every other reported figure. The
# ribbons are what makes the figure readable: they separate a real
# age-conditional divergence from tensor-smooth extrapolation in a sparse
# corner of the age x exposure plane.

age_cond_panels_with_bands <- lapply(seq_along(set_1_cells), function(i) {
  cell <- set_1_cells[[i]]
  spec <- model_templates[[cell$key]]
  x_limits  <- if (cell$col == 1) x_range_steps else x_range_enmo
  is_bottom <- cell$row == 3
  show_leg  <- (cell$row == 1 && cell$col == 2)
  p <- make_age_cond_amef_panel_with_bands(
    spec, age_cond_amef_data[[cell$key]],
    x_limits = x_limits, show_x_axis = is_bottom, show_legend = show_leg
  )
  p + ggplot2::labs(tag = cell$tag) +
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "plain", size = 12))
})

fig_supp_age_bands <- patchwork::wrap_plots(age_cond_panels_with_bands, ncol = 2)
ggplot2::ggsave(fig_supp_age_bands,
       file   = file.path(final_dir, "supp-fig-7-age-conditional-bands.pdf"),
       height = 5.4, width = 6.4)
cat("Saved supp-fig-7-age-conditional-bands.pdf\n")


# ---- Headline summary table per spec ----

# `summarize_one_fit` consumes curvature_draws + an explicit n_analytic
# (computed from prep_local_data so we don't have to load the full fit
# objects, which are ~5 MB each).
curvature_for_summary <- lapply(spec_keys, function(k) spec_draws(k, "curvature"))
names(curvature_for_summary) <- spec_keys

n_analytic_per_spec <- vapply(spec_keys, function(k) {
  spec <- model_templates[[k]]
  nrow(prep_local_data(spec, get(spec$data)))
}, integer(1))

summary_tbl <- do.call(rbind, lapply(spec_keys, function(k) {
  s <- summaries[[k]]
  summarize_one_fit(k, s$slope_draws, curvature_for_summary[[k]],
                    n_analytic = n_analytic_per_spec[[k]])
}))
write.csv(summary_tbl,
          file = file.path(tables_dir, "spec-summary.csv"),
          row.names = FALSE)
cat("Saved spec-summary.csv\n")
print(summary_tbl)

cat("\nDone: _final/figures.R\n")
