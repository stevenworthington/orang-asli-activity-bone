###############################################################################
# Manuscript figure assembly. Reads AERF / AMEF posterior draws from the
# targets cache, computes simul-bands + linear-projection per spec, and
# assembles publication-ready PDFs in the bone-turnover convention
# (graveyard fig-5 / fig-4 lineage): each cell is an AERF panel stacked
# over an AMEF panel sharing the exposure x-axis, with nested grey simul-
# bands (5/25/50/75/95%) + red-dashed linear-projection overlay.
#
# Outputs:
#   outputs/figures/final/fig-4-pa-bone.pdf      PA -> bone = manuscript Figure 4: 3 outcomes x 2 exposures
#                                                (rows: osteo / CTX / SOS;
#                                                 cols: steps / ENMO).
#                                                Column-major letters A-F.
#   outputs/figures/final/fig-3-urb.pdf          industrialization = manuscript Figure 3: 1 row x 3 outcomes
#                                                (cols: SOS / ENMO / steps).
#                                                Letters A-C.
#   outputs/tables/spec-summary.csv              Headline numbers per spec.
#
# Population-average summaries only -- the t2(age, exposure) tensor smooth
# is marginalized over `datagrid(model)` (age at its mean, other covariates
# at typical values). Age-conditional AMEFs (per Ian 2026-05-25)
# are a separate pass; tracked in TASKS.md.
###############################################################################


library(here)
source(here("code", "_startup", "init.R"))

set.seed(SEED)


# ---- Output dirs ----

final_dir  <- here("outputs", "figures", "final")
tables_dir <- here("outputs", "tables")
for (d in c(final_dir, tables_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}


# ---- Load all 9 specs' draws from the targets cache ----

# tar_map renames hyphens to dots in target names (`enmo-urb` -> `enmo.urb`),
# so the saved target keys are dotted.
load_spec_draws <- function(spec_key) {
  tar_key <- gsub("-", ".", spec_key)
  list(
    pred_draws  = targets::tar_read_raw(paste0("pred_draws_",  tar_key)),
    slope_draws = targets::tar_read_raw(paste0("slope_draws_", tar_key))
  )
}

spec_keys <- names(model_templates)
all_draws <- setNames(lapply(spec_keys, load_spec_draws), spec_keys)
cat("Loaded draws for", length(spec_keys), "specs\n")


# ---- Per-spec summaries: simul-bands AERF + AMEF, linear projection ----

# These are computed inside figures.R rather than as targets, on the principle
# that they're cheap to derive from the stored draws and only this figure
# consumer needs them. If a downstream consumer ever wants the simul-bands /
# lin-proj objects, promote them to per-spec tar_map() targets in _targets.R.

BANDS_LEVELS <- c(0.05, 0.25, 0.50, 0.75, 0.95)

build_spec_summaries <- function(key) {
  spec <- model_templates[[key]]
  d    <- all_draws[[key]]
  list(
    spec             = spec,
    pred_draws       = d$pred_draws,
    slope_draws      = d$slope_draws,
    simul_bands_aerf = simul_credible_bands(
      d$pred_draws,  exposure = !!rlang::sym(spec$exposure),
      levels = BANDS_LEVELS, function_type = "AERF", interval_type = "HPDI"
    ),
    simul_bands_amef = simul_credible_bands(
      d$slope_draws, exposure = !!rlang::sym(spec$exposure),
      levels = BANDS_LEVELS, function_type = "AMEF", interval_type = "HPDI"
    ),
    lin_proj         = linear_projection(
      d$pred_draws,  exposure = !!rlang::sym(spec$exposure),
      level = 0.95
    )
  )
}

summaries <- setNames(lapply(spec_keys, build_spec_summaries), spec_keys)
cat("Built simul-bands + linear projection for all", length(spec_keys), "specs\n")


# ---- PA -> bone manuscript figure (3 outcomes x 2 exposures) ----

# Visual layout (column-major letters per graveyard fig-5 convention):
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
# Shared-axes optimization (per wallace-energetic-constraints 99a9571
# pattern): the x-axis is unified per column (steps for col 1, ENMO for
# col 2) by pre-computing the union range across the 3 specs in each
# column and passing it as `x_limits` to every cell. Only the bottom row
# (C, F) shows x-axis text + title; rows 1+2 hide both via the AMEF
# panel's `show_x_axis = FALSE`. Y-axis stays per-row because outcomes
# differ.

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
ggsave(fig_1, file = file.path(final_dir, "fig-4-pa-bone.pdf"),
       height = 8.8, width = 5.4)
cat("Saved fig-4-pa-bone.pdf\n")


# ---- Industrialization manuscript figure ({SOS, ENMO, steps}) ----

# 1 row x 3 cols landscape. Order per Ian's 2026-05-23 email scope:
# "tibial sos, ENMO, and daily steps".

set_2_cells <- list(
  list(key = "sos-urb",   tag = "A"),
  list(key = "enmo-urb",  tag = "B"),
  list(key = "steps-urb", tag = "C")
)

set_2_panels <- lapply(set_2_cells, function(cell) {
  s    <- summaries[[cell$key]]
  aerf <- make_aerf_panel(s$spec, s$simul_bands_aerf, s$lin_proj)
  amef <- make_amef_panel(s$spec, s$simul_bands_amef, s$lin_proj,
                          slope_draws = s$slope_draws)
  stack_subfig(aerf, amef, tag = cell$tag)
})

fig_2 <- patchwork::wrap_plots(set_2_panels, ncol = 3)
# Sized so each panel matches Fig 1's panel size (Fig 1 is 5.4 x 8.8 for 2 cols
# x 3 rows -> ~2.7" wide x ~2.93" tall per cell). Fig 2 has 3 cols x 1 row, so
# 3 x 2.7" = 8.1" wide and 1 x 2.93" = 2.93" tall.
ggsave(fig_2, file = file.path(final_dir, "fig-3-urb.pdf"),
       height = 2.93, width = 8.1)
cat("Saved fig-3-urb.pdf\n")


# ---- Supplement figure: PA -> bone under the CONFOUNDER DAG (6 specs) ----

# Mirrors Fig 1's 3x2 grid layout, but uses the *-conf spec keys. The
# adjustment set adds fat_mass_kg_z + fat_free_mass_kg_z (operationalizing
# the bundled "Fat mass & lean body mass" DAG node as two separate
# regression covariates per the 2026-05-19 issue 1 decision). Confounder-DAG
# fits use the same prior structure as the mediator-DAG fits plus the two
# body-comp coefficients. Per-column x-range unified the same way Fig 1
# does; only the bottom row carries the x-axis text.

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
ggsave(fig_supp_conf,
       file   = file.path(final_dir, "supp-fig-6-pa-bone-conf.pdf"),
       height = 8.8, width = 5.4)
cat("Saved supp-fig-6-pa-bone-conf.pdf\n")


# ---- Age-conditional AMEF supplement figure (per Ian 2026-05-25) ----
#
# Mirrors the 3x2 Fig 4 grid (mediator-DAG PA -> bone) but each cell shows
# four median AMEF curves color-coded by age (25 / 35 / 50 / 65). Addresses
# Ian's concern that PA effects on bone might be largest pre-peak-bone-mass
# (~age 35) and that the older-adult share of the cohort could mask a real
# effect. The existing t2(age_years, exposure) tensor smooth already encodes
# age-varying exposure effects -- this figure just slices the AMEF at
# meaningful ages rather than marginalizing over the cohort's age distribution
# (which is what Fig 4 does).
#
# Uncertainty intentionally omitted: 4 overlapping ribbons would be unreadable.
# The population-average AMEFs in Fig 4 already carry the credible intervals.

age_slices <- c(25, 35, 50, 65)

age_cond_amef_for_key <- function(spec_key) {
  spec    <- model_templates[[spec_key]]
  tar_key <- gsub("-", ".", spec_key)
  fit     <- targets::tar_read_raw(paste0("fit_", tar_key))
  do.call(rbind, lapply(age_slices, function(a) {
    amef_at_age(fit, spec, age_value = a, dat_raw = get(spec$data))
  }))
}

cat("Computing age-conditional AMEFs (6 specs x 4 ages = 24 avg_slopes calls)...\n")
age_cond_amef_data <- setNames(
  lapply(vapply(set_1_cells, function(c) c$key, character(1)), age_cond_amef_for_key),
  vapply(set_1_cells, function(c) c$key, character(1))
)

# Age-conditional AMEFs WITH 95% HPDI bands (the sole age-conditional figure;
# the medians-only companion was dropped as redundant). Each age slice gets a
# posterior-uncertainty ribbon (alpha = 0.12) to expose where apparent age-
# conditional divergences are real signal vs. tensor-smooth extrapolation.

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
ggsave(fig_supp_age_bands,
       file   = file.path(final_dir, "supp-fig-7-age-conditional-bands.pdf"),
       height = 5.4, width = 6.4)
cat("Saved supp-fig-7-age-conditional-bands.pdf\n")


# ---- Headline summary table per spec ----

# `summarize_one_fit` consumes curvature_draws + an explicit n_analytic
# (computed from prep_local_data so we don't have to load the full fit
# objects, which are ~5 MB each).
curvature_for_summary <- lapply(spec_keys, function(k) {
  tar_key <- gsub("-", ".", k)
  targets::tar_read_raw(paste0("curvature_draws_", tar_key))
})
names(curvature_for_summary) <- spec_keys

n_analytic_per_spec <- vapply(spec_keys, function(k) {
  spec <- model_templates[[k]]
  nrow(prep_local_data(spec, get(spec$data)))
}, integer(1))

summary_tbl <- do.call(rbind, lapply(spec_keys, function(k) {
  s <- summaries[[k]]
  summarize_one_fit(k, s$pred_draws, s$slope_draws, curvature_for_summary[[k]],
                    n_analytic = n_analytic_per_spec[[k]])
}))
write.csv(summary_tbl,
          file = file.path(tables_dir, "spec-summary.csv"),
          row.names = FALSE)
cat("Saved spec-summary.csv\n")
print(summary_tbl)

cat("\nDone: _final/figures.R\n")
