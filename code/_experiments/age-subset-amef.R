###############################################################################
# Age-subset sensitivity: refit each PA → bone mediator spec separately on
# young (< 35) and old (≥ 35) subsets, then overlay the two subgroup AMEFs
# for comparison. Complement to the smooth-based age-conditional analysis in
# fig-supp-age-conditional-amef{,-with-bands}.pdf.
#
# Cutoff rationale: age 35 is the developmental-biology threshold for peak
# bone mass; it's also very close to the cohort's median age (36.5), so this
# is roughly a median split AND a substantively motivated cut. For each of
# the 6 mediator PA-bone specs:
#   - Refit on dat |> filter(age_years <  35)  =>  fit_<spec>_young
#   - Refit on dat |> filter(age_years >= 35)  =>  fit_<spec>_old
#
# Smoke mode (4 chains x 1500 iter, 500 warmup) is sufficient for the
# exploratory comparison — the inferential question is whether the two
# subgroups' AMEFs are visibly different from each other or from the
# all-cohort fit, not a high-precision per-subgroup point estimate. Full
# sampling can be a follow-up if the comparison suggests it would help.
#
# Outputs:
#   outputs/models/age-subset/<spec>-{young,old}.Rdata    (gitignored)
#   outputs/figures/final/supp-fig-8-age-subset.pdf    (committed)
###############################################################################


library(here)
source(here("code", "_startup", "init.R"))


# ---- Config -----------------------------------------------------------------

set.seed(SEED)

AGE_CUTOFF <- 35
SAMPLING <- list(
  warmup    = 500,
  iter      = 1500,
  chains    = 4,
  thin      = THIN
)

spec_keys <- c("osteo-steps", "ctx-steps",   "sos-steps",
               "osteo-enmo",  "ctx-enmo",    "sos-enmo")

out_dir_models <- here("outputs", "models", "age-subset")
out_dir_final  <- here("outputs", "figures", "final")
dir.create(out_dir_models, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_final,  recursive = TRUE, showWarnings = FALSE)


# ---- Fit one spec on a given age subset --------------------------------------

fit_subset <- function(spec_key, age_group) {
  spec <- model_templates[[spec_key]]
  out_path <- file.path(out_dir_models, sprintf("%s-%s.Rdata", spec_key, age_group))
  if (file.exists(out_path)) {                    # reuse cached refit (analytic data is stable)
    cat(sprintf("  Loading cached %s [%s]...\n", spec_key, age_group))
    load(out_path); return(invisible(fit))
  }
  dat_local <- prep_local_data(spec, dat)
  dat_sub <- if (age_group == "young") {
    dat_local |> dplyr::filter(age_years < AGE_CUTOFF)
  } else {
    dat_local |> dplyr::filter(age_years >= AGE_CUTOFF)
  }

  cat(sprintf("  Fitting %s [%s, n=%d]...\n", spec_key, age_group, nrow(dat_sub)))
  fit <- brm(
    spec$bf,
    data      = dat_sub,
    prior     = spec$priors,
    warmup    = SAMPLING$warmup,
    iter      = SAMPLING$iter,
    thin      = SAMPLING$thin,
    chains    = SAMPLING$chains,
    cores     = SAMPLING$chains,
    seed      = SEED,
    backend   = "cmdstanr",
    control   = BRMS_CONTROL,
    save_pars = save_pars(all = TRUE),
    refresh   = 0,
    silent    = 2
  )

  save(fit, file = out_path, compress = "gzip")
  invisible(fit)
}


# ---- AMEF draws for a subset fit, tagged with age_group ---------------------

amef_subset_draws <- function(fit, spec_key, age_group) {
  spec <- model_templates[[spec_key]]
  dat_local <- prep_local_data(spec, dat) |>
    dplyr::filter(
      if (age_group == "young") age_years < AGE_CUTOFF else age_years >= AGE_CUTOFF
    )
  grid <- seq(
    quantile(dat_local[[spec$exposure]], spec$grid_quantiles[1], na.rm = TRUE),
    quantile(dat_local[[spec$exposure]], spec$grid_quantiles[2], na.rm = TRUE),
    length.out = 51
  )
  # Counterfactual datagrid -- the fit is already restricted to this age
  # subgroup, so `model$data` is the subgroup data and the counterfactual
  # replication averages over the subgroup's observed covariate distribution.
  # Population-average AMEF within the subgroup. Matches the convention in
  # ~/.config/agents/CODING.md § Postestimation.
  dg_args <- list(model = fit, grid_type = "counterfactual")
  dg_args[[spec$exposure]] <- grid
  new_grid <- do.call(marginaleffects::datagrid, dg_args)

  slopes <- marginaleffects::avg_slopes(
    fit, variables = spec$exposure,
    newdata = new_grid, by = spec$exposure
  ) |>
    marginaleffects::posterior_draws()
  slopes$draw <- slopes$draw * spec$outcome_scale_factor
  slopes$age_group <- age_group
  slopes
}


# ---- Fit all 12 (6 specs x 2 subgroups) -------------------------------------

cat("Refitting 12 age-subset fits in smoke mode...\n")
t_start <- Sys.time()

subset_draws_list <- list()
for (k in spec_keys) {
  for (g in c("young", "old")) {
    fit <- fit_subset(k, g)
    subset_draws_list[[paste(k, g, sep = "-")]] <- amef_subset_draws(fit, k, g)
    rm(fit); gc()
  }
}

cat(sprintf("\nTotal time: %.1f min\n", as.numeric(Sys.time() - t_start, units = "mins")))


# ---- Build comparison figure ------------------------------------------------

# Combine young + old for each spec into one tibble, plot both subgroups on
# the same panel with HPDI ribbons (similar style to the with-bands age-
# conditional figure).

subset_amef <- function(spec_key) {
  dplyr::bind_rows(
    subset_draws_list[[paste(spec_key, "young", sep = "-")]],
    subset_draws_list[[paste(spec_key, "old",   sep = "-")]]
  )
}

set_1_cells <- list(
  list(key = "osteo-steps", tag = "A", col = 1, row = 1),
  list(key = "osteo-enmo",  tag = "D", col = 2, row = 1),
  list(key = "ctx-steps",   tag = "B", col = 1, row = 2),
  list(key = "ctx-enmo",    tag = "E", col = 2, row = 2),
  list(key = "sos-steps",   tag = "C", col = 1, row = 3),
  list(key = "sos-enmo",    tag = "F", col = 2, row = 3)
)

# Compute per-column x-limits using the union of young + old subset grids,
# so the two curves visually align within a column.
union_x_range_subset <- function(keys) {
  rng <- lapply(keys, function(k) {
    spec <- model_templates[[k]]
    d_y <- prep_local_data(spec, dat) |> dplyr::filter(age_years <  AGE_CUTOFF)
    d_o <- prep_local_data(spec, dat) |> dplyr::filter(age_years >= AGE_CUTOFF)
    yq <- quantile(d_y[[spec$exposure]], spec$grid_quantiles, na.rm = TRUE)
    oq <- quantile(d_o[[spec$exposure]], spec$grid_quantiles, na.rm = TRUE)
    c(min(yq[1], oq[1]), max(yq[2], oq[2]))
  })
  c(min(vapply(rng, `[`, numeric(1), 1)),
    max(vapply(rng, `[`, numeric(1), 2)))
}

x_range_steps <- union_x_range_subset(c("sos-steps", "ctx-steps", "osteo-steps"))
x_range_enmo  <- union_x_range_subset(c("sos-enmo",  "ctx-enmo",  "osteo-enmo"))


make_age_subset_panel <- function(spec, amef_subset_draws_df,
                                  x_limits = NULL, show_x_axis = TRUE,
                                  show_legend = FALSE) {
  exp_sym <- rlang::sym(spec$exposure)

  # 95% SIMULTANEOUS credible bands per subgroup (sup-norm calibrated),
  # matching the convention used by main-text Fig 3 / Fig 4 AMEF panels
  # and by the smooth-based age-conditional figure (functions.R
  # `make_age_cond_amef_panel_with_bands`).
  summary_df <- amef_subset_draws_df |>
    dplyr::group_split(age_group) |>
    purrr::map_dfr(function(df) {
      bands <- simul_credible_bands(
        df, exposure = !!exp_sym,
        value_col = "draw", levels = 0.95,
        function_type = "AMEF", interval_type = "HPDI"
      )
      bands$age_group <- df$age_group[1]
      bands
    }) |>
    dplyr::mutate(age_group = factor(age_group,
                                     levels = c("young", "old"),
                                     labels = c("< 35", "≥ 35")))

  y_label <- paste0("Slope of predicted\n", spec$outcome_label)
  x_axis_label <- spec$exposure_label

  grp_colors <- c("< 35" = "#0072B2", "≥ 35" = "#D55E00")  # Okabe-Ito blue + vermillion

  p <- ggplot2::ggplot(summary_df,
                       ggplot2::aes(x = !!exp_sym,
                                    ymin = lo, ymax = hi,
                                    fill = age_group)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
    ggplot2::geom_ribbon(alpha = 0.32, color = NA) +
    ggplot2::scale_fill_manual(name  = "Age group", values = grp_colors) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 1))) +
    x_scale_for(spec$exposure) +
    ggplot2::scale_y_continuous(labels = format_thousands) +
    ggplot2::labs(x = x_axis_label, y = y_label) +
    theme_pub() +
    ggplot2::theme(
      axis.title    = ggplot2::element_text(face = "plain", size = 9),
      legend.position = if (show_legend) "right" else "none"
    )

  if (!show_x_axis) {
    p <- p + ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  }
  if (!is.null(x_limits)) p <- p + ggplot2::coord_cartesian(xlim = x_limits)
  p
}

panels <- lapply(seq_along(set_1_cells), function(i) {
  cell <- set_1_cells[[i]]
  spec <- model_templates[[cell$key]]
  xl <- if (cell$col == 1) x_range_steps else x_range_enmo
  p <- make_age_subset_panel(spec, subset_amef(cell$key),
                             x_limits = xl,
                             show_x_axis = (cell$row == 3),
                             show_legend = (cell$row == 1 && cell$col == 2))
  p + ggplot2::labs(tag = cell$tag) +
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "plain", size = 12))
})

library(patchwork)
fig <- patchwork::wrap_plots(panels, ncol = 2)
ggsave(fig,
       file   = file.path(out_dir_final, "supp-fig-8-age-subset.pdf"),
       height = 5.4, width = 6.4,
       device = cairo_pdf)   # cairo_pdf for proper ≥ Unicode rendering
cat("\nSaved outputs/figures/final/supp-fig-8-age-subset.pdf\n")
cat("Done: code/_experiments/age-subset-amef.R\n")
