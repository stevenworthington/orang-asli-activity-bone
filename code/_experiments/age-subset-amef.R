###############################################################################
# Age-subset sensitivity: refit each PA → bone mediator spec separately on
# young (< 35) and old (≥ 35) subsets, then overlay the two subgroup AMEFs
# for comparison. Complement to the smooth-based age-conditional analysis in
# supp-fig-7-age-conditional-bands.pdf, which slices the tensor smooth at fixed
# ages instead of refitting on subsets.
#
# Cutoff rationale: age 35 is the developmental-biology threshold for peak
# bone mass; it's also very close to the cohort's median age (36.5), so this
# is roughly a median split AND a substantively motivated cut. For each of
# the 6 mediator PA-bone specs:
#   - Refit on dat |> filter(age_years <  35)  =>  fit_<spec>_young
#   - Refit on dat |> filter(age_years >= 35)  =>  fit_<spec>_old
#
# Sampling is the production configuration, so this figure carries the same
# posterior draws as every other reported quantity. A reduced run is tempting
# here, since the question the figure asks is whether two subgroup AMEFs look
# different rather than what either one is precisely. It is the wrong economy
# for this particular plot: the band drawn below is a 95% SIMULTANEOUS
# interval, and its endpoints sit at the extreme order statistics of whatever
# posterior sample it is given, so a thin posterior stops estimating them and
# simply reports the sample minimum and maximum.
#
# `SAMPLING_MODE=smoke` gives the fast run for exploratory work, the same
# switch the per-spec scripts honor. Any mode change requires clearing the
# cache below, which the filenames do not encode.
#
# Outputs:
#   outputs/models/age-subset/<spec>-{young,old}.Rdata    (gitignored)
#   outputs/figures/final/supp-fig-8-age-subset.pdf    (committed)
###############################################################################


source(here::here("code", "_startup", "init.R"))


# ---- Config -----------------------------------------------------------------

set.seed(SEED)

AGE_CUTOFF <- 35
SMOKE_MODE <- identical(Sys.getenv("SAMPLING_MODE"), "smoke")
SAMPLING <- if (SMOKE_MODE) {
  list(warmup = SMOKE_WARMUP, iter = SMOKE_ITER, chains = SMOKE_CHAINS, thin = THIN)
} else {
  list(warmup = WARMUP, iter = ITER, chains = CHAINS, thin = THIN)
}

spec_keys <- c("osteo-steps", "ctx-steps",   "sos-steps",
               "osteo-enmo",  "ctx-enmo",    "sos-enmo")

out_dir_models <- here::here("outputs", "models", "age-subset")
out_dir_final  <- here::here("outputs", "figures", "final")
dir.create(out_dir_models, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_final,  recursive = TRUE, showWarnings = FALSE)


# ---- Fit one spec on a given age subset --------------------------------------

fit_subset <- function(spec_key, age_group) {
  spec <- model_templates[[spec_key]]
  out_path <- file.path(out_dir_models, sprintf("%s-%s.Rdata", spec_key, age_group))
  # The filename encodes the spec name and the age group and nothing else, so on
  # its own it cannot tell a production fit from a smoke fit, or a fit of the
  # current registry from one of a superseded registry. Rather than trusting it,
  # the cached object is checked against what would be fitted now.
  dat_local <- prep_local_data(spec, dat)
  dat_sub <- if (age_group == "young") {
    dat_local |> dplyr::filter(age_years < AGE_CUTOFF)
  } else {
    dat_local |> dplyr::filter(age_years >= AGE_CUTOFF)
  }

  pri <- build_priors(spec, dat_sub)
  expected_draws <- SAMPLING$chains * (SAMPLING$iter - SAMPLING$warmup) / SAMPLING$thin

  # A brms object carries the data, formula and priors it was fitted with, so the
  # cache can be checked rather than trusted, and no sidecar file can go missing.
  cache_mismatch <- function(fit) {
    if (!inherits(fit, "brmsfit")) return("not a brmsfit")
    if (!base::setequal(names(fit$data), base::intersect(names(dat_sub), names(fit$data))))
      return("stored data has different columns")
    if (!isTRUE(all.equal(fit$data, as.data.frame(dat_sub)[, names(fit$data), drop = FALSE],
                          check.attributes = FALSE)))
      return("stored data differs from the current analytic subset")
    if (!identical(paste(deparse(stats::formula(fit)$formula), collapse = " "),
                   paste(deparse(spec$bf$formula), collapse = " ")))
      return("stored formula differs from the registry")
    k <- c("prior", "class", "coef", "group")
    if (!isTRUE(all.equal(fit$prior[order(fit$prior$class, fit$prior$coef), k],
                          pri[order(pri$class, pri$coef), k],
                          check.attributes = FALSE)))
      return("stored priors differ from build_priors()")
    if (brms::ndraws(fit) != expected_draws)
      return(sprintf("stored fit has %d draws, expected %d (smoke vs production?)",
                     brms::ndraws(fit), expected_draws))
    NA_character_
  }

  if (file.exists(out_path)) {
    load(out_path)
    why <- cache_mismatch(fit)
    if (is.na(why)) {
      cat(sprintf("  Loading cached %s [%s]...\n", spec_key, age_group))
      return(invisible(fit))
    }
    cat(sprintf("  Cached %s [%s] rejected (%s); refitting.\n", spec_key, age_group, why))
  }

  cat(sprintf("  Fitting %s [%s, n=%d]...\n", spec_key, age_group, nrow(dat_sub)))
  # Priors autoscale to the ANALYTIC data, so they are built from the SUBSET,
  # not the full cohort -- the subset has its own outcome SD and its own
  # design-matrix column SDs.
  fit <- brms::brm(
    spec$bf,
    data      = dat_sub,
    prior     = build_priors(spec, dat_sub),
    warmup    = SAMPLING$warmup,
    iter      = SAMPLING$iter,
    thin      = SAMPLING$thin,
    chains    = SAMPLING$chains,
    cores     = SAMPLING$chains,
    seed      = SEED,
    backend   = "cmdstanr",
    control   = BRMS_CONTROL,
    save_pars = brms::save_pars(all = TRUE),
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
  # spec_grid() on the SUBSET frame: each subgroup gets its own support, which
  # is deliberate (the panels union the two ranges for display below).
  grid <- spec_grid(spec, dat_local)
  # Counterfactual datagrid -- the fit is already restricted to this age
  # subgroup, so `model$data` is the subgroup data and the counterfactual
  # replication averages over the subgroup's observed covariate distribution.
  # Population-average AMEF within the subgroup, the same convention the
  # main-text panels use.
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

cat(sprintf("Refitting 12 age-subset fits in %s mode (%d chains x %d iter, %d warmup, thin %d)...\n",
            if (SMOKE_MODE) "smoke" else "production",
            SAMPLING$chains, SAMPLING$iter, SAMPLING$warmup, SAMPLING$thin))
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

  # 95% SIMULTANEOUS credible bands per subgroup. `interval_type` is left at
  # its default, so these are built the same way as the main-text Fig 3 / Fig 4
  # AMEF panels and the smooth-based age-conditional figure (functions.R
  # `make_age_cond_amef_panel_with_bands`), and follow BAND_TYPE with them.
  summary_df <- amef_subset_draws_df |>
    dplyr::group_split(age_group) |>
    purrr::map_dfr(function(df) {
      bands <- simul_credible_bands(
        df, exposure = !!exp_sym,
        value_col = "draw", levels = 0.95,
        function_type = "AMEF"
      )
      bands$age_group <- df$age_group[1]
      bands
    }) |>
    dplyr::mutate(age_group = factor(age_group,
                                     levels = c("young", "old"),
                                     labels = c("< 35", "≥ 35")))

  # Presentation boundary. The cached draws are per ONE modelled exposure unit,
  # so the reported contrast is applied here -- the same place, and the same
  # way, as the two shared AMEF builders in functions.R. The y-label carries
  # that contrast too: a slope axis reading "Slope of predicted Tibial SOS
  # (m/s)" states no denominator, so the numbers on it mean nothing on their
  # own.
  cu <- spec$contrast_units
  if (!is.numeric(cu) || length(cu) != 1L || !is.finite(cu) || cu <= 0)
    stop("spec$contrast_units must be a positive number; got ", format(cu))
  summary_df <- summary_df |>
    dplyr::mutate(dplyr::across(dplyr::any_of(c("m", "sd", "med", "lo", "hi")),
                                function(v) v * cu))

  y_label <- paste0("Slope of predicted\n", spec$outcome_label,
                    "\n", spec$contrast_label)
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

fig <- patchwork::wrap_plots(panels, ncol = 2)
ggplot2::ggsave(fig,
       file   = file.path(out_dir_final, "supp-fig-8-age-subset.pdf"),
       height = 5.4, width = 6.4,
       device = cairo_pdf)   # cairo_pdf for proper ≥ Unicode rendering
cat("\nSaved outputs/figures/final/supp-fig-8-age-subset.pdf\n")
cat("Done: code/_experiments/age-subset-amef.R\n")
