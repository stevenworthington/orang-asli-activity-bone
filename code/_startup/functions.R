###############################################################################
# Shared helper functions for the bone-turnover analyses: package loading
# wrapper, publication ggplot themes, scaling helpers, posterior-probability
# summaries, six-panel posterior-predictive checks (overall + grouped), the
# AERF/AMEF/curvature plotting helper, and the flatness / linearity /
# simultaneous-band diagnostics. Sourced by `_startup/init.R` so every analysis
# script gets these helpers for free.
###############################################################################


# ---- Package loading wrapper ----
#
# Strictly LOADS packages (no install). If a package isn't in the project
# library, library() errors and the fix is `renv::restore()`. Installs go
# through renv (via the pak backend), never through this wrapper. See
# ~/.config/agents/CODING.md § R "Package installation and loading".
load_pkgs <- function(...) {
  pkgs <- c(...)
  invisible(lapply(pkgs, library, character.only = TRUE))
}


# ---- ggplot themes ----

theme_pub <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, color = "grey90"),
      strip.background = element_rect(fill = "grey95", color = NA),
      plot.title       = element_text(size = base_size + 1, face = "bold"),
      axis.title       = element_text(size = base_size),
      legend.position  = "none"
    )
}

theme_pub_leg <- function(base_size = 10) {
  theme_pub(base_size) +
    theme(legend.position = "top")
}


# ---- Scaling and shape ----

# z-score (preserves NAs)
scale_this <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Fisher-Pearson skewness (NIST definition)
skew <- function(y) {
  y <- na.omit(y)
  n <- length(y)
  diff <- y - mean(y)
  (sqrt(n - 1) / (n - 2)) * n * (sum(diff^3) / (sum(diff^2)^1.5))
}


# ---- Posterior probability summaries ----

# Posterior threshold probabilities. For a vector of values, returns inside-
# /outside-ROPE probabilities (treating each value as ± a symmetric ROPE
# bound) plus directional P(draws > x) and P(draws < x).
calc_probs <- function(draws, seq_values, rope_lower, rope_upper) {
  tibble(
    values        = seq_values,
    inside_rope   = map_dbl(seq_values, ~ mean(draws >= -abs(.x) & draws <= abs(.x))),
    outside_rope  = map_dbl(seq_values, ~ mean(draws < -abs(.x) | draws > abs(.x))),
    greater_rope  = map_dbl(seq_values, ~ mean(draws > abs(.x))),
    less_rope     = map_dbl(seq_values, ~ mean(draws < -abs(.x))),
    greater_point = map_dbl(seq_values, ~ mean(draws > .x)),
    less_point    = map_dbl(seq_values, ~ mean(draws < .x)),
    rope_lower    = rope_lower,
    rope_upper    = rope_upper,
    prob_in_rope  = mean(draws >= rope_lower & draws <= rope_upper)
  )
}

# Directional probabilities + median + HDI for a single draw vector.
calc_probs_direction <- function(draws) {
  hdi <- HDInterval::hdi(draws)
  tibble(
    prob_positive = mean(draws > 0),
    prob_negative = mean(draws < 0),
    median        = median(draws),
    hdi_lower     = unname(hdi[1]),
    hdi_upper     = unname(hdi[2])
  )
}


# ---- Curvature draws (finite-difference second derivative) ----
#
# Kept in the helper set even though curvature isn't part of the bone-turnover
# figure convention (graveyard fig-5 shows AERF + AMEF only; flatness /
# linearity probabilities carry the headline). The pipeline target
# `curvature_draws_<spec>` continues to produce these in case a downstream
# sanity-check ever wants the second derivative.

compute_curvature_draws <- function(slope_draws, exposure) {
  exp_sym <- rlang::sym(exposure)
  slope_draws |>
    dplyr::arrange(drawid, !!exp_sym) |>
    dplyr::group_by(drawid) |>
    dplyr::mutate(
      x_lag     = dplyr::lag(!!exp_sym),
      x_lead    = dplyr::lead(!!exp_sym),
      draw_lag  = dplyr::lag(draw),
      draw_lead = dplyr::lead(draw),
      curvature = (draw_lead - draw_lag) / (x_lead - x_lag)
    ) |>
    dplyr::filter(!is.na(curvature)) |>
    dplyr::ungroup() |>
    dplyr::transmute(drawid, !!exp_sym := !!exp_sym, draw = curvature)
}


# ---- linear_projection: per-draw best linear approximation of the AERF ----
#
# For each posterior draw of the AERF, fit y = intercept + beta * x and
# return the slope summary plus a band on the fitted line over the exposure
# grid. The band uses sup-norm calibration via simul_credible_bands() so it
# reads as "joint coverage of the entire line at `level`."
#
# Returns a list with:
#   $beta_draws : per-draw slope tibble (drawid, beta)
#   $beta_hpdi  : tibble(lo, mid, hi) with the slope's `level` HPDI
#   $band       : tibble(<exposure>, lo, mid (median), hi, level, ...) — the
#                 projection band for the AERF panel
#
# Used by make_aerf_panel() (band overlay) and make_amef_panel() (slope HPDI
# rendered as red-dashed horizontals + filled rect across the exposure range).

linear_projection <- function(pred_draws, exposure, level = 0.95,
                              value_col = "draw") {

  stopifnot(is.data.frame(pred_draws))

  exposure_sym <- rlang::ensym(exposure)
  exposure_nm  <- rlang::as_name(exposure_sym)

  req_cols <- c("drawid", exposure_nm, value_col)
  if (!all(req_cols %in% names(pred_draws))) {
    stop("pred_draws must contain drawid, exposure, and value columns")
  }
  if (!is.numeric(level) || length(level) != 1 || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("level must be a single number strictly between 0 and 1.")
  }

  x_grid <- pred_draws |>
    dplyr::distinct(.data[[exposure_nm]]) |>
    dplyr::arrange(.data[[exposure_nm]]) |>
    dplyr::pull(.data[[exposure_nm]])

  # Best linear projection of the AERF per posterior draw, computed as the
  # closed-form OLS coefficients (cov(x,y)/var(x)) -- a deterministic linear
  # functional of the Bayesian posterior, NOT a model fit. Numerically identical
  # to lm() coefficients; kept lm-free so no frequentist call appears anywhere.
  lin_summ <- pred_draws |>
    dplyr::filter(!is.na(.data[[value_col]]), !is.na(.data[[exposure_nm]])) |>
    dplyr::group_by(drawid) |>
    dplyr::summarize(
      beta = {
        xc <- .data[[exposure_nm]] - mean(.data[[exposure_nm]])
        sum(xc * (.data[[value_col]] - mean(.data[[value_col]]))) / sum(xc^2)
      },
      intercept = mean(.data[[value_col]]) - beta * mean(.data[[exposure_nm]]),
      .groups   = "drop"
    )

  # draw-level projected line values on the exposure grid
  lin_proj_draws <- lin_summ |>
    dplyr::select(drawid, intercept, beta) |>
    tidyr::crossing(tibble::tibble(!!exposure_sym := x_grid)) |>
    dplyr::mutate(draw = intercept + beta * .data[[exposure_nm]])

  # slope summary (HPDI on beta)
  beta_summary <- lin_summ |>
    ggdist::median_hdi(beta, .width = level) |>
    dplyr::rename(lo = .lower, mid = beta, hi = .upper) |>
    dplyr::select(lo, mid, hi)

  # simul-band on the AERF projection
  band <- simul_credible_bands(
    draws_df      = lin_proj_draws,
    exposure      = !!exposure_sym,
    value_col     = "draw",
    levels        = level,
    function_type = "AERF",
    interval_type = "HPDI"
  )

  list(
    level      = level,
    beta_draws = lin_summ |> dplyr::select(drawid, beta),
    beta_hpdi  = beta_summary,
    band       = band
  )
}


# ---- Nested simul-band ribbons (5/25/50/75/95% in greys) ----
#
# Helper used by make_aerf_panel / make_amef_panel. Takes the long-format
# tibble that `simul_credible_bands(levels = c(0.05, 0.25, 0.50, 0.75, 0.95))`
# returns (one row per exposure × level) and lays five nested geom_ribbons
# in increasing alpha. Returns a list of geom_ribbon layers ready to add to
# a ggplot.

nested_simul_ribbons <- function(bands) {
  list(
    ggplot2::geom_ribbon(data = ~ dplyr::filter(.x, level == 0.95),
                         ggplot2::aes(ymin = lo, ymax = hi),
                         fill = "grey20", alpha = 0.10),
    ggplot2::geom_ribbon(data = ~ dplyr::filter(.x, level == 0.75),
                         ggplot2::aes(ymin = lo, ymax = hi),
                         fill = "grey20", alpha = 0.12),
    ggplot2::geom_ribbon(data = ~ dplyr::filter(.x, level == 0.50),
                         ggplot2::aes(ymin = lo, ymax = hi),
                         fill = "grey20", alpha = 0.14),
    ggplot2::geom_ribbon(data = ~ dplyr::filter(.x, level == 0.25),
                         ggplot2::aes(ymin = lo, ymax = hi),
                         fill = "grey20", alpha = 0.15),
    ggplot2::geom_ribbon(data = ~ dplyr::filter(.x, level == 0.05),
                         ggplot2::aes(ymin = lo, ymax = hi),
                         fill = "grey20", alpha = 0.15)
  )
}


# ---- Panel + stack builders (bone-turnover convention) ----
#
# Match graveyard fig-5: each "cell" is a vertical stack of an AERF panel
# (top) over an AMEF panel (bottom), sharing the exposure x-axis. The AERF
# panel hides its x-axis (carried by the AMEF below). Both panels show
# nested grey simul-bands + a red-dashed linear-projection overlay. The
# AMEF panel additionally shows a zero line and a rug of observed exposure
# values at the bottom.

# Compact axis-label formatter: 12500 -> "12.5k", 5000 -> "5k", 30 -> "30",
# 0.15 -> "0.15". Dynamic decimal precision so tick labels are
# distinguishable. Used as the default y-axis labeller on AERF / AMEF
# panels (and also wired into the x-axis for `steps_1k` exposures).
format_thousands <- function(x) {
  if (!is.numeric(x)) stop("Input must be numeric")
  if (length(x) == 0) return(character(0))

  big <- x[!is.na(x) & abs(x) >= 1000]
  digits <- 1L
  if (length(big) >= 2) {
    diffs_k <- diff(sort(unique(big))) / 1000
    diffs_k <- diffs_k[diffs_k > 0]
    if (length(diffs_k) > 0) {
      digits <- max(0L, as.integer(ceiling(-log10(min(diffs_k)))))
    }
  }

  sapply(x, function(y) {
    if (is.na(y)) return(NA_character_)
    if (abs(y) < 1000) return(as.character(y))
    paste0(formatC(round(y / 1000, digits), format = "fg"), "k")
  })
}

# X-axis tick formatting: steps in /1000 units get a "k" suffix (5k, 10k, ...),
# other exposures use raw numeric.
x_scale_for <- function(exposure_nm) {
  if (grepl("steps_1k", exposure_nm)) {
    ggplot2::scale_x_continuous(labels = function(x) paste0(x, "k"))
  } else {
    ggplot2::scale_x_continuous()
  }
}

make_aerf_panel <- function(spec, simul_bands_aerf, lin_proj,
                            y_label = NULL, x_axis_label = NULL,
                            x_limits = NULL) {
  exp_sym <- rlang::sym(spec$exposure)
  bands   <- simul_bands_aerf |> dplyr::arrange(dplyr::desc(level))
  # Two-line y-axis label (matches graveyard fig-5 convention) so the title
  # fits inside the cell without clipping at 2.6"-wide subfigures.
  if (is.null(y_label)) y_label <- paste0("Predicted\n", spec$outcome_label)

  p <- ggplot2::ggplot(bands, ggplot2::aes(x = !!exp_sym)) +
    nested_simul_ribbons(bands) +
    # red-dashed linear-projection envelope on the AERF
    ggplot2::geom_ribbon(data = lin_proj$band,
                         ggplot2::aes(x = !!exp_sym, ymin = lo, ymax = hi),
                         fill = "red", alpha = 0.10, inherit.aes = FALSE) +
    ggplot2::geom_line(data = lin_proj$band,
                       ggplot2::aes(x = !!exp_sym, y = lo),
                       color = "red", linewidth = 0.5, linetype = "dashed",
                       inherit.aes = FALSE) +
    ggplot2::geom_line(data = lin_proj$band,
                       ggplot2::aes(x = !!exp_sym, y = hi),
                       color = "red", linewidth = 0.5, linetype = "dashed",
                       inherit.aes = FALSE) +
    x_scale_for(spec$exposure) +
    ggplot2::scale_y_continuous(labels = format_thousands) +
    ggplot2::labs(x = x_axis_label, y = y_label) +
    theme_pub() +
    ggplot2::theme(axis.title = ggplot2::element_text(face = "plain", size = 9))

  # Apply per-column x-limit (used by shared-axes figure assembly to align
  # cells in the same exposure column). coord_cartesian rather than
  # scale_x_continuous(limits=) because the latter clips data outside the
  # range, which would drop simul-band edges; coord_cartesian only changes
  # the view.
  if (!is.null(x_limits)) p <- p + ggplot2::coord_cartesian(xlim = x_limits)
  p
}

make_amef_panel <- function(spec, simul_bands_amef, lin_proj,
                            slope_draws = NULL, y_label = NULL,
                            x_axis_label = NULL,
                            x_limits = NULL,
                            show_x_axis = TRUE) {
  exp_sym <- rlang::sym(spec$exposure)
  bands   <- simul_bands_amef |> dplyr::arrange(dplyr::desc(level))
  if (is.null(y_label))      y_label      <- paste0("Slope of predicted\n", spec$outcome_label)
  if (is.null(x_axis_label)) x_axis_label <- spec$exposure_label

  # Linear-projection slope band: bound to the simul-band grid (1st-99th
  # percentile of observed exposure), matching the AERF panel's envelope
  # extent. Using annotate("segment") / annotate("rect") with explicit grid
  # endpoints rather than -Inf / +Inf so the band doesn't extend beyond the
  # data when coord_cartesian(xlim = ...) widens the plot view (shared-axes
  # layout).
  grid_lo <- min(bands[[spec$exposure]], na.rm = TRUE)
  grid_hi <- max(bands[[spec$exposure]], na.rm = TRUE)

  p <- ggplot2::ggplot(bands, ggplot2::aes(x = !!exp_sym)) +
    nested_simul_ribbons(bands) +
    # red-dashed linear-projection slope HPDI: rect + 2 horizontals, all
    # bounded to the simul-band grid range.
    ggplot2::annotate("rect",
                      xmin = grid_lo, xmax = grid_hi,
                      ymin = lin_proj$beta_hpdi$lo, ymax = lin_proj$beta_hpdi$hi,
                      fill = "red", alpha = 0.10) +
    ggplot2::annotate("segment",
                      x = grid_lo, xend = grid_hi,
                      y = lin_proj$beta_hpdi$lo, yend = lin_proj$beta_hpdi$lo,
                      color = "red", linewidth = 0.5, linetype = "dashed") +
    ggplot2::annotate("segment",
                      x = grid_lo, xend = grid_hi,
                      y = lin_proj$beta_hpdi$hi, yend = lin_proj$beta_hpdi$hi,
                      color = "red", linewidth = 0.5, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
    x_scale_for(spec$exposure) +
    ggplot2::scale_y_continuous(labels = format_thousands) +
    ggplot2::labs(x = x_axis_label, y = y_label) +
    theme_pub() +
    ggplot2::theme(axis.title = ggplot2::element_text(face = "plain", size = 9))

  # Rug at the bottom of the AMEF panel from observed exposure values
  # (trimmed to the 1st-99th percentile to match the prediction grid).
  if (!is.null(slope_draws)) {
    grid_lo <- min(slope_draws[[spec$exposure]], na.rm = TRUE)
    grid_hi <- max(slope_draws[[spec$exposure]], na.rm = TRUE)
    rug_df <- get(spec$data) |>
      dplyr::filter(!is.na(.data[[spec$exposure]])) |>
      dplyr::filter(.data[[spec$exposure]] >= grid_lo,
                    .data[[spec$exposure]] <= grid_hi)
    p <- p + ggplot2::geom_rug(
      data = rug_df,
      ggplot2::aes(x = !!exp_sym),
      sides = "b", alpha = 0.3, length = ggplot2::unit(0.02, "npc"),
      inherit.aes = FALSE
    )
  }

  # Hide the x-axis text / title / ticks when this AMEF panel is not on the
  # bottom row of its column (shared-axes assembly).
  if (!show_x_axis) {
    p <- p + ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  }

  # Apply per-column x-limit (shared-axes assembly). coord_cartesian rather
  # than scale_x_continuous(limits=) so simul-band edges don't get clipped.
  if (!is.null(x_limits)) p <- p + ggplot2::coord_cartesian(xlim = x_limits)

  p
}

# Age-conditional AMEF panel. Overlays median AMEF curves at multiple age
# slices (color-coded by age) for a single spec. Used by the age-conditional AMEF
# supplement figure that addresses Ian's 2026-05-25 question about whether
# PA effects vary by age. Uncertainty is omitted at the panel level (4
# overlapping ribbons would be visually unreadable); readers should refer
# to the population-average AMEF in the corresponding Figure 4 panel for
# the credible-interval extent.

make_age_cond_amef_panel <- function(spec, amef_age_draws,
                                     y_label = NULL, x_axis_label = NULL,
                                     x_limits = NULL, show_x_axis = TRUE,
                                     show_legend = FALSE) {
  exp_sym <- rlang::sym(spec$exposure)

  summary_df <- amef_age_draws |>
    dplyr::group_by(age_years, !!exp_sym) |>
    dplyr::summarize(med = stats::median(draw), .groups = "drop") |>
    dplyr::mutate(age_years = factor(age_years))

  if (is.null(y_label))      y_label      <- paste0("Slope of predicted\n", spec$outcome_label)
  if (is.null(x_axis_label)) x_axis_label <- spec$exposure_label

  p <- ggplot2::ggplot(summary_df,
                       ggplot2::aes(x = !!exp_sym, y = med, color = age_years)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
    ggplot2::geom_line(linewidth = 0.8, key_glyph = "rect") +
    ggplot2::scale_color_manual(
      name   = "Age (years)",
      values = c("25" = "#0072B2",   # Okabe-Ito blue
                 "35" = "#009E73",   # bluish green
                 "50" = "#E69F00",   # orange
                 "65" = "#CC79A7")   # reddish purple
    ) +
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


# Age-conditional AMEF panel WITH 95% HPDI bands. Same overall structure
# as make_age_cond_amef_panel but additionally renders posterior uncertainty
# ribbons for each age slice. Useful for judging when apparent age-conditional
# divergences are tensor-smooth artifacts (large HPDIs in sparse data corners)
# vs. signal (tight HPDIs that nevertheless differ between ages). Visually
# busier than the medians-only version; both are kept side by side as
# supp-fig-7-age-conditional-bands.pdf (the medians-only companion was dropped).

make_age_cond_amef_panel_with_bands <- function(spec, amef_age_draws,
                                                y_label = NULL, x_axis_label = NULL,
                                                x_limits = NULL, show_x_axis = TRUE,
                                                show_legend = FALSE) {
  exp_sym <- rlang::sym(spec$exposure)

  # Compute 95% SIMULTANEOUS credible bands per age slice (not pointwise
  # HPDIs). Matches the convention used by main-text Fig 3 / Fig 4 AMEF
  # panels: sup-norm calibrated for joint coverage of the entire curve,
  # not pointwise per-x coverage. The bands are computed per age slice
  # because each slice has its own posterior of AMEF curves; pooling
  # across ages would average over the conditioning the figure is meant
  # to expose.
  summary_df <- amef_age_draws |>
    dplyr::group_split(age_years) |>
    purrr::map_dfr(function(df) {
      bands <- simul_credible_bands(
        df, exposure = !!exp_sym,
        value_col = "draw", levels = 0.95,
        function_type = "AMEF", interval_type = "HPDI"
      )
      bands$age_years <- df$age_years[1]
      bands
    }) |>
    dplyr::mutate(age_years = factor(age_years))

  if (is.null(y_label))      y_label      <- paste0("Slope of predicted\n", spec$outcome_label)
  if (is.null(x_axis_label)) x_axis_label <- spec$exposure_label

  age_colors <- c("25" = "#0072B2",   # Okabe-Ito blue
                  "35" = "#009E73",   # bluish green
                  "50" = "#E69F00",   # orange
                  "65" = "#CC79A7")   # reddish purple

  p <- ggplot2::ggplot(summary_df,
                       ggplot2::aes(x = !!exp_sym,
                                    ymin = lo, ymax = hi,
                                    fill = age_years)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
    ggplot2::geom_ribbon(alpha = 0.22, color = NA) +
    ggplot2::scale_fill_manual(name  = "Age (years)", values = age_colors) +
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


# Stack the AERF panel over the AMEF panel into a single subfigure. The
# AERF panel loses its x-axis (carried by the AMEF below); both panels are
# aligned vertically. Optional letter tag goes to the top-left of the AERF.

stack_subfig <- function(aerf, amef, tag = NULL) {
  aerf <- aerf + ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.text.x  = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  )
  if (!is.null(tag)) aerf <- aerf + ggplot2::labs(tag = tag)
  patchwork::wrap_plots(aerf, amef, ncol = 1, heights = c(1, 1)) &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "plain", size = 12))
}


# ---- aerf_amef_plot: legacy 3-panel (AERF / AMEF / curvature) ----
#
# Retained for backwards-compat with any caller that still expects the
# energetic-constraints-style 3-panel composite. New bone-turnover figures
# use the AERF + AMEF 2-panel stack via make_aerf_panel / make_amef_panel /
# stack_subfig above (graveyard fig-5 convention).

aerf_amef_plot <- function(spec, pred_draws, slope_draws,
                           curvature_draws = NULL, truth_df = NULL) {

  ex <- spec$exposure
  exp_sym <- rlang::sym(ex)

  if (is.null(curvature_draws)) {
    curvature_draws <- compute_curvature_draws(slope_draws, ex)
  }

  # X-axis ticks: rescaled-step columns (`ad_steps_1k`) get a "k" suffix.
  x_scale <- if (grepl("steps_1k", ex)) {
    ggplot2::scale_x_continuous(labels = function(x) paste0(x, "k"))
  } else {
    ggplot2::scale_x_continuous()
  }

  # Strip text combines panel role + y-axis quantity so each facet
  # carries both the name (AERF / AMEF / Curvature) and what's plotted
  # in one place. y-axis label is dropped (the strip carries it).
  panel_levels <- c(
    paste0("AERF: ",                       spec$outcome_label),
    paste0("AMEF: ",                       spec$amef_label),
    paste0("Curvature (slope of AMEF): ",  spec$curvature_label)
  )

  combined <- dplyr::bind_rows(
    pred_draws      |> dplyr::mutate(panel = panel_levels[1]),
    slope_draws     |> dplyr::mutate(panel = panel_levels[2]),
    curvature_draws |> dplyr::mutate(panel = panel_levels[3])
  ) |>
    dplyr::mutate(panel = factor(panel, levels = panel_levels))

  summary_df <- combined |>
    dplyr::group_by(panel, !!exp_sym) |>
    ggdist::median_hdci(draw, .width = c(0.5, 0.8, 0.95)) |>
    dplyr::ungroup()

  # Zero-reference line only on the AMEF and Curvature panels.
  zero_lines <- tibble::tibble(
    panel = factor(panel_levels[2:3], levels = panel_levels),
    yint  = 0
  )

  p <- ggplot2::ggplot(summary_df,
                       ggplot2::aes(x = !!exp_sym, y = draw,
                                    ymin = .lower, ymax = .upper)) +
    ggplot2::geom_hline(data = zero_lines,
                        ggplot2::aes(yintercept = yint),
                        linetype = "dashed", color = "grey50",
                        inherit.aes = FALSE) +
    ggdist::geom_lineribbon(alpha = 0.6) +
    ggplot2::scale_fill_brewer(palette = "Oranges") +
    x_scale +
    ggplot2::facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    ggplot2::labs(x = spec$exposure_label, y = NULL, fill = "HPDI") +
    theme_pub_leg() +
    ggplot2::theme(
      strip.text     = ggplot2::element_text(face = "bold", hjust = 0,
                                             size = ggplot2::rel(0.95)),
      panel.spacing  = ggplot2::unit(0.6, "lines")
    )

  if (!is.null(truth_df)) {
    truth_df <- truth_df |>
      dplyr::mutate(panel = factor(panel, levels = panel_levels))
    p <- p + ggplot2::geom_line(
      data = truth_df,
      ggplot2::aes(x = !!exp_sym, y = true),
      color = "black", linetype = "dashed", linewidth = 0.7,
      inherit.aes = FALSE
    )
  }
  p
}


# ---- Posterior predictive checks ----

# Six-panel PP check: min, mean, median, max, SD, skew.
pp_check_stats <- function(model, ndraws = 100) {
  p1 <- pp_check(model, type = "stat", stat = "min",    ndraws = ndraws) + ggtitle("Min")    + theme_pub()
  p2 <- pp_check(model, type = "stat", stat = "mean",   ndraws = ndraws) + ggtitle("Mean")   + theme_pub()
  p3 <- pp_check(model, type = "stat", stat = "median", ndraws = ndraws) + ggtitle("Median") + theme_pub()
  p4 <- pp_check(model, type = "stat", stat = "max",    ndraws = ndraws) + ggtitle("Max")    + theme_pub()
  p5 <- pp_check(model, type = "stat", stat = "sd",     ndraws = ndraws) + ggtitle("SD")     + theme_pub()
  p6 <- pp_check(model, type = "stat", stat = function(x) skew(x), ndraws = ndraws) + ggtitle("Skew") + theme_pub()
  (p1 + p2 + p3) / (p4 + p5 + p6)
}

# Same six panels, but stratified by a grouping factor (passed as a string).
pp_check_stats_grouped <- function(model, group, ndraws = 100) {

  yrep <- brms::posterior_predict(model, ndraws = ndraws)   # ndraws x n
  resp <- as.character(model$formula$formula[[2]])
  y    <- model$data[[resp]]
  g    <- model$data[[group]]
  if (is.null(y)) stop("Response column '", resp, "' not found in model$data")

  stats <- list(
    Min    = function(x) min(x,    na.rm = TRUE),
    Mean   = function(x) mean(x,   na.rm = TRUE),
    Median = function(x) median(x, na.rm = TRUE),
    Max    = function(x) max(x,    na.rm = TRUE),
    SD     = function(x) sd(x,     na.rm = TRUE),
    Skew   = function(x) skew(x)
  )

  group_levels <- sort(unique(as.character(g)))
  stat_levels  <- names(stats)

  rep_rows <- list()
  obs_rows <- list()
  for (s in stat_levels) {
    f <- stats[[s]]
    for (lv in group_levels) {
      idx <- which(as.character(g) == lv)
      if (length(idx) == 0) next
      v <- vapply(seq_len(nrow(yrep)),
                  function(i) f(yrep[i, idx]),
                  numeric(1))
      rep_rows[[length(rep_rows) + 1]] <- data.frame(
        stat = factor(s,  levels = stat_levels),
        grp  = factor(lv, levels = group_levels),
        val  = as.numeric(v)
      )
      obs_rows[[length(obs_rows) + 1]] <- data.frame(
        stat = factor(s,  levels = stat_levels),
        grp  = factor(lv, levels = group_levels),
        val  = f(y[idx])
      )
    }
  }
  yrep_df <- do.call(rbind, rep_rows)
  obs_df  <- do.call(rbind, obs_rows)

  ggplot(yrep_df, aes(x = val)) +
    geom_histogram(bins = 25, fill = "#A6BEE0", color = "#5C7BA0",
                   linewidth = 0.15) +
    geom_vline(data = obs_df, aes(xintercept = val),
               color = "#1F2937", linewidth = 0.7) +
    facet_grid(stat ~ grp, scales = "free", switch = "y") +
    labs(x = NULL, y = NULL,
         title = paste0("Posterior-predictive stats, faceted by ", group),
         subtitle = "histogram = yrep across draws; dark line = observed") +
    theme_pub() +
    theme(
      strip.placement      = "outside",
      strip.background     = element_rect(fill = "grey95", color = NA),
      strip.text.y.left    = element_text(angle = 0, face = "bold", size = 9),
      strip.text.x.top     = element_text(face = "bold", size = 9),
      axis.text.x          = element_text(size = 7),
      axis.text.y          = element_blank(),
      axis.ticks.y         = element_blank(),
      panel.grid.minor     = element_blank(),
      panel.spacing.x      = unit(0.4, "lines"),
      panel.spacing.y      = unit(0.3, "lines")
    )
}


# ---- approx_flat: posterior probability the AMEF is practically flat ----
#
# For each posterior draw of the AMEF (slope curve), record the largest
# |slope| anywhere on the exposure grid. p_flat is the fraction of draws
# whose maximum |slope| stays below `eps_flat`. Returns one row per eps
# threshold so the user can sweep across "what counts as negligible."

approx_flat <- function(
    slope_draws,
    exposure,
    eps_flat,
    per_x_units = 1,
    value_col = "draw"
) {

  stopifnot(is.data.frame(slope_draws))

  exposure_sym <- rlang::ensym(exposure)
  exposure_nm  <- rlang::as_name(exposure_sym)

  req_cols <- c("drawid", exposure_nm, value_col)
  if (!all(req_cols %in% names(slope_draws))) {
    stop("slope_draws must contain drawid, exposure, and value columns")
  }

  if (!is.numeric(per_x_units) || length(per_x_units) != 1 || !is.finite(per_x_units) || per_x_units <= 0) {
    stop("per_x_units must be a single positive number (e.g., 1, 10, 100).")
  }

  if (!is.numeric(eps_flat) || length(eps_flat) < 1 || any(!is.finite(eps_flat)) || any(eps_flat < 0)) {
    stop("eps_flat must be a numeric scalar or vector of finite, non-negative values.")
  }

  by_draw_max_abs <- slope_draws |>
    dplyr::filter(!is.na(.data[[value_col]])) |>
    dplyr::group_by(drawid) |>
    dplyr::summarize(
      max_abs = max(abs(.data[[value_col]]), na.rm = TRUE),
      .groups = "drop"
    )

  by_draw_max_change <- by_draw_max_abs |>
    dplyr::mutate(max_change = max_abs * per_x_units)

  curve <- tibble::tibble(eps_flat = eps_flat) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      per_x_units       = per_x_units,
      p_approx_flat     = mean(by_draw_max_change$max_change < eps_flat, na.rm = TRUE),
      p_approx_not_flat = 1 - p_approx_flat
    ) |>
    dplyr::ungroup()

  list(
    curve              = curve,
    by_draw_max_abs    = by_draw_max_abs,
    by_draw_max_change = by_draw_max_change,
    per_x_units        = per_x_units
  )
}


# ---- approx_linearity: posterior probability the AERF is approximately linear ----
#
# For each posterior draw of the AERF, fit the best straight line on the
# exposure grid and record the maximum residual (largest deviation from
# linearity). p_approx_linear is the fraction of draws whose max residual
# stays below tol_percent% of an outcome scale (sd / mad / range).

approx_linearity <- function(
    pred_draws,
    model,
    exposure,
    resp,
    outcome,
    tol_percent,
    metric = c("sd", "mad", "range"),
    outcome_multiplier = 1,
    value_col = "draw"
) {

  stopifnot(is.data.frame(pred_draws))

  exposure_sym <- rlang::ensym(exposure)
  exposure_nm  <- rlang::as_name(exposure_sym)

  req_cols <- c("drawid", exposure_nm, value_col)
  if (!all(req_cols %in% names(pred_draws))) {
    stop("pred_draws must contain drawid, exposure, and value columns")
  }

  max_dev_draws <- pred_draws |>
    dplyr::filter(!is.na(.data[[value_col]])) |>
    dplyr::group_by(drawid) |>
    dplyr::summarize(
      # Max deviation of the AERF from its closed-form OLS projection (per draw):
      # an lm-free posterior functional, identical to abs(resid(lm(...))).
      max_dev = {
        x <- .data[[exposure_nm]]; y <- .data[[value_col]]
        xc <- x - mean(x); b <- sum(xc * (y - mean(y))) / sum(xc^2)
        a <- mean(y) - b * mean(x); max(abs(y - (a + b * x)), na.rm = TRUE)
      },
      .groups = "drop"
    )

  mf <- model.frame(model, resp = resp)

  if (!(outcome %in% names(mf))) {
    stop("`outcome` must be a column in model.frame(model, resp = resp).")
  }

  y_used_fit <- mf[[outcome]]

  metric  <- unique(tolower(metric))
  allowed <- c("sd", "mad", "range")
  if (!all(metric %in% allowed)) {
    stop("metric must be one or more of: 'sd', 'mad', 'range'")
  }

  scale_tbl <- tibble::tibble(metric = metric) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      scale_fit = dplyr::case_when(
        metric == "sd"    ~ stats::sd(y_used_fit, na.rm = TRUE),
        metric == "mad"   ~ stats::mad(y_used_fit, na.rm = TRUE),
        metric == "range" ~ max(y_used_fit, na.rm = TRUE) - min(y_used_fit, na.rm = TRUE)
      ),
      scale_pred = scale_fit * outcome_multiplier
    ) |>
    dplyr::ungroup()

  if (any(!is.finite(scale_tbl$scale_fit) | scale_tbl$scale_fit <= 0)) {
    stop("computed outcome scale is not finite/positive.")
  }

  curve <- tidyr::expand_grid(
    scale_tbl,
    tol_percent = tol_percent
  ) |>
    dplyr::mutate(
      tol_abs = (tol_percent / 100) * scale_pred
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      p_approx_linear     = mean(max_dev_draws$max_dev < tol_abs, na.rm = TRUE),
      p_approx_non_linear = 1 - p_approx_linear
    ) |>
    dplyr::ungroup() |>
    dplyr::select(
      metric,
      tol_percent,
      tol_abs,
      p_approx_linear,
      p_approx_non_linear
    )

  list(
    curve              = curve,
    scale              = scale_tbl,
    max_dev_draws      = max_dev_draws,
    outcome_multiplier = outcome_multiplier
  )
}


# ---- simul_credible_bands: simultaneous credible bands (vector of levels) ----
#
# Standardized sup-norm bands with empirical-quantile calibration of c*.
# The band width at each x is c* * pointwise_sd(x); c* is chosen so that the
# fraction of draws whose entire curve lies inside the band is ~ level.
# Returns a long-format tibble with one row per (exposure x level) for
# plotting nested bands (e.g. 50/80/95%) on a single panel.

simul_credible_bands <- function(
    draws_df,
    exposure,
    value_col     = "draw",
    levels        = c(0.95),
    function_type = c("AERF", "AMEF"),
    interval_type = c("perc", "HPDI")
) {

  function_type <- match.arg(function_type)
  interval_type <- match.arg(interval_type)

  stopifnot(is.data.frame(draws_df))
  if (!("drawid" %in% names(draws_df))) {
    stop("draws_df must contain a column named 'drawid'.")
  }
  if (!(value_col %in% names(draws_df))) {
    stop(sprintf("draws_df must contain the value column '%s'.", value_col))
  }

  if (!is.numeric(levels) || any(!is.finite(levels)) || any(levels <= 0) || any(levels >= 1)) {
    stop("levels must be numeric values strictly between 0 and 1 (e.g., c(0.5, 0.8, 0.95)).")
  }

  exposure_sym <- rlang::ensym(exposure)
  exposure_nm  <- rlang::as_name(exposure_sym)
  if (!(exposure_nm %in% names(draws_df))) {
    stop(sprintf("draws_df must contain the exposure column '%s'.", exposure_nm))
  }

  summ <- draws_df |>
    dplyr::group_by(!!exposure_sym) |>
    dplyr::summarize(
      m  = mean(.data[[value_col]], na.rm = TRUE),
      sd = stats::sd(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::distinct(.data[[exposure_nm]], .keep_all = TRUE)

  if (any(!is.finite(summ$sd)) || any(summ$sd <= 0, na.rm = TRUE)) {
    stop("Some exposure grid points have sd <= 0 across draws; cannot form standardized sup-norm band.")
  }

  med_tbl <- draws_df |>
    dplyr::group_by(!!exposure_sym) |>
    dplyr::summarize(
      med = stats::median(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::distinct(.data[[exposure_nm]], .keep_all = TRUE)

  summ <- summ |>
    dplyr::left_join(med_tbl, by = exposure_nm)

  t_by_draw <- draws_df |>
    dplyr::left_join(summ, by = exposure_nm) |>
    dplyr::mutate(z = (.data[[value_col]] - .data[["m"]]) / .data[["sd"]]) |>
    dplyr::group_by(.data[["drawid"]]) |>
    dplyr::summarize(t = max(abs(.data[["z"]]), na.rm = TRUE), .groups = "drop")

  t_sorted <- sort(t_by_draw$t)
  n_t <- length(t_sorted)
  if (n_t == 0) stop("No finite t values found; cannot construct bands.")

  levels <- sort(unique(levels))

  c_tbl <- tibble::tibble(level = levels) |>
    dplyr::mutate(
      idx    = pmax(1L, pmin(n_t, ceiling(level * n_t))),
      c_star = t_sorted[idx]
    ) |>
    dplyr::select(level, c_star)

  bands <- summ |>
    tidyr::crossing(c_tbl) |>
    dplyr::mutate(
      lo = .data[["m"]] - .data[["c_star"]] * .data[["sd"]],
      hi = .data[["m"]] + .data[["c_star"]] * .data[["sd"]]
    ) |>
    dplyr::mutate(
      function_type = function_type,
      interval_type = interval_type
    ) |>
    dplyr::select(
      !!exposure_sym, m, sd, med, lo, hi, level, function_type, interval_type
    )

  bands
}
