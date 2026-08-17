###############################################################################
# Walking-experiment reverse-ROPE analysis.
#
# Asks whether an acute bout of walking produces a short-term rise in
# osteocalcin, using the reverse-ROPE approach applied throughout this project:
# fit a Bayesian model, then for each post-exercise timepoint report U95 -- the
# 95% posterior UPPER bound on the magnitude of the osteocalcin change --
# interpreted against the SD-fraction (0.2 SD landmark) and the assay
# resolvability floor. A small U95 is positive evidence for a practical null; a
# large U95 means the non-detection is uninformative (power artifact).
#
# The question is framed this way, rather than as a significance test, because a
# non-significant result at this sample size cannot separate a genuine null from
# a power artifact. U95 distinguishes the two directly.
#
# Design: within-subject, n = 10, osteocalcin (ng/mL) at pre / t0 (immediately
# post) / t4 (4 h post). Identification is by design (the walking protocol IS
# the intervention) -- no DAG / adjustment set; age & sex are constant within
# person and absorbed by the participant term.
#
# Reported study inference is Bayesian (brms). The a priori resolution check
# (does n = 10 even have the resolution to support a practical-null claim?)
# uses a frequentist proxy, which this project uses for calibration and power
# only, never for reported study inference.
#
# Outputs:
#   outputs/models/walking-experiment/model.Rdata          (gitignored)
#   outputs/models/walking-experiment/reverse-rope.csv     (contrast summaries)
#   outputs/models/walking-experiment/calibration.csv      (resolution check)
#   outputs/figures/working/walking-experiment/walking-experiment-rope.pdf
###############################################################################


source(here::here("code", "_startup", "init.R"))


# ---- Config ----

set.seed(SEED)

# Assay CV proxy: main-study R-PLEX Human Osteocalcin assay, same laboratory.
# Used only to place the resolvability floor; a CV specific to this experiment
# was not separately available. The resolvability floor follows the same
# convention as the main analyses: total analytical CV (intra ⊕ inter, in
# quadrature) at the sample median.
INTRA_CV <- 0.020
INTER_CV <- 0.028
TOTAL_CV <- sqrt(INTRA_CV^2 + INTER_CV^2)   # 3.44%, matching the main analyses

SESOI_SD_FRAC <- 0.2   # "small effect" landmark, matching the main paper
NSIM          <- 2000  # calibration sims per scenario (matches power-sim convention)

out_dir_models <- here::here("outputs", "models", "walking-experiment")
out_dir_fig    <- here::here("outputs", "figures", "working", "walking-experiment")
dir.create(out_dir_models, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_fig,    recursive = TRUE, showWarnings = FALSE)


# ---- Data ----

dat_exp <- readr::read_csv(
  here::here("data", "processed", "walking-experiment-data.csv"),
  show_col_types = FALSE
) |>
  janitor::clean_names(case = "snake") |>
  dplyr::mutate(
    time_point     = factor(time_point, levels = c("pre", "t0", "t4")),
    participant_id = factor(participant_id)
  )

# Baseline (pre) reference scales.
pre_vals    <- dat_exp |> dplyr::filter(time_point == "pre") |> dplyr::pull(osteocalcin_ng_ml)
mean_pre    <- mean(pre_vals)
median_pre  <- median(pre_vals)
sd_pre      <- stats::sd(pre_vals)                 # between-person SD (ng/mL)
sesoi       <- SESOI_SD_FRAC * sd_pre       # 0.2 SD landmark (ng/mL)
noise_floor <- median_pre * TOTAL_CV        # assay resolvability floor (ng/mL)

cat(sprintf(
  "\nBaseline osteocalcin: mean = %.1f, median = %.1f, between-person SD = %.1f ng/mL\n  0.2 SD landmark = %.1f ng/mL | assay noise floor (total CV %.2f%% at median) = %.1f ng/mL\n",
  mean_pre, median_pre, sd_pre, sesoi, 100 * TOTAL_CV, noise_floor
))


# ---- Fit: hierarchical lognormal ----

fit_path <- file.path(out_dir_models, "model.Rdata")

# A bare file.exists() check reuses whatever fit happens to be on disk, which is
# wrong the moment the data or the sampling configuration change -- and the
# synthetic example generator writes to the same input filename the real data
# uses, so "fit on synthetic, then obtain the real data" is a path a reader can
# actually take. A brmsfit carries the data it was fitted with, so the cached
# object is checked against the current inputs rather than trusted.
expected_draws <- CHAINS * (ITER - WARMUP) / THIN
cache_mismatch <- function(fit) {
  if (!inherits(fit, "brmsfit")) return("not a brmsfit")
  cols <- names(fit$data)
  if (!all(cols %in% names(dat_exp))) return("stored data has unexpected columns")
  if (!isTRUE(all.equal(fit$data, as.data.frame(dat_exp)[, cols, drop = FALSE],
                        check.attributes = FALSE)))
    return("stored data differs from the current walking-experiment data")
  if (brms::ndraws(fit) != expected_draws)
    return(sprintf("stored fit has %d draws, expected %d", brms::ndraws(fit), expected_draws))
  NA_character_
}

cached_ok <- FALSE
if (file.exists(fit_path)) {
  load(fit_path)
  why <- cache_mismatch(model_fit)
  cached_ok <- is.na(why)
  if (!cached_ok) cat(sprintf("Cached fit rejected (%s); refitting.\n", why))
}

if (cached_ok) {
  cat("Loading cached model...\n")
} else {
  cat("Fitting hierarchical lognormal (osteocalcin ~ time_point + (1|participant))...\n")
  model_fit <- brms::brm(
    brms::bf(osteocalcin_ng_ml ~ time_point + (1 | participant_id),
       family = brms::lognormal(link = "identity", link_sigma = "log")),
    data    = dat_exp,
    prior   = brms::set_prior("student_t(3, 0, 2.5)", class = "b"),
    warmup  = WARMUP,
    iter    = ITER,
    thin    = THIN,
    chains  = CHAINS,
    cores   = CHAINS,
    seed    = SEED,
    backend = "cmdstanr",
    control = BRMS_CONTROL,
    refresh = 0,
    silent  = 2
  )
  save(model_fit, file = fit_path, compress = "gzip")
}

# Convergence sanity check.
rhat_max <- max(brms::rhat(model_fit), na.rm = TRUE)
# brms::nuts_params rather than rstan::get_sampler_params: backend-agnostic, so
# it works under cmdstanr without attaching rstan. The two read the same table.
np       <- brms::nuts_params(model_fit)
ndiv     <- sum(np$Value[np$Parameter == "divergent__"])
cat(sprintf("Convergence: max Rhat = %.3f | divergent transitions = %d\n", rhat_max, ndiv))


# ---- Contrasts vs baseline (population-average, response scale) ----

# avg_comparisons averages each non-reference timepoint vs the reference (pre)
# over the observed sample -> population-average change in ng/mL (g-computation).
cmp <- marginaleffects::avg_comparisons(
  model_fit,
  variables = list(time_point = "reference"),
  type      = "response"
)
cmp_draws <- marginaleffects::posterior_draws(cmp)

contrast_levels <- unique(cmp_draws$contrast)

draws_by_contrast <- lapply(contrast_levels, function(lv) {
  cmp_draws$draw[cmp_draws$contrast == lv]
})
names(draws_by_contrast) <- contrast_levels


# ---- Reverse-ROPE summaries (U95 = 95% upper bound on |change|) ----

summarize_contrast <- function(d, label) {
  # ggdist::hdci, the project's one scalar interval primitive. Not
  # interchangeable with HDInterval::hdi, which spans floor(n * q) + 1 order
  # statistics and so carries slightly more than q, where hdci interpolates to
  # carry exactly q. The reported U95 below is a plain quantile either way.
  hdi <- ggdist::hdci(d, .width = 0.95)
  u95 <- as.numeric(quantile(abs(d), 0.95))
  tibble::tibble(
    contrast        = label,
    p_increase      = mean(d > 0),
    hpdi_lo_ngml    = unname(hdi[1]),
    hpdi_hi_ngml    = unname(hdi[2]),
    u95_ngml        = u95,
    u95_frac_sd     = u95 / sd_pre,
    u95_pct_base    = 100 * u95 / mean_pre,
    below_0p2_sd    = u95 < sesoi,
    above_noise     = u95 > noise_floor
  )
}

rope_tbl <- purrr::map2_dfr(draws_by_contrast, names(draws_by_contrast), summarize_contrast)
readr::write_csv(rope_tbl, file.path(out_dir_models, "reverse-rope.csv"))

cat("\n==== Reverse-ROPE summary (Bayesian) ====\n")
print(as.data.frame(rope_tbl), digits = 3)


# ---- A priori resolution check (frequentist proxy; calibration only) ----

# Within-person change noise (ng/mL): SD of observed post-minus-pre differences,
# pooled over t0 and t4. Treated as additive noise on the natural scale.
wide <- dat_exp |>
  dplyr::select(participant_id, time_point, osteocalcin_ng_ml) |>
  tidyr::pivot_wider(names_from = time_point, values_from = osteocalcin_ng_ml)
diffs_pooled <- c(wide$t0 - wide$pre, wide$t4 - wide$pre)
sigma_w      <- stats::sd(diffs_pooled)
n_part       <- length(pre_vals)

cat(sprintf("\nWithin-person change SD (pooled t0/t4) = %.1f ng/mL | n = %d\n", sigma_w, n_part))

# For one simulated dataset under a true mean change `true_mu`, apply the paired
# frequentist analysis and the reverse-ROPE decision rules.
simulate_scenario <- function(true_mu, nsim = NSIM) {
  detect    <- logical(nsim)
  prac_null <- logical(nsim)
  u95v      <- numeric(nsim)
  tcrit95   <- qt(0.975, df = n_part - 1)
  tcrit90   <- qt(0.95,  df = n_part - 1)
  for (i in seq_len(nsim)) {
    delta <- true_mu + rnorm(n_part, 0, sigma_w)   # per-person change
    m     <- mean(delta)
    se    <- stats::sd(delta) / sqrt(n_part)
    ci_lo <- m - tcrit95 * se
    ci_hi <- m + tcrit95 * se
    detect[i]    <- (ci_lo > 0) || (ci_hi < 0)     # 95% CI excludes 0
    u95          <- abs(m) + tcrit90 * se          # one-sided 95% upper bound on |effect|
    u95v[i]      <- u95
    prac_null[i] <- u95 < sesoi                    # rules out an effect as large as 0.2 SD
  }
  tibble::tibble(
    scenario         = NA_character_,
    true_mu_ngml     = true_mu,
    mean_u95_ngml    = mean(u95v),
    p_detect         = mean(detect),
    p_practical_null = mean(prac_null)
  )
}

set.seed(SEED)   # deterministic calibration regardless of upstream RNG / model-cache state
calib <- dplyr::bind_rows(
  simulate_scenario(0)     |> dplyr::mutate(scenario = "true null"),
  simulate_scenario(sesoi) |> dplyr::mutate(scenario = "true 0.2 SD effect")
)
readr::write_csv(calib, file.path(out_dir_models, "calibration.csv"))

cat("\n==== A priori resolution check (frequentist proxy) ====\n")
print(as.data.frame(calib), digits = 3)
cat(sprintf(
  "\nInterpretation: under a true null the design concludes 'practical null' %.0f%% of the time\n  (false-detect rate %.1f%%); against a true 0.2 SD effect power is %.0f%%\n  (wrongly-null %.0f%%).\n",
  100 * calib$p_practical_null[1], 100 * calib$p_detect[1],
  100 * calib$p_detect[2], 100 * calib$p_practical_null[2]
))


# ---- Figures ----

okabe <- c("t0 - pre" = "#0072B2", "t4 - pre" = "#D55E00")
lab_for <- function(lv) c("t0 - pre" = "t0 (immediately post)",
                          "t4 - pre" = "t4 (4 h post)")[lv]

# Panel A: per-participant trajectories.
pA <- ggplot2::ggplot(dat_exp,
             ggplot2::aes(x = time_point, y = osteocalcin_ng_ml, group = participant_id)) +
  ggplot2::geom_line(alpha = 0.45, linewidth = 0.4, color = "grey40") +
  ggplot2::geom_point(alpha = 0.7, size = 1.4, color = "grey20") +
  ggplot2::scale_x_discrete(labels = c(pre = "pre", t0 = "t0", t4 = "t4")) +
  ggplot2::labs(x = NULL, y = "Osteocalcin (ng/mL)", tag = "A",
       title = "Per-participant trajectories") +
  theme_pub() +
  ggplot2::theme(plot.tag = ggplot2::element_text(face = "plain", size = 12))

# Panel B: reverse-ROPE effect-size-vs-probability curves with U95 markers.
x_hi  <- max(vapply(draws_by_contrast, function(d) quantile(abs(d), 0.999), numeric(1)))
xgrid <- seq(0, x_hi, length.out = 250)
curve_df <- purrr::map2_dfr(draws_by_contrast, names(draws_by_contrast), function(d, lv) {
  tibble::tibble(contrast = lv, x = xgrid,
                 p = vapply(xgrid, function(t) mean(abs(d) > t), numeric(1)))
})
u95_df <- rope_tbl |> dplyr::transmute(contrast, u95_ngml)

pB <- ggplot2::ggplot(curve_df, ggplot2::aes(x = x, y = p, color = contrast)) +
  # 0.05, not 0.95. The y axis is the SURVIVAL probability P(|delta| > x), and
  # U95 is the 95th percentile of |delta| -- so each curve passes through 0.05,
  # not 0.95, at its own U95 vertical.
  ggplot2::geom_hline(yintercept = 0.05, color = "grey60", linetype = "dotted", linewidth = 0.3) +
  ggplot2::geom_vline(xintercept = noise_floor, color = "#56B4E9", linetype = "dashed", linewidth = 0.4) +
  ggplot2::geom_vline(xintercept = sesoi,       color = "grey45",  linetype = "dashed", linewidth = 0.4) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_vline(data = u95_df, ggplot2::aes(xintercept = u95_ngml, color = contrast),
             linetype = "solid", linewidth = 0.5, show.legend = FALSE) +
  ggplot2::annotate("text", x = noise_floor, y = 0.05, label = "assay\nnoise", hjust = -0.05,
           size = 2.5, color = "#56B4E9") +
  ggplot2::annotate("text", x = sesoi, y = 0.55, label = "0.2 SD", hjust = -0.1,
           size = 2.5, color = "grey45") +
  ggplot2::scale_color_manual(values = okabe, labels = lab_for, name = NULL) +
  ggplot2::labs(x = "Candidate effect size  |Δ osteocalcin|  (ng/mL)",
       y = "P( |Δ| > x )", tag = "B",
       title = "Reverse-ROPE: U95 = 95% upper bound on |change|") +
  theme_pub() +
  ggplot2::theme(legend.position = "top",
        plot.tag = ggplot2::element_text(face = "plain", size = 12))

fig <- patchwork::wrap_plots(pA, pB, ncol = 2, widths = c(1, 1.4))
ggplot2::ggsave(file.path(out_dir_fig, "walking-experiment-rope.pdf"),
       fig, width = 9, height = 4.0, device = cairo_pdf)

cat("\nSaved figure: outputs/figures/working/walking-experiment/walking-experiment-rope.pdf\n")
cat("Done: code/_experiments/walking-experiment-rope.R\n")
