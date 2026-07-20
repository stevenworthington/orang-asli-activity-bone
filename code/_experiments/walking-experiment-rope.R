###############################################################################
# Walking-experiment reverse-ROPE analysis (Supplementary Material 11).
#
# Re-analyzes the coauthors' acute walking experiment, which their original
# write-up summarized with Wilcoxon signed-rank tests (P >= 0.49). A
# non-significant Wilcoxon cannot separate a genuine null from a power
# artifact. This script answers that question with the project's reverse-ROPE
# machinery: fit a Bayesian model, then for each post-exercise timepoint report
# U95 -- the 95% posterior UPPER bound on the magnitude of the osteocalcin
# change -- interpreted against the SD-fraction (0.2 SD landmark) and the assay
# resolvability floor. A small U95 is positive evidence for a practical null; a
# large U95 means the non-detection is uninformative (power artifact).
#
# Design: within-subject, n = 10, osteocalcin (ng/mL) at pre / t0 (immediately
# post) / t4 (4 h post). Identification is by design (the walking protocol IS
# the intervention) -- no DAG / adjustment set; age & sex are constant within
# person and absorbed by the participant term.
#
# Reported study inference is Bayesian (brms). The a priori resolution check
# (does n = 10 even have the resolution to support a practical-null claim?)
# uses the project's sanctioned frequentist proxy for calibration / power only.
#
# Outputs:
#   outputs/models/walking-experiment/model.Rdata          (gitignored)
#   outputs/models/walking-experiment/reverse-rope.csv     (contrast summaries)
#   outputs/models/walking-experiment/calibration.csv      (resolution check)
#   outputs/figures/working/walking-experiment/walking-experiment-rope.pdf
###############################################################################


library(here)
source(here("code", "_startup", "init.R"))


# ---- Config ----

set.seed(SEED)

# Assay CV proxy: main-study R-PLEX Human Osteocalcin assay (same lab, Duke).
# Used only to place the resolvability floor; confirm the experiment's own CV
# with the coauthors when available. The resolvability floor follows the SM6
# convention: total analytical CV (intra ⊕ inter, in quadrature) at the sample
# median.
INTRA_CV <- 0.020
INTER_CV <- 0.028
TOTAL_CV <- sqrt(INTRA_CV^2 + INTER_CV^2)   # 3.44%, matching Supplementary Material 6

SESOI_SD_FRAC <- 0.2   # "small effect" landmark, matching the main paper
NSIM          <- 2000  # calibration sims per scenario (matches power-sim convention)

out_dir_models <- here("outputs", "models", "walking-experiment")
out_dir_fig    <- here("outputs", "figures", "working", "walking-experiment")
dir.create(out_dir_models, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_fig,    recursive = TRUE, showWarnings = FALSE)


# ---- Data ----

dat_exp <- readr::read_csv(
  here("data", "processed", "compiled-walking-experiment-18june2026.csv"),
  show_col_types = FALSE
) |>
  janitor::clean_names(case = "snake") |>
  mutate(
    time_point     = factor(time_point, levels = c("pre", "t0", "t4")),
    participant_id = factor(participant_id),
    sex            = factor(sex)
  )

# Baseline (pre) reference scales.
pre_vals    <- dat_exp |> dplyr::filter(time_point == "pre") |> dplyr::pull(osteocalcin_ng_ml)
mean_pre    <- mean(pre_vals)
median_pre  <- median(pre_vals)
sd_pre      <- sd(pre_vals)                 # between-person SD (ng/mL)
sesoi       <- SESOI_SD_FRAC * sd_pre       # 0.2 SD landmark (ng/mL)
noise_floor <- median_pre * TOTAL_CV        # assay resolvability floor (ng/mL), SM6 convention

cat(sprintf(
  "\nBaseline osteocalcin: mean = %.1f, median = %.1f, between-person SD = %.1f ng/mL\n  0.2 SD landmark = %.1f ng/mL | assay noise floor (total CV %.2f%% at median) = %.1f ng/mL\n",
  mean_pre, median_pre, sd_pre, sesoi, 100 * TOTAL_CV, noise_floor
))


# ---- Fit: hierarchical lognormal ----

fit_path <- file.path(out_dir_models, "model.Rdata")
if (file.exists(fit_path)) {
  cat("Loading cached model...\n")
  load(fit_path)
} else {
  cat("Fitting hierarchical lognormal (osteocalcin ~ time_point + (1|participant))...\n")
  model_fit <- brm(
    bf(osteocalcin_ng_ml ~ time_point + (1 | participant_id),
       family = lognormal(link = "identity", link_sigma = "log")),
    data    = dat_exp,
    prior   = set_prior("student_t(3, 0, 2.5)", class = "b"),
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
ndiv     <- sum(sapply(rstan::get_sampler_params(model_fit$fit, inc_warmup = FALSE),
                       function(x) sum(x[, "divergent__"])))
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
  hdi <- HDInterval::hdi(d, credMass = 0.95)
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
sigma_w      <- sd(diffs_pooled)
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
    se    <- sd(delta) / sqrt(n_part)
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
pA <- ggplot(dat_exp,
             aes(x = time_point, y = osteocalcin_ng_ml, group = participant_id)) +
  geom_line(alpha = 0.45, linewidth = 0.4, color = "grey40") +
  geom_point(alpha = 0.7, size = 1.4, color = "grey20") +
  scale_x_discrete(labels = c(pre = "pre", t0 = "t0", t4 = "t4")) +
  labs(x = NULL, y = "Osteocalcin (ng/mL)", tag = "A",
       title = "Per-participant trajectories") +
  theme_pub() +
  theme(plot.tag = element_text(face = "plain", size = 12))

# Panel B: reverse-ROPE effect-size-vs-probability curves with U95 markers.
x_hi  <- max(vapply(draws_by_contrast, function(d) quantile(abs(d), 0.999), numeric(1)))
xgrid <- seq(0, x_hi, length.out = 250)
curve_df <- purrr::map2_dfr(draws_by_contrast, names(draws_by_contrast), function(d, lv) {
  tibble::tibble(contrast = lv, x = xgrid,
                 p = vapply(xgrid, function(t) mean(abs(d) > t), numeric(1)))
})
u95_df <- rope_tbl |> dplyr::transmute(contrast, u95_ngml)

pB <- ggplot(curve_df, aes(x = x, y = p, color = contrast)) +
  geom_hline(yintercept = 0.95, color = "grey60", linetype = "dotted", linewidth = 0.3) +
  geom_vline(xintercept = noise_floor, color = "#56B4E9", linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = sesoi,       color = "grey45",  linetype = "dashed", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_vline(data = u95_df, aes(xintercept = u95_ngml, color = contrast),
             linetype = "solid", linewidth = 0.5, show.legend = FALSE) +
  annotate("text", x = noise_floor, y = 0.05, label = "assay\nnoise", hjust = -0.05,
           size = 2.5, color = "#56B4E9") +
  annotate("text", x = sesoi, y = 0.55, label = "0.2 SD", hjust = -0.1,
           size = 2.5, color = "grey45") +
  scale_color_manual(values = okabe, labels = lab_for, name = NULL) +
  labs(x = "Candidate effect size  |Δ osteocalcin|  (ng/mL)",
       y = "P( |Δ| > x )", tag = "B",
       title = "Reverse-ROPE: U95 = 95% upper bound on |change|") +
  theme_pub() +
  theme(legend.position = "top",
        plot.tag = element_text(face = "plain", size = 12))

fig <- patchwork::wrap_plots(pA, pB, ncol = 2, widths = c(1, 1.4))
ggsave(file.path(out_dir_fig, "walking-experiment-rope.pdf"),
       fig, width = 9, height = 4.0, device = cairo_pdf)

cat("\nSaved figure: outputs/figures/working/walking-experiment/walking-experiment-rope.pdf\n")
cat("Done: code/_experiments/walking-experiment-rope.R\n")
