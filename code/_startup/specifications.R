###############################################################################
# Central specifications registry. Sourced by `_startup/init.R` AFTER `dat` is
# loaded (specs reference data objects by string name, so the data must exist
# in the environment first). Consumed by every Tier-2 analysis script so the
# formula + family + priors + grid bounds for each model live in exactly one
# place. See `docs/briefs/2026-05-21-specifications-registry-architecture.md`.
#
# Project scope (confirmed with Ian 2026-05-25; see
# `docs/briefs/2026-05-23-ian-clean-dataset-email.md`,
# `docs/briefs/2026-05-25-ian-scope-confirmation-reply.md`):
#
#   Set 1 -- causal PA → bone. Two exposures (ENMO, daily steps) × three
#            outcomes (tibial SOS, CTX-1, osteocalcin) = 6 analyses. DAG-
#            derived adjustment set: age (in the tensor smooth) + sex +
#            functional status + pregnancy/lactation.
#
#   Set 2 -- industrialization → {tibial SOS, ENMO, daily steps}, age + sex
#            only = 3 analyses.
#
# = 9 analyses total.
###############################################################################


# ---- Shared MCMC config ----

WARMUP        <- 2000
ITER          <- 7000
THIN          <- 5
CHAINS        <- 10
SEED          <- 2138
BRMS_CONTROL  <- list(adapt_delta = 0.999, max_treedepth = 20)


# ---- Shared priors ----

# Outcomes used in PA → bone are rescaled in `data.R` so the prior scales are
# interpretable on the rescaled scale: tibial SOS / 1000 (intercept ~ 3.5),
# osteocalcin / 10000 (lognormal: log-intercept ~ 1.6), CTX-1 raw ng/mL
# (lognormal: log-intercept ~ -1.5). The student_t(3, 0, 2.5) family is wide
# enough to be uninformative on every one of these scales and matches the
# project's established prior convention.

# PA -> bone, mediator DAG (primary, main text). 5-var MSAS.
brm_priors_pa_bone_med <- c(
  set_prior("student_t(3, 0, 2.5)", class = "sds"),
  set_prior("student_t(3, 0, 2.5)", class = "b", coef = "sexmale"),
  set_prior("student_t(3, 0, 2.5)", class = "b", coef = "pregnant_or_breastfeeding_n_y_0_1"),
  set_prior("student_t(3, 0, 2.5)", class = "b", coef = "smoking_binary_n_y_0_1"),
  set_prior("student_t(3, 0, 2.5)", class = "b", coef = "alcohol_binary_n_y_0_1"),
  set_prior("student_t(3, 0, 2.5)", class = "b", coef = "functional_status_n_y_0_1")
)

# PA -> bone, confounder DAG (supplement). 6-var MSAS (mediator + body comp).
# Fat mass & lean body mass operationalized as fat_mass_kg_z + fat_free_mass_kg_z
# (per the 2026-05-19 issue 1 decision: bundled DAG node, two regression covariates).
brm_priors_pa_bone_conf <- c(
  brm_priors_pa_bone_med,
  set_prior("student_t(3, 0, 2.5)", class = "b", coef = "fat_mass_kg_z"),
  set_prior("student_t(3, 0, 2.5)", class = "b", coef = "fat_free_mass_kg_z")
)

# Industrialization analyses adjust for age + sex only (per Ian 2026-05-23).
brm_priors_urb <- c(
  set_prior("student_t(3, 0, 2.5)", class = "sds"),
  set_prior("student_t(3, 0, 2.5)", class = "b", coef = "sexmale")
)


# ---- Per-template model specifications ----
#
# Keyed by kebab-case script basename. Family lives INSIDE bf() so each spec
# is a single self-contained object. Data is referenced by string ("dat")
# rather than by value -- consumers do `get(spec$data)`. Each spec carries
# enough metadata (outcome_label, exposure_label, amef_label, curvature_label,
# outcome_scale_factor) for the shared `aerf_amef_plot()` helper to drive
# strip text + axis ticks without any per-script customization.
#
# `outcome_scale_factor` is the multiplier that brings pred_draws / slope_draws
# back to the natural reporting scale (m/s for SOS, raw step count for the
# step-count outcome, pg/mL for osteocalcin). Analyses multiply their
# pred_draws$draw and slope_draws$draw by this value before save so every
# downstream artifact is in natural units.

model_templates <- list(

  # -------------------------------------------------------------------------
  # Set 1 -- causal PA → bone (n = 6)
  # -------------------------------------------------------------------------

  "sos-steps" = list(
    bf = bf(
      tibia_sos_1k ~ t2(age_years, ad_steps_1k, k = c(5, 5)) +
                       sex + pregnant_or_breastfeeding_n_y_0_1 +
                       smoking_binary_n_y_0_1 +
                       alcohol_binary_n_y_0_1 +
                       functional_status_n_y_0_1 +
                       village_id,
      family = gaussian(link = "identity")
    ),
    outcome              = "tibia_sos_1k",
    outcome_label        = "Tibial SOS (m/s)",
    outcome_scale_factor = 1000,
    exposure             = "ad_steps_1k",
    exposure_label       = "Daily step count",
    amef_label           = "dSOS / dSteps",
    curvature_label      = "d^2 SOS / dSteps^2",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    priors               = brm_priors_pa_bone_med
  ),

  "sos-enmo" = list(
    bf = bf(
      tibia_sos_1k ~ t2(age_years, ad_mean_enmo_mg_0_24hr, k = c(5, 5)) +
                       sex + pregnant_or_breastfeeding_n_y_0_1 +
                       smoking_binary_n_y_0_1 +
                       alcohol_binary_n_y_0_1 +
                       functional_status_n_y_0_1 +
                       village_id,
      family = gaussian(link = "identity")
    ),
    outcome              = "tibia_sos_1k",
    outcome_label        = "Tibial SOS (m/s)",
    outcome_scale_factor = 1000,
    exposure             = "ad_mean_enmo_mg_0_24hr",
    exposure_label       = "Mean daily ENMO (mg)",
    amef_label           = "dSOS / dENMO",
    curvature_label      = "d^2 SOS / dENMO^2",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    priors               = brm_priors_pa_bone_med
  ),

  "ctx-steps" = list(
    bf = bf(
      ctx1_ng_ml ~ t2(age_years, ad_steps_1k, k = c(5, 5)) +
                     sex + pregnant_or_breastfeeding_n_y_0_1 +
                     smoking_binary_n_y_0_1 +
                     alcohol_binary_n_y_0_1 +
                     functional_status_n_y_0_1 +
                       village_id,
      family = lognormal(link = "identity", link_sigma = "log")
    ),
    outcome              = "ctx1_ng_ml",
    outcome_label        = "CTX-1 (ng/mL)",
    outcome_scale_factor = 1,
    exposure             = "ad_steps_1k",
    exposure_label       = "Daily step count",
    amef_label           = "dCTX / dSteps",
    curvature_label      = "d^2 CTX / dSteps^2",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    priors               = brm_priors_pa_bone_med
  ),

  "ctx-enmo" = list(
    bf = bf(
      ctx1_ng_ml ~ t2(age_years, ad_mean_enmo_mg_0_24hr, k = c(5, 5)) +
                     sex + pregnant_or_breastfeeding_n_y_0_1 +
                     smoking_binary_n_y_0_1 +
                     alcohol_binary_n_y_0_1 +
                     functional_status_n_y_0_1 +
                       village_id,
      family = lognormal(link = "identity", link_sigma = "log")
    ),
    outcome              = "ctx1_ng_ml",
    outcome_label        = "CTX-1 (ng/mL)",
    outcome_scale_factor = 1,
    exposure             = "ad_mean_enmo_mg_0_24hr",
    exposure_label       = "Mean daily ENMO (mg)",
    amef_label           = "dCTX / dENMO",
    curvature_label      = "d^2 CTX / dENMO^2",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    priors               = brm_priors_pa_bone_med
  ),

  "osteo-steps" = list(
    bf = bf(
      osteocalcin_pg_ml_10k ~ t2(age_years, ad_steps_1k, k = c(5, 5)) +
                                sex + pregnant_or_breastfeeding_n_y_0_1 +
                                smoking_binary_n_y_0_1 +
                                alcohol_binary_n_y_0_1 +
                                functional_status_n_y_0_1 +
                       village_id,
      family = lognormal(link = "identity", link_sigma = "log")
    ),
    outcome              = "osteocalcin_pg_ml_10k",
    outcome_label        = "Osteocalcin (pg/mL)",
    outcome_scale_factor = 10000,
    exposure             = "ad_steps_1k",
    exposure_label       = "Daily step count",
    amef_label           = "dOsteocalcin / dSteps",
    curvature_label      = "d^2 Osteocalcin / dSteps^2",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    priors               = brm_priors_pa_bone_med
  ),

  "osteo-enmo" = list(
    bf = bf(
      osteocalcin_pg_ml_10k ~ t2(age_years, ad_mean_enmo_mg_0_24hr, k = c(5, 5)) +
                                sex + pregnant_or_breastfeeding_n_y_0_1 +
                                smoking_binary_n_y_0_1 +
                                alcohol_binary_n_y_0_1 +
                                functional_status_n_y_0_1 +
                       village_id,
      family = lognormal(link = "identity", link_sigma = "log")
    ),
    outcome              = "osteocalcin_pg_ml_10k",
    outcome_label        = "Osteocalcin (pg/mL)",
    outcome_scale_factor = 10000,
    exposure             = "ad_mean_enmo_mg_0_24hr",
    exposure_label       = "Mean daily ENMO (mg)",
    amef_label           = "dOsteocalcin / dENMO",
    curvature_label      = "d^2 Osteocalcin / dENMO^2",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    priors               = brm_priors_pa_bone_med
  ),

  # -------------------------------------------------------------------------
  # Set 1b -- causal PA → bone, CONFOUNDER DAG (n = 6, supplement variants).
  # MSAS = mediator + Fat mass & lean body mass (fat_mass_kg_z + fat_free_mass_kg_z).
  # Per the 2026-05-20 dual-DAG settlement with Ian.
  # -------------------------------------------------------------------------

  "sos-steps-conf" = list(
    bf = bf(
      tibia_sos_1k ~ t2(age_years, ad_steps_1k, k = c(5, 5)) +
                       sex + pregnant_or_breastfeeding_n_y_0_1 +
                       smoking_binary_n_y_0_1 +
                       alcohol_binary_n_y_0_1 +
                       functional_status_n_y_0_1 +
                       fat_mass_kg_z + fat_free_mass_kg_z +
                       village_id,
      family = gaussian(link = "identity")
    ),
    outcome              = "tibia_sos_1k",
    outcome_label        = "Tibial SOS (m/s)",
    outcome_scale_factor = 1000,
    exposure             = "ad_steps_1k",
    exposure_label       = "Daily step count",
    amef_label           = "dSOS / dSteps",
    curvature_label      = "d^2 SOS / dSteps^2",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    priors               = brm_priors_pa_bone_conf
  ),

  "sos-enmo-conf" = list(
    bf = bf(
      tibia_sos_1k ~ t2(age_years, ad_mean_enmo_mg_0_24hr, k = c(5, 5)) +
                       sex + pregnant_or_breastfeeding_n_y_0_1 +
                       smoking_binary_n_y_0_1 +
                       alcohol_binary_n_y_0_1 +
                       functional_status_n_y_0_1 +
                       fat_mass_kg_z + fat_free_mass_kg_z +
                       village_id,
      family = gaussian(link = "identity")
    ),
    outcome              = "tibia_sos_1k",
    outcome_label        = "Tibial SOS (m/s)",
    outcome_scale_factor = 1000,
    exposure             = "ad_mean_enmo_mg_0_24hr",
    exposure_label       = "Mean daily ENMO (mg)",
    amef_label           = "dSOS / dENMO",
    curvature_label      = "d^2 SOS / dENMO^2",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    priors               = brm_priors_pa_bone_conf
  ),

  "ctx-steps-conf" = list(
    bf = bf(
      ctx1_ng_ml ~ t2(age_years, ad_steps_1k, k = c(5, 5)) +
                     sex + pregnant_or_breastfeeding_n_y_0_1 +
                     smoking_binary_n_y_0_1 +
                     alcohol_binary_n_y_0_1 +
                     functional_status_n_y_0_1 +
                     fat_mass_kg_z + fat_free_mass_kg_z +
                       village_id,
      family = lognormal(link = "identity", link_sigma = "log")
    ),
    outcome              = "ctx1_ng_ml",
    outcome_label        = "CTX-1 (ng/mL)",
    outcome_scale_factor = 1,
    exposure             = "ad_steps_1k",
    exposure_label       = "Daily step count",
    amef_label           = "dCTX / dSteps",
    curvature_label      = "d^2 CTX / dSteps^2",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    priors               = brm_priors_pa_bone_conf
  ),

  "ctx-enmo-conf" = list(
    bf = bf(
      ctx1_ng_ml ~ t2(age_years, ad_mean_enmo_mg_0_24hr, k = c(5, 5)) +
                     sex + pregnant_or_breastfeeding_n_y_0_1 +
                     smoking_binary_n_y_0_1 +
                     alcohol_binary_n_y_0_1 +
                     functional_status_n_y_0_1 +
                     fat_mass_kg_z + fat_free_mass_kg_z +
                       village_id,
      family = lognormal(link = "identity", link_sigma = "log")
    ),
    outcome              = "ctx1_ng_ml",
    outcome_label        = "CTX-1 (ng/mL)",
    outcome_scale_factor = 1,
    exposure             = "ad_mean_enmo_mg_0_24hr",
    exposure_label       = "Mean daily ENMO (mg)",
    amef_label           = "dCTX / dENMO",
    curvature_label      = "d^2 CTX / dENMO^2",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    priors               = brm_priors_pa_bone_conf
  ),

  "osteo-steps-conf" = list(
    bf = bf(
      osteocalcin_pg_ml_10k ~ t2(age_years, ad_steps_1k, k = c(5, 5)) +
                                sex + pregnant_or_breastfeeding_n_y_0_1 +
                                smoking_binary_n_y_0_1 +
                                alcohol_binary_n_y_0_1 +
                                functional_status_n_y_0_1 +
                                fat_mass_kg_z + fat_free_mass_kg_z +
                       village_id,
      family = lognormal(link = "identity", link_sigma = "log")
    ),
    outcome              = "osteocalcin_pg_ml_10k",
    outcome_label        = "Osteocalcin (pg/mL)",
    outcome_scale_factor = 10000,
    exposure             = "ad_steps_1k",
    exposure_label       = "Daily step count",
    amef_label           = "dOsteocalcin / dSteps",
    curvature_label      = "d^2 Osteocalcin / dSteps^2",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    priors               = brm_priors_pa_bone_conf
  ),

  "osteo-enmo-conf" = list(
    bf = bf(
      osteocalcin_pg_ml_10k ~ t2(age_years, ad_mean_enmo_mg_0_24hr, k = c(5, 5)) +
                                sex + pregnant_or_breastfeeding_n_y_0_1 +
                                smoking_binary_n_y_0_1 +
                                alcohol_binary_n_y_0_1 +
                                functional_status_n_y_0_1 +
                                fat_mass_kg_z + fat_free_mass_kg_z +
                       village_id,
      family = lognormal(link = "identity", link_sigma = "log")
    ),
    outcome              = "osteocalcin_pg_ml_10k",
    outcome_label        = "Osteocalcin (pg/mL)",
    outcome_scale_factor = 10000,
    exposure             = "ad_mean_enmo_mg_0_24hr",
    exposure_label       = "Mean daily ENMO (mg)",
    amef_label           = "dOsteocalcin / dENMO",
    curvature_label      = "d^2 Osteocalcin / dENMO^2",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    priors               = brm_priors_pa_bone_conf
  ),

  # -------------------------------------------------------------------------
  # Set 2 -- industrialization (n = 3). Adjustment set: age + sex only.
  # -------------------------------------------------------------------------

  "sos-urb" = list(
    bf = bf(
      tibia_sos_1k ~ t2(age_years, industrial_index, k = c(5, 5)) + sex,
      family = gaussian(link = "identity")
    ),
    outcome              = "tibia_sos_1k",
    outcome_label        = "Tibial SOS (m/s)",
    outcome_scale_factor = 1000,
    exposure             = "industrial_index",
    exposure_label       = "Industrialization index",
    amef_label           = "dSOS / dIndustrialization",
    curvature_label      = "d^2 SOS / dIndustrialization^2",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    priors               = brm_priors_urb
  ),

  "steps-urb" = list(
    bf = bf(
      ad_steps_1k ~ t2(age_years, industrial_index, k = c(5, 5)) + sex,
      family = gaussian(link = "identity")
    ),
    outcome              = "ad_steps_1k",
    outcome_label        = "Daily step count",
    outcome_scale_factor = 1000,
    exposure             = "industrial_index",
    exposure_label       = "Industrialization index",
    amef_label           = "dSteps / dIndustrialization",
    curvature_label      = "d^2 Steps / dIndustrialization^2",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    priors               = brm_priors_urb
  ),

  "enmo-urb" = list(
    bf = bf(
      ad_mean_enmo_mg_0_24hr ~ t2(age_years, industrial_index, k = c(5, 5)) + sex,
      family = gaussian(link = "identity")
    ),
    outcome              = "ad_mean_enmo_mg_0_24hr",
    outcome_label        = "Mean daily\nENMO (mg)",
    outcome_scale_factor = 1,
    exposure             = "industrial_index",
    exposure_label       = "Industrialization index",
    amef_label           = "dENMO / dIndustrialization",
    curvature_label      = "d^2 ENMO / dIndustrialization^2",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    priors               = brm_priors_urb
  )
)
