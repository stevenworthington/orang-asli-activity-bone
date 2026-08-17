###############################################################################
# Central specifications registry: the formula, family, prior settings and grid
# bounds for each model, in exactly one place. Sourced by `_startup/init.R` and,
# without data.R, by `_targets.R`.
#
# It does NOT require `dat` to exist. Formula columns are symbols inside
# brms::bf() and are not evaluated until a fit; the data is named by the string
# "dat", which consumers resolve with get(spec$data).
#
# Project scope:
#
#   Causal PA → bone. Two exposures (ENMO, daily steps) × three
#            outcomes (tibial SOS, CTX-1, osteocalcin) = 6 analyses. DAG-
#            derived adjustment set: age (in the tensor smooth) + sex +
#            pregnancy/lactation + smoking + alcohol + functional status,
#            plus community identity as a fixed effect.
#
#   Industrialization → {tibial SOS, ENMO, daily steps}, age + sex
#            only = 3 analyses.
#
# = 9 reported analyses. The registry also carries 6 confounder-DAG sensitivity
# variants of the PA → bone specs (below), so 15 specs are fitted in total.
#
#=============================================================================
# LIKELIHOODS. One per outcome, set in the `bf()` call of each spec below. For
# the PA -> bone analyses the family was selected per outcome on posterior
# predictive grounds rather than left Gaussian by default:
#   tibial SOS     Student-t   heavy lower tail
#   CTX-1          Gamma(log)  positive support, right skew
#   osteocalcin    lognormal
# The industrialization likelihood is NOT a posterior-predictive choice:
#   industrializ.  Gaussian    identity link, on estimand grounds -- see below
#
# INDUSTRIALIZATION FORM. One-stage with a community random intercept and
# ADDITIVE smooths. Three reasons: the index is community-constant, so an age x
# index interaction is informed by only ~25 distinct index values -- identified,
# but too weakly for `t2` to buy anything the design can support; the additive
# multilevel form is the one whose operating characteristics were measured; and
# the identity link matches the reported estimand (steps per index UNIT -- a log
# link would make it multiplicative). It is deliberately not chosen for fit: a
# log link scores better on the posterior predictive battery, at the cost of
# changing what the reported effect means.
#
# PRIORS. Autoscaled at fit time by `build_priors()` in `_startup/functions.R`
# -- there is no static `priors` field, because the scales depend on the
# analytic data. Each spec declares only `prior_kappa` and, for the
# community-level specs, `prior_sd_community`.
#
#   PA -> bone            kappa 0.25, community indicators FLAT
#   industrialization     kappa 0.5,  community SD prior mean 0.33 * sd_lp
#
# Both settings have measured operating characteristics. The harnesses that
# measure them are `_experiments/calibration-onestage-production.R` (whose
# defaults are the settings used here) and `_experiments/calibration-pa-bone.R`.
# Each simulates on the real study design and reports false-positive rate,
# interval coverage and power for the linear-projection estimand the analyses
# report. Run them rather than relying on a number quoted here: a figure copied
# into a comment cannot be re-derived and goes stale silently.
#
# SCALING. One derived scaling per variable, set in `data.R`: tibial SOS /200,
# daily steps /5,000, mean ENMO /10, osteocalcin /10,000. `contrast_units` is
# how many of those units make the reported contrast, so it is 1 wherever the
# divisor was chosen to equal the contrast.
#=============================================================================
###############################################################################


# ---- Shared MCMC config ----

WARMUP        <- 2000
ITER          <- 7000
THIN          <- 5
CHAINS        <- 10
SEED          <- 2138
BRMS_CONTROL  <- list(adapt_delta = 0.999, max_treedepth = 20)

# Smoke-test sampling: plumbing only, never a reported number. Defined once and
# used by BOTH paths -- fit_one_spec(mode = "smoke") and the nine analysis
# scripts' SAMPLING_MODE=smoke branch -- so the two cannot drift apart.
SMOKE_WARMUP  <- 500
SMOKE_ITER    <- 1500
SMOKE_CHAINS  <- 4


# ---- Prior settings ----

# Consumed by build_priors(); see its header for the full scheme and for why
# the PA community indicators are left flat.
KAPPA_PA           <- 0.25
KAPPA_URB          <- 0.5
SD_COMMUNITY_URB   <- 0.33


# ---- Per-template model specifications ----
#
# Keyed by kebab-case script basename. Family lives INSIDE bf() so each spec
# is a single self-contained object. Data is referenced by string ("dat")
# rather than by value -- consumers do `get(spec$data)`. Each spec carries
# enough metadata (outcome_label, exposure_label, amef_label, curvature_label,
# outcome_scale_factor, contrast_units, contrast_label) for the shared plotting
# and table helpers to drive strip text + axis ticks + reported contrasts
# without any per-script customization.
#
# `outcome_scale_factor` is the multiplier that brings pred_draws / slope_draws
# back to the natural reporting scale (m/s for SOS, raw step count for the
# step-count outcome, pg/mL for osteocalcin). Analyses multiply their
# pred_draws$draw and slope_draws$draw by this value before save so every
# downstream artifact is in natural units.
#
# `contrast_units` is the width of the reported contrast in EXPOSURE units.
# Multiply a per-unit slope by it to get the reported effect. Because data.R
# scales daily steps by 5,000 and ENMO by 10 -- the two reported contrasts --
# it is 1 for every PA spec. The industrialization index is unscaled and
# reported per 10 index units, so it is 10 there. A slope multiplied by the
# wrong value is wrong by a factor of 5 or 10 with no error raised, so this
# field exists to keep the number in one place.
#
# WHERE IT IS APPLIED. Every SAVED artifact -- the targets draws,
# `linear_projection()`'s beta, `summarize_one_fit()`'s inputs -- stays per ONE
# exposure unit. The conversion happens at the presentation boundary:
# `make_amef_panel()` and `make_age_cond_amef_panel_with_bands()` for figures,
# the two `_final/` table scripts at their own call sites. One conversion per
# consumer, so nothing squares the contrast by applying it twice.
#
# `re_formula` is passed through to postestimation when a spec sets it. The
# three industrialization specs set NA, which makes predictions marginal of the
# community random intercept -- the across-community estimand the manuscript
# claims. Without it marginaleffects takes brms's default (NULL), conditioning
# on the fitted community offsets: the reported quantities agree to 12
# significant figures either way, but the BAND is visibly narrower, because the
# population intercept and the mean community offset are anticorrelated. The PA
# specs have no random effects (community enters as fixed indicators) and so
# leave the field unset.

model_templates <- list(

  # -------------------------------------------------------------------------
  # Causal PA → bone, MEDIATOR DAG (n = 6; primary, main text)
  # -------------------------------------------------------------------------

  "sos-steps" = list(
    bf = brms::bf(
      tibia_sos_200 ~ t2(age_years, ad_steps_5k, k = c(5, 5)) +
                       sex + pregnant_or_breastfeeding_n_y_0_1 +
                       smoking_binary_n_y_0_1 +
                       alcohol_binary_n_y_0_1 +
                       functional_status_n_y_0_1 +
                       village_id,
      family = brms::student()
    ),
    outcome              = "tibia_sos_200",
    outcome_label        = "Tibial SOS (m/s)",
    outcome_scale_factor = 200,
    exposure             = "ad_steps_5k",
    exposure_label       = "Daily step count",
    amef_label           = "dSOS / dSteps",
    curvature_label      = "d^2 SOS / dSteps^2",
    contrast_units       = 1,
    contrast_label       = "per 5,000 steps",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    prior_kappa          = KAPPA_PA
  ),

  "sos-enmo" = list(
    bf = brms::bf(
      tibia_sos_200 ~ t2(age_years, enmo_10, k = c(5, 5)) +
                       sex + pregnant_or_breastfeeding_n_y_0_1 +
                       smoking_binary_n_y_0_1 +
                       alcohol_binary_n_y_0_1 +
                       functional_status_n_y_0_1 +
                       village_id,
      family = brms::student()
    ),
    outcome              = "tibia_sos_200",
    outcome_label        = "Tibial SOS (m/s)",
    outcome_scale_factor = 200,
    exposure             = "enmo_10",
    exposure_label       = "Mean daily ENMO (mg)",
    amef_label           = "dSOS / dENMO",
    curvature_label      = "d^2 SOS / dENMO^2",
    contrast_units       = 1,
    contrast_label       = "per 10 mg ENMO",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    prior_kappa          = KAPPA_PA
  ),

  "ctx-steps" = list(
    bf = brms::bf(
      ctx1_ng_ml ~ t2(age_years, ad_steps_5k, k = c(5, 5)) +
                     sex + pregnant_or_breastfeeding_n_y_0_1 +
                     smoking_binary_n_y_0_1 +
                     alcohol_binary_n_y_0_1 +
                     functional_status_n_y_0_1 +
                     village_id,
      family = Gamma(link = "log")
    ),
    outcome              = "ctx1_ng_ml",
    outcome_label        = "CTX-1 (ng/mL)",
    outcome_scale_factor = 1,
    exposure             = "ad_steps_5k",
    exposure_label       = "Daily step count",
    amef_label           = "dCTX / dSteps",
    curvature_label      = "d^2 CTX / dSteps^2",
    contrast_units       = 1,
    contrast_label       = "per 5,000 steps",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    prior_kappa          = KAPPA_PA
  ),

  "ctx-enmo" = list(
    bf = brms::bf(
      ctx1_ng_ml ~ t2(age_years, enmo_10, k = c(5, 5)) +
                     sex + pregnant_or_breastfeeding_n_y_0_1 +
                     smoking_binary_n_y_0_1 +
                     alcohol_binary_n_y_0_1 +
                     functional_status_n_y_0_1 +
                     village_id,
      family = Gamma(link = "log")
    ),
    outcome              = "ctx1_ng_ml",
    outcome_label        = "CTX-1 (ng/mL)",
    outcome_scale_factor = 1,
    exposure             = "enmo_10",
    exposure_label       = "Mean daily ENMO (mg)",
    amef_label           = "dCTX / dENMO",
    curvature_label      = "d^2 CTX / dENMO^2",
    contrast_units       = 1,
    contrast_label       = "per 10 mg ENMO",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    prior_kappa          = KAPPA_PA
  ),

  "osteo-steps" = list(
    bf = brms::bf(
      osteocalcin_pg_ml_10k ~ t2(age_years, ad_steps_5k, k = c(5, 5)) +
                                sex + pregnant_or_breastfeeding_n_y_0_1 +
                                smoking_binary_n_y_0_1 +
                                alcohol_binary_n_y_0_1 +
                                functional_status_n_y_0_1 +
                                village_id,
      family = brms::lognormal(link = "identity", link_sigma = "log")
    ),
    outcome              = "osteocalcin_pg_ml_10k",
    outcome_label        = "Osteocalcin (pg/mL)",
    outcome_scale_factor = 10000,
    exposure             = "ad_steps_5k",
    exposure_label       = "Daily step count",
    amef_label           = "dOsteocalcin / dSteps",
    curvature_label      = "d^2 Osteocalcin / dSteps^2",
    contrast_units       = 1,
    contrast_label       = "per 5,000 steps",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    prior_kappa          = KAPPA_PA
  ),

  "osteo-enmo" = list(
    bf = brms::bf(
      osteocalcin_pg_ml_10k ~ t2(age_years, enmo_10, k = c(5, 5)) +
                                sex + pregnant_or_breastfeeding_n_y_0_1 +
                                smoking_binary_n_y_0_1 +
                                alcohol_binary_n_y_0_1 +
                                functional_status_n_y_0_1 +
                                village_id,
      family = brms::lognormal(link = "identity", link_sigma = "log")
    ),
    outcome              = "osteocalcin_pg_ml_10k",
    outcome_label        = "Osteocalcin (pg/mL)",
    outcome_scale_factor = 10000,
    exposure             = "enmo_10",
    exposure_label       = "Mean daily ENMO (mg)",
    amef_label           = "dOsteocalcin / dENMO",
    curvature_label      = "d^2 Osteocalcin / dENMO^2",
    contrast_units       = 1,
    contrast_label       = "per 10 mg ENMO",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    prior_kappa          = KAPPA_PA
  ),

  # -------------------------------------------------------------------------
  # Causal PA → bone, CONFOUNDER DAG (n = 6, supplementary variants).
  # MSAS = mediator + Fat mass & lean body mass (fat_mass_kg_z + fat_free_mass_kg_z).
  # -------------------------------------------------------------------------

  "sos-steps-conf" = list(
    bf = brms::bf(
      tibia_sos_200 ~ t2(age_years, ad_steps_5k, k = c(5, 5)) +
                       sex + pregnant_or_breastfeeding_n_y_0_1 +
                       smoking_binary_n_y_0_1 +
                       alcohol_binary_n_y_0_1 +
                       functional_status_n_y_0_1 +
                       fat_mass_kg_z + fat_free_mass_kg_z +
                       village_id,
      family = brms::student()
    ),
    outcome              = "tibia_sos_200",
    outcome_label        = "Tibial SOS (m/s)",
    outcome_scale_factor = 200,
    exposure             = "ad_steps_5k",
    exposure_label       = "Daily step count",
    amef_label           = "dSOS / dSteps",
    curvature_label      = "d^2 SOS / dSteps^2",
    contrast_units       = 1,
    contrast_label       = "per 5,000 steps",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    prior_kappa          = KAPPA_PA
  ),

  "sos-enmo-conf" = list(
    bf = brms::bf(
      tibia_sos_200 ~ t2(age_years, enmo_10, k = c(5, 5)) +
                       sex + pregnant_or_breastfeeding_n_y_0_1 +
                       smoking_binary_n_y_0_1 +
                       alcohol_binary_n_y_0_1 +
                       functional_status_n_y_0_1 +
                       fat_mass_kg_z + fat_free_mass_kg_z +
                       village_id,
      family = brms::student()
    ),
    outcome              = "tibia_sos_200",
    outcome_label        = "Tibial SOS (m/s)",
    outcome_scale_factor = 200,
    exposure             = "enmo_10",
    exposure_label       = "Mean daily ENMO (mg)",
    amef_label           = "dSOS / dENMO",
    curvature_label      = "d^2 SOS / dENMO^2",
    contrast_units       = 1,
    contrast_label       = "per 10 mg ENMO",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    prior_kappa          = KAPPA_PA
  ),

  "ctx-steps-conf" = list(
    bf = brms::bf(
      ctx1_ng_ml ~ t2(age_years, ad_steps_5k, k = c(5, 5)) +
                     sex + pregnant_or_breastfeeding_n_y_0_1 +
                     smoking_binary_n_y_0_1 +
                     alcohol_binary_n_y_0_1 +
                     functional_status_n_y_0_1 +
                     fat_mass_kg_z + fat_free_mass_kg_z +
                     village_id,
      family = Gamma(link = "log")
    ),
    outcome              = "ctx1_ng_ml",
    outcome_label        = "CTX-1 (ng/mL)",
    outcome_scale_factor = 1,
    exposure             = "ad_steps_5k",
    exposure_label       = "Daily step count",
    amef_label           = "dCTX / dSteps",
    curvature_label      = "d^2 CTX / dSteps^2",
    contrast_units       = 1,
    contrast_label       = "per 5,000 steps",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    prior_kappa          = KAPPA_PA
  ),

  "ctx-enmo-conf" = list(
    bf = brms::bf(
      ctx1_ng_ml ~ t2(age_years, enmo_10, k = c(5, 5)) +
                     sex + pregnant_or_breastfeeding_n_y_0_1 +
                     smoking_binary_n_y_0_1 +
                     alcohol_binary_n_y_0_1 +
                     functional_status_n_y_0_1 +
                     fat_mass_kg_z + fat_free_mass_kg_z +
                     village_id,
      family = Gamma(link = "log")
    ),
    outcome              = "ctx1_ng_ml",
    outcome_label        = "CTX-1 (ng/mL)",
    outcome_scale_factor = 1,
    exposure             = "enmo_10",
    exposure_label       = "Mean daily ENMO (mg)",
    amef_label           = "dCTX / dENMO",
    curvature_label      = "d^2 CTX / dENMO^2",
    contrast_units       = 1,
    contrast_label       = "per 10 mg ENMO",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    prior_kappa          = KAPPA_PA
  ),

  "osteo-steps-conf" = list(
    bf = brms::bf(
      osteocalcin_pg_ml_10k ~ t2(age_years, ad_steps_5k, k = c(5, 5)) +
                                sex + pregnant_or_breastfeeding_n_y_0_1 +
                                smoking_binary_n_y_0_1 +
                                alcohol_binary_n_y_0_1 +
                                functional_status_n_y_0_1 +
                                fat_mass_kg_z + fat_free_mass_kg_z +
                                village_id,
      family = brms::lognormal(link = "identity", link_sigma = "log")
    ),
    outcome              = "osteocalcin_pg_ml_10k",
    outcome_label        = "Osteocalcin (pg/mL)",
    outcome_scale_factor = 10000,
    exposure             = "ad_steps_5k",
    exposure_label       = "Daily step count",
    amef_label           = "dOsteocalcin / dSteps",
    curvature_label      = "d^2 Osteocalcin / dSteps^2",
    contrast_units       = 1,
    contrast_label       = "per 5,000 steps",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    prior_kappa          = KAPPA_PA
  ),

  "osteo-enmo-conf" = list(
    bf = brms::bf(
      osteocalcin_pg_ml_10k ~ t2(age_years, enmo_10, k = c(5, 5)) +
                                sex + pregnant_or_breastfeeding_n_y_0_1 +
                                smoking_binary_n_y_0_1 +
                                alcohol_binary_n_y_0_1 +
                                functional_status_n_y_0_1 +
                                fat_mass_kg_z + fat_free_mass_kg_z +
                                village_id,
      family = brms::lognormal(link = "identity", link_sigma = "log")
    ),
    outcome              = "osteocalcin_pg_ml_10k",
    outcome_label        = "Osteocalcin (pg/mL)",
    outcome_scale_factor = 10000,
    exposure             = "enmo_10",
    exposure_label       = "Mean daily ENMO (mg)",
    amef_label           = "dOsteocalcin / dENMO",
    curvature_label      = "d^2 Osteocalcin / dENMO^2",
    contrast_units       = 1,
    contrast_label       = "per 10 mg ENMO",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    prior_kappa          = KAPPA_PA
  ),

  # -------------------------------------------------------------------------
  # Industrialization (n = 3). Adjustment set: age + sex only.
  #
  # One-stage, community random intercept, ADDITIVE smooths, Gaussian identity
  # link -- see the specification block at the top of this file. The exposure
  # smooth gets k = 4 rather than 5: the index varies only between communities,
  # so it carries ~25 distinct values, not n.
  # -------------------------------------------------------------------------

  "sos-urb" = list(
    bf = brms::bf(
      tibia_sos_200 ~ s(age_years, k = 5) + s(industrial_index, k = 4) +
                        sex + (1 | village_id),
      family = gaussian(link = "identity")
    ),
    outcome              = "tibia_sos_200",
    outcome_label        = "Tibial SOS (m/s)",
    outcome_scale_factor = 200,
    exposure             = "industrial_index",
    exposure_label       = "Industrialization index",
    amef_label           = "dSOS / dIndustrialization",
    curvature_label      = "d^2 SOS / dIndustrialization^2",
    contrast_units       = 10,
    contrast_label       = "per 10 index units",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    # The index is a COMMUNITY attribute, constant within village_id, so its
    # quantiles must be taken over the distinct community values. See
    # spec_grid() in pipeline-helpers.R for what goes wrong otherwise.
    grid_cluster         = "village_id",
    grid_n               = 41,
    prior_kappa          = KAPPA_URB,
    prior_sd_community   = SD_COMMUNITY_URB,
    re_formula           = NA
  ),

  "steps-urb" = list(
    bf = brms::bf(
      ad_steps_5k ~ s(age_years, k = 5) + s(industrial_index, k = 4) +
                      sex + (1 | village_id),
      family = gaussian(link = "identity")
    ),
    outcome              = "ad_steps_5k",
    outcome_label        = "Daily step count",
    outcome_scale_factor = 5000,
    exposure             = "industrial_index",
    exposure_label       = "Industrialization index",
    amef_label           = "dSteps / dIndustrialization",
    curvature_label      = "d^2 Steps / dIndustrialization^2",
    contrast_units       = 10,
    contrast_label       = "per 10 index units",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    # The index is a COMMUNITY attribute, constant within village_id, so its
    # quantiles must be taken over the distinct community values. See
    # spec_grid() in pipeline-helpers.R for what goes wrong otherwise.
    grid_cluster         = "village_id",
    grid_n               = 41,
    prior_kappa          = KAPPA_URB,
    prior_sd_community   = SD_COMMUNITY_URB,
    re_formula           = NA
  ),

  "enmo-urb" = list(
    bf = brms::bf(
      enmo_10 ~ s(age_years, k = 5) + s(industrial_index, k = 4) +
                  sex + (1 | village_id),
      family = gaussian(link = "identity")
    ),
    outcome              = "enmo_10",
    outcome_label        = "Mean daily\nENMO (mg)",
    outcome_scale_factor = 10,
    exposure             = "industrial_index",
    exposure_label       = "Industrialization index",
    amef_label           = "dENMO / dIndustrialization",
    curvature_label      = "d^2 ENMO / dIndustrialization^2",
    contrast_units       = 10,
    contrast_label       = "per 10 index units",
    data                 = "dat",
    grid_quantiles       = c(0.01, 0.99),
    # The index is a COMMUNITY attribute, constant within village_id, so its
    # quantiles must be taken over the distinct community values. See
    # spec_grid() in pipeline-helpers.R for what goes wrong otherwise.
    grid_cluster         = "village_id",
    grid_n               = 41,
    prior_kappa          = KAPPA_URB,
    prior_sd_community   = SD_COMMUNITY_URB,
    re_formula           = NA
  )
)


# ---- Natural units per unit of each modelled exposure column ----
#
# Used by `x_scale_for()` in `_startup/functions.R` to draw x-axis ticks in the
# units a reader thinks in, whatever divisor `data.R` applied. It lives beside
# `contrast_units` because both answer "how does this column relate to the
# quantity we report", and both must be updated together when a scaling moves.
#
# An explicit lookup rather than a name-matching test inside the plotting
# helper: a `grepl()` on the column name goes silently FALSE the moment a
# divisor changes and the column is renamed with it, and the axis then falls
# through to raw modelling units -- a step axis reading "1, 2, 3" where it means
# 5,000 / 10,000 / 15,000 steps, with nothing raised. So the mapping is stated
# here and `x_scale_for()` STOPS on an unregistered exposure rather than
# defaulting: a one-line chore when an exposure is added, against a mislabelled
# published figure.
EXPOSURE_AXIS_SCALE <- c(
  ad_steps_5k      = 5000,   # 1 unit = 5,000 steps
  enmo_10          = 10,     # 1 unit = 10 mg
  industrial_index = 1       # already in index units
)


# ---- Which specs are community-level ----

URB_KEYS <- c("sos-urb", "steps-urb", "enmo-urb")

is_urb_spec <- function(key) key %in% URB_KEYS
