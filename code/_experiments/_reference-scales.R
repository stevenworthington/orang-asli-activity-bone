###############################################################################
# Reference scales the effect-magnitude bounds are read against.
#
# A bound in raw outcome units means nothing on its own -- "we can rule out
# effects larger than 15 m/s" is only interpretable beside something else.
# These are the two yardsticks the manuscript uses:
#
#   NOISE_FLOOR    the measurement-resolution floor, so a bound can be compared
#                  with what the instrument could detect at all
#   analytic_sd()  the between-person SD of the spec's own analytic sample, so a
#                  bound can be quoted as a fraction of natural variation
#
# Sourced by the three scripts that report those bounds --
# `pa-contrast-effects.R`, `conf-slope-summaries.R` and
# `effect-size-probability-fig.R`. Defined once here because two copies of a
# reference scale is how two scripts end up quoting different denominators for
# the same bound.
#
# The constants below need nothing else loaded, and `analytic_sd()` touches the
# data only when called -- which lets the one consumer that runs without the
# startup set still share the floors.
#
# Deliberately NOT in `_startup/specifications.R`: `_targets.R` takes a
# `format = "file"` dependency on that file, so anything added to it invalidates
# all 15 fits. Nothing here is a modelling choice.
###############################################################################


# Measurement-resolution floors, in each outcome's natural units. External
# constants from the Methods -- for tibial SOS the IN-STUDY test-retest SEM,
# derived from the measured repeatability ICC, the total analytical CV (intra +
# inter, in quadrature) at the assay median for the two biomarkers, and the
# accelerometer step-count accuracy. Mean daily ENMO has no clean instrument
# floor, hence NA. None of these is derivable from the data.
#
# The tibial SOS floor is NOT the quantitative-ultrasound manufacturer's stated
# accuracy. That is a different and smaller quantity. Do not "correct" it back.
NOISE_FLOOR <- c(sos = 31.6, ctx = 0.008, osteo = 990, steps = 1200, enmo = NA)

UNIT_LABEL  <- c(sos = "m/s", ctx = "ng/mL", osteo = "pg/mL",
                 steps = "steps", enmo = "mg")

# The natural-unit column each outcome's SD is taken from.
OUTCOME_COL <- c(sos = "tibia_sos", ctx = "ctx1_ng_ml", osteo = "osteocalcin_pg_ml",
                 steps = "ad_tot_step_count_0_24hr", enmo = "ad_mean_enmo_mg_0_24hr")

# Outcome key from a registry spec key: "osteo-steps-conf" -> "osteo".
outcome_key <- function(spec_key) sub("-.*", "", spec_key)

# Between-person SD on the spec's OWN ANALYTIC SAMPLE -- the rows the estimate is
# computed from, after listwise deletion on that spec's formula. Computed rather
# than hard-coded, so it tracks the data file.
#
# Analytic rather than cohort: the bound is estimated on the analytic sample and
# the G-computation averages over that sample's covariate distribution, so the
# yardstick has to be drawn from the same population the estimate describes.
#
# Spec-keyed rather than outcome-keyed, because the mediator and confounder
# variants of one outcome delete different rows, and a single per-outcome
# constant cannot express that. Keeping it to one denominator is the point: two
# copies of a reference scale is how two scripts end up quoting different
# fractions for the same bound.
#
# Being a function, this needs `dat`, `model_templates` and `prep_local_data()`
# only when it is CALLED -- so this file can still be sourced for the floors and
# labels above without the startup set loaded.
analytic_sd <- function(spec_key) {
  spec <- model_templates[[spec_key]]
  d    <- prep_local_data(spec, dat)
  stats::sd(d[[OUTCOME_COL[[outcome_key(spec_key)]]]], na.rm = TRUE)
}
