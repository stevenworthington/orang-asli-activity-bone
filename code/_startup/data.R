###############################################################################
# Read + prep the working dataset (`dat`). The OA HeLP physical-activity / bone
# CSV: n = 1,007, with one canonical `industrial_index` column. Factors `sex`
# and `village_id`.
#
# SCALING. Exactly one derived scaling per variable, used wherever that variable
# appears — as an outcome in one analysis and as an exposure in another. Tibial
# SOS /200, daily steps /5,000, mean ENMO /10, osteocalcin /10,000. CTX-1 and
# the industrialization index are left raw.
#
# Why these divisors. Priors on the coefficient block are autoscaled
# (normal(0, kappa * sd_lp / sd(x_j))), but the `sds` prior on smooth wiggliness
# is not: brms offers no autoscaling for it. A single exponential(2) on `sds` is
# only comparable across outcomes when they sit on a common scale.
#
# That argument applies to the IDENTITY-LINK outcomes only: tibial SOS /200,
# steps /5,000 and ENMO /10 each land near sd ~ 1 on the scale the smooth is
# built on. It does not apply to osteocalcin. Under the lognormal likelihood the
# linear predictor is log(y), and log(y / 10,000) = log(y) - log(10,000): the
# divisor shifts the intercept and leaves sd(log y) at ~0.50 whatever it is.
# Osteocalcin /10,000 is kept for numerical conditioning of the intercept, not
# for a common smooth-prior scale. CTX-1, also log-link, is left raw for the
# same reason -- there is nothing to gain.
#
# Why one scaling per variable rather than one per role. A linear rescaling of
# an exposure inside a penalized smooth is very nearly free: mgcv::smoothCon
# rescales the penalty relative to the model matrix (scale.penalty = TRUE), so
# halving or doubling the divisor leaves the reported bounds effectively
# unchanged.
#
# The raw source columns (tibia_sos, ad_tot_step_count_0_24hr,
# ad_mean_enmo_mg_0_24hr, osteocalcin_pg_ml) remain in `dat` because they come
# from the CSV; nothing models them directly.
###############################################################################


# The derivation itself lives in `prep_dat()` in `_startup/functions.R`, called
# from here and from the `dat` target in `_targets.R`, so both paths derive the
# columns identically.
dat <- prep_dat(here::here("data", "processed", "orang-asli-pa-bone-data.csv"))
