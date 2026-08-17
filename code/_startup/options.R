###############################################################################
# Global options for the session: brms / marginaleffects defaults on top of
# base-R options, plus the one presentation switch the figure scripts read.
###############################################################################


# `mc.cores` is deliberately not set. Every brm() call passes `cores`
# explicitly, so a global would only supply a default nothing uses, and
# parallel::detectCores() returns NA on some machines.
options(
  scipen                              = 20,
  brms.backend                        = "cmdstanr",
  marginaleffects_posterior_interval  = "hdi",
  marginaleffects_posterior_center    = "median"
)


# ---- Simultaneous-band construction ----
#
# Read by `simul_credible_bands()` in `_startup/functions.R`, which documents
# both constructions:
#
#   HPDI       nested pointwise highest-density intervals, with their common
#              mass calibrated upward until whole curves are contained.
#              Asymmetric wherever the posterior is skewed. The default, and
#              what the figures in the paper use.
#   symmetric  mean +/- c* x pointwise SD, sup-norm calibrated, and therefore
#              symmetric about the mean whatever the posterior looks like.
#
# `BAND_TYPE=symmetric` in the environment rebuilds every simultaneous band the
# symmetric way without touching code. Validated here so a typo stops at session
# start rather than part-way through a figure build.
#
# Deliberately not in specifications.R: the targets pipeline hashes that file,
# so putting a plotting switch there would invalidate all 15 fits whenever it
# changed, and no fit depends on it.
local({
  band_type <- Sys.getenv("BAND_TYPE", "HPDI")
  if (!band_type %in% c("HPDI", "symmetric"))
    stop("BAND_TYPE must be 'HPDI' or 'symmetric', not '", band_type, "'.")
  options(bone.band_type = band_type)
})
