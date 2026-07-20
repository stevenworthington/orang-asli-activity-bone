###############################################################################
# Global options for the session. Paper-shape defaults for brms /
# marginaleffects on top of base-R options.
###############################################################################


options(
  scipen                              = 20,
  mc.cores                            = parallel::detectCores(),
  brms.backend                        = "cmdstanr",
  marginaleffects_posterior_interval  = "hdi",
  marginaleffects_posterior_center    = "median"
)
