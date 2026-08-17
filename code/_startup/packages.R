###############################################################################
# Package loading. Uses the load_pkgs() wrapper defined in functions.R.
#
# Attaching is for interactive convenience and for the operators and S3 methods
# that need it (patchwork's `+` and `/`, ggplot2's `+`). Every call in the
# codebase is namespaced, so nothing here is load-bearing for correctness — but
# what is attached decides what masks what, which the check in
# `code/_checks/verify-startup.R` enforces.
#
# The Stan backend is cmdstanr, set in options.R.
###############################################################################


load_pkgs(
  # Data
  "tidyverse", "janitor",
  # Bayesian modeling
  "brms", "cmdstanr",
  # Posterior summaries
  "tidybayes", "posterior", "bayesplot", "bayestestR", "performance",
  # Marginal effects
  "marginaleffects",
  # Visualization
  "patchwork", "ggdist", "ggridges", "gghalves", "ggokabeito", "scales",
  # Utilities
  "splines", "mgcv", "HDInterval", "dagitty", "yaml"
)
