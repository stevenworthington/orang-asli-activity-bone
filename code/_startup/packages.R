###############################################################################
# Package loading. Uses the load_pkgs() wrapper defined in functions.R.
###############################################################################


load_pkgs(
  # Data
  "tidyverse", "janitor",
  # Bayesian modeling
  "brms", "cmdstanr", "rstan",
  # Posterior summaries
  "tidybayes", "posterior", "bayesplot", "bayestestR", "performance",
  # Marginal effects
  "marginaleffects",
  # Visualization
  "patchwork", "ggdist", "ggridges", "gghalves", "ggokabeito", "scales",
  # Utilities
  "splines", "mgcv", "HDInterval", "dagitty", "yaml"
)
