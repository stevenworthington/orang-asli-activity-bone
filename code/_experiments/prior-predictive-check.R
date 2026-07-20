###############################################################################
# Prior predictive checks for the Bayesian GAM specs.
#
# Restores the prior-predictive workflow that lived in the pre-narrowing matrix
# (.graveyard/bayes-osteo-mins.Rmd § "Prior Predictive Check") and re-points it
# at the current spec registry, so the manuscript's prior-predictive-validation
# claim (Methods / Supplementary Material) has a live, re-runnable artifact.
#
# For each representative spec it refits the SAME bf()/family/priors as the
# production pipeline but with sample_prior = "only" (no likelihood), then draws
# the prior predictive distribution of the outcome mean. A sensible prior puts
# the observed outcome mean comfortably inside the bulk of that distribution
# without being so tight that it drives the posterior -- the student_t(3, 0, 2.5)
# scales are deliberately broad on the rescaled outcomes (see specifications.R).
#
# Run (pty wrapper, sandbox off -- see CODING.md "Running R from an agent"):
#   script -q /dev/null Rscript code/_experiments/prior-predictive-check.R < /dev/null
###############################################################################


library(here)
source(here("code", "_startup", "init.R"))


# ---- Representative specs: every outcome family x prior set ----

ppc_keys <- c(
  "sos-steps",   # Gaussian outcome (SOS),   mediator priors
  "ctx-steps",   # lognormal outcome (CTX),  mediator priors
  "osteo-steps", # lognormal outcome (osteo),mediator priors
  "sos-urb",     # Gaussian outcome (SOS),   age+sex priors
  "steps-urb",   # Gaussian outcome (steps), age+sex priors
  "enmo-urb"     # Gaussian outcome (ENMO),  age+sex priors
)

# Prior-only sampling needs far fewer draws than posterior inference -- the goal
# is to visualize the prior predictive spread, not to characterize a posterior.
PPC_WARMUP <- 500
PPC_ITER   <- 1500
PPC_CHAINS <- 4


# ---- Prior-only refit for one spec ----

fit_prior_only <- function(spec_key, dat_raw) {
  spec      <- model_templates[[spec_key]]
  dat_local <- prep_local_data(spec, dat_raw)
  # sample_prior = "only" needs a PROPER prior on every parameter. The production
  # fit leaves the tensor-smooth's fixed-effect coefficients (class "b" on the
  # t2() terms) flat/improper -- fine for posterior sampling, but it makes the
  # prior alone improper and unsamplable. We add the project-standard weak prior
  # student_t(3, 0, 2.5) on those remaining b coefficients so the prior predictive
  # is well-defined; the named covariate coefficients keep their registry priors.
  ppc_prior <- c(spec$priors, set_prior("student_t(3, 0, 2.5)", class = "b"))
  brm(
    spec$bf,
    data         = dat_local,
    prior        = ppc_prior,
    sample_prior = "only",
    warmup       = PPC_WARMUP,
    iter         = PPC_ITER,
    thin         = 1,
    chains       = PPC_CHAINS,
    cores        = PPC_CHAINS,
    seed         = SEED,
    backend      = "cmdstanr",
    control      = BRMS_CONTROL,
    refresh      = 0,
    silent       = 2
  )
}


# ---- Run + render one prior-predictive mean check per spec ----

out_dir <- here("outputs", "figures", "working", "prior-predictive")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_pdf <- file.path(out_dir, "prior-predictive-checks.pdf")

pdf(out_pdf, width = 9, height = 6)
for (k in ppc_keys) {
  spec <- model_templates[[k]]
  message("prior-predictive: ", k)
  fit  <- fit_prior_only(k, dat)
  p <- pp_check(fit, type = "stat", stat = "mean", ndraws = 500) +
    ggtitle(sprintf("%s  (%s ~ %s)", k, spec$outcome, spec$exposure))
  print(p)
}
invisible(dev.off())
cat("saved", out_pdf, "\n")
