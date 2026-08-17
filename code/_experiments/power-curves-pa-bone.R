###############################################################################
# A PRIORI power curves for the six physical-activity -> bone analyses, using the
# FREQUENTIST PROXY. A priori, NOT post-hoc: power at pre-specified true effect
# sizes swept across a grid (in units of the outcome's between-person SD), with the
# resolution-floor and 0.2-SD landmarks marked. N_SIM replicates per grid point.
#
# Estimator (frequentist twin of the reported Bayesian GAM): the reported method
# fits a Bayesian GAM with t2(age, exposure) and reports the AERF's linear-projection
# slope. For the power of that slope we use the fast linear twin --
#   lm(Yscale ~ ns(age_c, 4) + exposure + male + preg + smoke + alcohol + functional
#               + village_id)
# -- which carries flexible age, the full mediator-DAG adjustment set, and the
# within-community (community fixed-effects) identification, and tests the
# exposure coefficient directly. Each twin matches its spec's finalized
# likelihood, so the simulated data carries the tails the likelihood was chosen
# for: tibial SOS via mgcv::scat (Student-t, nu read off the data); CTX-1 via a
# Gamma(log) GLM on the RAW response; osteocalcin log-transformed (the lognormal
# twin). detect = exposure coefficient's 95% Wald CI excludes 0, for all six.
#
# DGP: a baseline fit WITHOUT exposure supplies the age/sex/covariate/community
# structure and residual SD; a true per-contrast effect of size X (in SD units,
# over a 5,000-step or 10-mg ENMO contrast) is added as a linear exposure term, and
# Gaussian residuals are drawn. Real design throughout (real covariates, communities,
# exposure distribution).
#
# Outputs: outputs/_experiments/power-curves/pa-bone-power.rds + .csv
# Run: Rscript code/_experiments/power-curves-pa-bone.R
# Smoke: POW_SMOKE=1
###############################################################################

suppressMessages({ library(here); library(splines); library(mgcv) })
set.seed(2138)

out_dir <- here("outputs", "_experiments", "power-curves")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

SMOKE <- nzchar(Sys.getenv("POW_SMOKE"))
N_SIM <- if (SMOKE) 40 else 2000
NGRID <- if (SMOKE) 4  else 11

d <- read.csv(here("data", "processed", "orang-asli-pa-bone-data.csv"))

COVS  <- c("male", "pregnant_or_breastfeeding_n_y_0_1", "smoking_binary_n_y_0_1",
           "alcohol_binary_n_y_0_1", "functional_status_n_y_0_1")
STEPS <- "ad_tot_step_count_0_24hr"; ENMO <- "ad_mean_enmo_mg_0_24hr"

# analysis registry. `fam` is the production likelihood, and it selects the
# frequentist twin used below -- one twin per family, because a Gaussian lm is
# the right twin only for the two lognormal arms.
#
#   lognormal -> lm on log(y)                  exactly the twin
#   gamma     -> glm(family = Gamma(log))      CTX-1
#   student   -> gam(family = scat())          tibial SOS, nu read off the data
#
# `logy` means only "measure effect sizes in SD units of log(y)", which is the
# right scale whenever the linear predictor is a log scale -- so it is TRUE for
# Gamma as well as lognormal.
analyses <- list(
  list(key = "sos-steps",  y = "tibia_sos",         fam = "student",   logy = FALSE, x = STEPS, contrast = 5000, res = 32),
  list(key = "sos-enmo",   y = "tibia_sos",         fam = "student",   logy = FALSE, x = ENMO,  contrast = 10,   res = 32),
  list(key = "ctx-steps",  y = "ctx1_ng_ml",        fam = "gamma",     logy = TRUE,  x = STEPS, contrast = 5000, res = 0.043), # interassay CV
  list(key = "ctx-enmo",   y = "ctx1_ng_ml",        fam = "gamma",     logy = TRUE,  x = ENMO,  contrast = 10,   res = 0.043),
  list(key = "osteo-steps",y = "osteocalcin_pg_ml", fam = "lognormal", logy = TRUE,  x = STEPS, contrast = 5000, res = 0.028),
  list(key = "osteo-enmo", y = "osteocalcin_pg_ml", fam = "lognormal", logy = TRUE,  x = ENMO,  contrast = 10,   res = 0.028)
)

# Wald interval on the exposure coefficient, used for EVERY family so the six
# arms are decided the same way. At n = 265-650 it differs from a t-based
# confint() by < 1%, and it is the one rule that all three twins support.
wald_excludes_zero <- function(est, se) {
  if (!is.finite(est) || !is.finite(se) || se <= 0) return(NA_real_)
  as.numeric(abs(est) > stats::qnorm(0.975) * se)
}

power_curve <- function(a) {
  need <- c(a$y, a$x, "age_years", "sex", "village_id",
            "pregnant_or_breastfeeding_n_y_0_1", "smoking_binary_n_y_0_1",
            "alcohol_binary_n_y_0_1", "functional_status_n_y_0_1")
  s <- d[stats::complete.cases(d[, need]) & is.finite(d[[a$y]]) & is.finite(d[[a$x]]), ]
  s$village_id <- droplevels(factor(s$village_id))
  s$age_c <- as.numeric(scale(s$age_years))
  s$male  <- as.integer(s$sex == "male")
  s$expo  <- s[[a$x]]
  s$Yobs  <- if (a$logy) log(s[[a$y]]) else s[[a$y]]
  sd_y    <- sd(s$Yobs)

  covterms <- paste(COVS, collapse = " + ")
  # Gamma models the response on its own scale, so its base fit takes RAW y --
  # `Yobs` is log(y) for the log-scale arms and would be negative for CTX-1,
  # which Gamma rejects. `sd_y` above stays sd(log(y)) either way: that is the
  # scale the effect grid is expressed in, and it is a log scale under Gamma too.
  s$Yraw   <- s[[a$y]]
  resp     <- if (a$fam == "gamma") "Yraw" else "Yobs"
  base_f <- as.formula(paste0(resp, " ~ ns(age_c, 4) + ", covterms, " + village_id"))
  full_f <- as.formula(paste0("ysim ~ ns(age_c, 4) + expo + ", covterms, " + village_id"))
  # Nuisance structure and noise model, both taken from a fit of the SAME family
  # the production model uses, so the simulated data has the tail behaviour the
  # likelihood was chosen for. Simulating Gaussian noise and fitting a Student-t
  # twin would drive nu to infinity, collapsing it back to a Gaussian and
  # measuring the power of a model this study does not use.
  if (a$fam == "student") {
    base <- mgcv::gam(base_f, data = s, family = mgcv::scat(), method = "REML")
    th   <- base$family$getTheta(TRUE)            # c(nu, sigma) on the natural scale
    NU   <- as.numeric(th[1]); SG <- as.numeric(th[2])
    nuis <- as.numeric(fitted(base)); sig <- SG
  } else if (a$fam == "gamma") {
    base  <- stats::glm(base_f, data = s, family = stats::Gamma(link = "log"))
    nuis  <- as.numeric(stats::predict(base, type = "link"))   # eta, log scale
    SHAPE <- 1 / summary(base)$dispersion
    sig   <- NA_real_
  } else {
    base <- lm(base_f, data = s)
    nuis <- fitted(base); sig <- summary(base)$sigma
  }

  # grid of per-contrast effects in SD units (informative range: landmarks sit at
  # ~0.05-0.2 SD, and power saturates by ~0.5 SD); map to per-unit slope
  res_sd   <- a$res / sd_y                                   # resolution floor in SD units
  # include the exact landmark effect sizes (resolution floor + 0.2 SD) as grid
  # points so landmark power is a direct simulated readout, not interpolated
  x_sd     <- sort(unique(c(seq(0, 0.5, length.out = NGRID), res_sd, 0.2)))
  detect1  <- function(xsd) {
    slope <- (xsd * sd_y) / a$contrast                       # per-unit-exposure slope
    if (a$fam == "student") {
      # scat residuals: sigma * t_nu, i.e. the fitted model's own noise law
      s$ysim <- nuis + slope * s$expo + sig * stats::rt(nrow(s), df = NU)
      gm <- tryCatch(mgcv::gam(full_f, data = s, family = mgcv::scat(), method = "REML"),
                     error = function(e) NULL)
      if (is.null(gm)) return(NA_real_)
      pt <- summary(gm)$p.table
      if (!("expo" %in% rownames(pt))) return(NA_real_)
      wald_excludes_zero(pt["expo", 1], pt["expo", 2])
    } else if (a$fam == "gamma") {
      mu <- exp(nuis + slope * s$expo)
      s$ysim <- stats::rgamma(nrow(s), shape = SHAPE, rate = SHAPE / mu)
      gl <- tryCatch(stats::glm(full_f, data = s, family = stats::Gamma(link = "log")),
                     error = function(e) NULL)
      if (is.null(gl)) return(NA_real_)
      co <- summary(gl)$coefficients
      if (!("expo" %in% rownames(co))) return(NA_real_)
      wald_excludes_zero(co["expo", 1], co["expo", 2])
    } else {
      s$ysim <- nuis + slope * s$expo + rnorm(nrow(s), 0, sig)
      fit <- lm(full_f, data = s)
      co  <- summary(fit)$coefficients
      wald_excludes_zero(co["expo", 1], co["expo", 2])
    }
  }
  pw <- vapply(x_sd, function(xsd) mean(replicate(N_SIM, detect1(xsd)), na.rm = TRUE), numeric(1))
  cat(sprintf("\n%s [%s]  (n=%d, %d communities; modeling-scale SD=%.3g; resolution=%.2f SD)\n",
              a$key, a$fam, nrow(s), nlevels(s$village_id), sd_y, res_sd))
  for (i in seq_along(x_sd)) cat(sprintf("  effect=%.2f SD  power=%.3f\n", x_sd[i], pw[i]))
  data.frame(key = a$key, family = a$fam, outcome = a$y, exposure = a$x,
             effect_sd = x_sd, power = pw,
             resolution_sd = res_sd, n = nrow(s), n_villages = nlevels(s$village_id))
}

if (SMOKE) analyses <- analyses[c(1, 3, 5)]   # one arm per family
res <- do.call(rbind, lapply(analyses, power_curve))
saveRDS(res, file.path(out_dir, "pa-bone-power.rds"))
write.csv(res, file.path(out_dir, "pa-bone-power.csv"), row.names = FALSE)
cat("\nsaved", file.path(out_dir, "pa-bone-power.csv"), "\n")
