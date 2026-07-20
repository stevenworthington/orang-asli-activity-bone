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
# within-village (village fixed-effects) identification, and tests the exposure
# coefficient directly. Outcomes on their modeling scale: tibial SOS native;
# osteocalcin and CTX-1 log-transformed (the lognormal twin).  detect = exposure
# coefficient's 95% CI excludes 0.
#
# DGP: a baseline fit WITHOUT exposure supplies the age/sex/covariate/village
# structure and residual SD; a true per-contrast effect of size X (in SD units,
# over a 5,000-step or 10-mg ENMO contrast) is added as a linear exposure term, and
# Gaussian residuals are drawn. Real design throughout (real covariates, villages,
# exposure distribution).
#
# Outputs: outputs/_experiments/power-curves/pa-bone-power.rds + .csv
# Run: script -q /dev/null Rscript code/_experiments/power-curves-pa-bone.R < /dev/null
# Smoke: POW_SMOKE=1
###############################################################################

suppressMessages({ library(here); library(splines) })
set.seed(2138)

out_dir <- here("outputs", "_experiments", "power-curves")
if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE)
  try(system2("xattr", c("-w", "com.dropbox.ignored", "1", out_dir)), silent = TRUE) }

SMOKE <- nzchar(Sys.getenv("POW_SMOKE"))
N_SIM <- if (SMOKE) 40 else 2000
NGRID <- if (SMOKE) 4  else 11

d <- read.csv(here("data", "processed", "Orang-Asli-pa-vs-bone-60126.csv"))

COVS  <- c("male", "pregnant_or_breastfeeding_n_y_0_1", "smoking_binary_n_y_0_1",
           "alcohol_binary_n_y_0_1", "functional_status_n_y_0_1")
STEPS <- "ad_tot_step_count_0_24hr"; ENMO <- "ad_mean_enmo_mg_0_24hr"

# analysis registry: outcome col, log?, exposure col, contrast, resolution floor (modeling-scale units)
analyses <- list(
  list(key = "sos-steps",  y = "tibia_sos",         logy = FALSE, x = STEPS, contrast = 5000, res = 32),
  list(key = "sos-enmo",   y = "tibia_sos",         logy = FALSE, x = ENMO,  contrast = 10,   res = 32),
  list(key = "ctx-steps",  y = "ctx1_ng_ml",        logy = TRUE,  x = STEPS, contrast = 5000, res = 0.043), # interassay CV
  list(key = "ctx-enmo",   y = "ctx1_ng_ml",        logy = TRUE,  x = ENMO,  contrast = 10,   res = 0.043),
  list(key = "osteo-steps",y = "osteocalcin_pg_ml", logy = TRUE,  x = STEPS, contrast = 5000, res = 0.028),
  list(key = "osteo-enmo", y = "osteocalcin_pg_ml", logy = TRUE,  x = ENMO,  contrast = 10,   res = 0.028)
)

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
  base_f <- as.formula(paste0("Yobs ~ ns(age_c, 4) + ", covterms, " + village_id"))
  full_f <- as.formula(paste0("ysim ~ ns(age_c, 4) + expo + ", covterms, " + village_id"))
  base   <- lm(base_f, data = s)
  nuis   <- fitted(base); sig <- summary(base)$sigma

  # grid of per-contrast effects in SD units (informative range: landmarks sit at
  # ~0.05-0.2 SD, and power saturates by ~0.5 SD); map to per-unit slope
  res_sd   <- a$res / sd_y                                   # resolution floor in SD units
  # include the exact landmark effect sizes (resolution floor + 0.2 SD) as grid
  # points so landmark power is a direct simulated readout, not interpolated
  x_sd     <- sort(unique(c(seq(0, 0.5, length.out = NGRID), res_sd, 0.2)))
  detect1  <- function(xsd) {
    slope <- (xsd * sd_y) / a$contrast                       # per-unit-exposure slope
    s$ysim <- nuis + slope * s$expo + rnorm(nrow(s), 0, sig)
    fit <- lm(full_f, data = s)
    ci  <- suppressWarnings(confint(fit, "expo"))
    as.numeric(ci[1] > 0 | ci[2] < 0)
  }
  pw <- vapply(x_sd, function(xsd) mean(replicate(N_SIM, detect1(xsd)), na.rm = TRUE), numeric(1))
  cat(sprintf("\n%s  (n=%d, %d villages; modeling-scale SD=%.3g; resolution=%.2f SD)\n",
              a$key, nrow(s), nlevels(s$village_id), sd_y, res_sd))
  for (i in seq_along(x_sd)) cat(sprintf("  effect=%.2f SD  power=%.3f\n", x_sd[i], pw[i]))
  data.frame(key = a$key, outcome = a$y, exposure = a$x, effect_sd = x_sd, power = pw,
             resolution_sd = res_sd, n = nrow(s), n_villages = nlevels(s$village_id))
}

if (SMOKE) analyses <- analyses[c(1, 3)]
res <- do.call(rbind, lapply(analyses, power_curve))
saveRDS(res, file.path(out_dir, "pa-bone-power.rds"))
write.csv(res, file.path(out_dir, "pa-bone-power.csv"), row.names = FALSE)
cat("\nsaved", file.path(out_dir, "pa-bone-power.csv"), "\n")
