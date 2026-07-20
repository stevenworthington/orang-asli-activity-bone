# Age -> outcome rates for the biological "age-equivalent" anchor.
# NON-PARAMETRIC: fit a GAM smooth s(age) + sex (mgcv) and read the rate off the
# fitted curve -- no global linearity assumption. Report (a) the average post-peak
# rate as the SECANT of the smooth from age 35 (peak bone mass) to 70, (b) a LOCAL
# slope of the smooth at age 55, and (c) the linear lm(age>=35) slope for comparison.
suppressMessages(library(mgcv))
dat <- read.csv("data/processed/Orang-Asli-pa-vs-bone-60126.csv")
outcomes <- c("Tibial SOS (m/s)"     = "tibia_sos",
              "Osteocalcin (pg/mL)"  = "osteocalcin_pg_ml",
              "CTX-1 (ng/mL)"        = "ctx1_ng_ml",
              "Daily steps"          = "ad_tot_step_count_0_24hr",
              "Mean daily ENMO (mg)" = "ad_mean_enmo_mg_0_24hr")
rs <- names(which.max(table(dat$sex)))   # reference sex (additive, so secant is sex-invariant)

for (nm in names(outcomes)) {
  cn <- outcomes[[nm]]
  d  <- dat[is.finite(dat[[cn]]) & is.finite(dat$age_years) & !is.na(dat$sex), ]
  d$y <- d[[cn]]
  m  <- gam(y ~ s(age_years) + sex, data = d)
  pr <- function(a) as.numeric(predict(m, data.frame(age_years = a, sex = rs)))
  sec <- (pr(70) - pr(35)) / (70 - 35)          # avg rate over post-peak span (non-parametric)
  loc <- (pr(56) - pr(54)) / 2                  # local slope of the smooth at ~55
  tot <- pr(70) - pr(35)                        # total smoothed change 35 -> 70
  lin <- unname(coef(lm(y ~ age_years + sex, data = d[d$age_years >= 35, ]))["age_years"])
  cat(sprintf("%-22s | smooth secant 35->70: %10.3f /yr | local @55: %10.3f /yr | linear(>=35): %10.3f /yr | total 35->70: %10.1f\n",
              nm, sec, loc, lin, tot))
}
