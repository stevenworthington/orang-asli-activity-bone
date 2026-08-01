# Outcome natural variability (mean / median / SD / CV) for the noise-anchor table.
# Base R only; reads the analytic CSV. Pairs with instrument-accuracy noise (assay
# CVs from Methods; qUS manufacturer accuracy) to give two reference scales per outcome.
dat <- read.csv("data/processed/Orang-Asli-pa-vs-bone-60126.csv")
cols <- c("Tibial SOS (m/s)"      = "tibia_sos",
          "Osteocalcin (pg/mL)"   = "osteocalcin_pg_ml",
          "CTX-1 (ng/mL)"         = "ctx1_ng_ml",
          "Daily steps"           = "ad_tot_step_count_0_24hr",
          "Mean daily ENMO (mg)"  = "ad_mean_enmo_mg_0_24hr")
for (nm in names(cols)) {
  cn <- cols[[nm]]
  if (!cn %in% names(dat)) { cat(sprintf("%-22s  COLUMN '%s' NOT FOUND\n", nm, cn)); next }
  v <- dat[[cn]]; v <- v[is.finite(v)]
  cat(sprintf("%-22s n=%4d | mean=%11.3f | median=%11.3f | SD=%11.3f | CV=%5.1f%%\n",
              nm, length(v), mean(v), median(v), sd(v), 100 * sd(v) / mean(v)))
}
