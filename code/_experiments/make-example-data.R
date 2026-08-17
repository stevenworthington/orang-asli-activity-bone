###############################################################################
# Example-data generator — lets the pipeline run without the real data.
#
# The analytical data are not distributed with this repository (see
# data/README.md). That makes the code unrunnable for anyone who has not been
# granted access, so this script writes a SYNTHETIC dataset with the same
# columns, types and design shape, allowing the whole pipeline to be executed
# end to end as a working example.
#
# NOTHING HERE IS DERIVED FROM THE REAL DATA. Every value is drawn from a
# distribution written out in full below. No model was fitted to the real
# participants, nothing was resampled from them, and nothing was perturbed. The
# only quantities taken from the study are its design shape — the number of
# participants, the number of communities, and how many communities carry each
# measurement, all of which appear in the paper — together with the approximate
# scale of each variable, given below as round stand-ins rather than as the
# study's summary statistics. Reading this file is sufficient to confirm that.
#
# THE RESULTS IT PRODUCES ARE MEANINGLESS AS SCIENCE. The effects planted below
# are deliberately NOT the paper's: the physical-activity effect on bone is set
# clearly positive here, where the study reports a null. Anything computed from
# this dataset describes this generator, not the Orang Asli.
#
# Output: data/processed/orang-asli-pa-bone-data.csv  (gitignored)
#         data/processed/walking-experiment-data.csv
#
# IT WILL NOT OVERWRITE EITHER FILE. Those are the same paths the real data are
# placed at, so a reader who has been granted access and later runs the example
# would otherwise silently lose them. Set EXAMPLE_DATA_FORCE=1 to overwrite
# deliberately, or delete the files first.
#
# Run: Rscript code/_experiments/make-example-data.R
###############################################################################

FORCE <- nzchar(Sys.getenv("EXAMPLE_DATA_FORCE"))

# Synthetic files identify themselves, via a marker listing what this script
# wrote. That is what makes the guard usable: re-running the example overwrites
# its own output freely, while a file this script did not write -- the real data,
# placed at exactly these paths -- stops it. Deleting the marker makes the
# generator treat everything beside it as real, which is the safe direction.
#
# Refuse before generating, not after writing: a warning printed underneath a
# destroyed file is not a safeguard.
MARKER <- file.path("data", "processed", ".synthetic-example")
synthetic_already <- function(path)
  file.exists(MARKER) && basename(path) %in% readLines(MARKER, warn = FALSE)

WRITTEN <- character(0)

write_guarded <- function(x, path, what) {
  if (file.exists(path) && !synthetic_already(path) && !FORCE)
    stop("refusing to overwrite ", path, "\n",
         "  That path already holds a file this script did not write. If it is the\n",
         "  real ", what, ", writing synthetic data over it would destroy it.\n",
         "  Move it aside, or set EXAMPLE_DATA_FORCE=1 to overwrite deliberately.",
         call. = FALSE)
  write.csv(x, path, row.names = FALSE)
  WRITTEN <<- c(WRITTEN, basename(path))
}

set.seed(20260817)


# ---- Design shape (published study facts; no participant values involved) ----

N            <- 1007   # participants
N_COMMUNITY  <- 31     # communities represented in the file
N_SOS_COMM   <- 24     # communities with tibial SOS in the analytic sample
N_BIO_COMM   <- 9      # communities with biomarkers in the analytic sample
N_ACCEL      <- 867    # participants with accelerometry
N_SOS        <- 813    # participants with a tibial SOS measurement
N_BIO        <- 289    # participants with biomarker assays


# ---- Variable scales (rounded, of the order the paper reports) ---------------
#
# Deliberately round numbers rather than the published summary statistics, so
# that nothing here can be mistaken for a reconstruction of the real sample.

SOS_MEAN   <- 3840;  SOS_SD   <- 180     # m/s
STEPS_MEAN <- 11000; STEPS_SD <- 4500    # steps/day
ENMO_MEAN  <- 33;    ENMO_SD  <- 10      # mg
OSTEO_MED  <- 29000; OSTEO_LSD <- 0.55   # pg/mL, lognormal
CTX_MED    <- 0.13;  CTX_LSD  <- 0.60    # ng/mL, lognormal
INDEX_LO   <- 5;     INDEX_HI <- 42      # industrialization index

# Planted effects. NOT the study's findings -- see the header.
BETA_PA_SOS    <-  40      # m/s per 5,000 steps: clearly positive (study: null)
BETA_INDEX_STEP <- -90     # steps per index unit: a decline
BETA_AGE_SOS   <- -4.5     # m/s per year, post-peak


# ---- Communities -------------------------------------------------------------

# Unequal sizes, because the community-level analyses are sensitive to how
# lopsided the clusters are. Drawn lognormal and clipped to a plausible spread,
# then forced to sum to N.
sizes <- round(exp(rnorm(N_COMMUNITY, log(30), 0.55)))
sizes <- pmax(6L, pmin(80L, sizes))
sizes[N_COMMUNITY] <- sizes[N_COMMUNITY] + (N - sum(sizes))
while (sizes[N_COMMUNITY] < 6L) {                 # nudge if the fix-up went low
  j <- which.max(sizes); sizes[j] <- sizes[j] - 1L
  sizes[N_COMMUNITY] <- sizes[N_COMMUNITY] + 1L
}
stopifnot(sum(sizes) == N)

village_id <- rep(seq_len(N_COMMUNITY), times = sizes)
index_by_c <- round(seq(INDEX_LO, INDEX_HI, length.out = N_COMMUNITY) +
                      rnorm(N_COMMUNITY, 0, 1.2), 2)
industrial_index <- index_by_c[village_id]        # constant within community


# ---- Participants ------------------------------------------------------------

age_years <- pmin(91, pmax(18, round(18 + rgamma(N, shape = 2.2, scale = 11), 1)))
sex       <- sample(c("female", "male"), N, replace = TRUE, prob = c(0.62, 0.38))

functional_status_n_y_0_1        <- rbinom(N, 1, 0.875)
smoking_binary_n_y_0_1           <- rbinom(N, 1, 0.30)
alcohol_binary_n_y_0_1           <- rbinom(N, 1, 0.18)
pregnant_or_breastfeeding_n_y_0_1 <- ifelse(sex == "female" & age_years < 45,
                                            rbinom(N, 1, 0.22), 0L)

fat_mass_kg      <- pmax(2,  round(rnorm(N, 15, 6), 1))
fat_free_mass_kg <- pmax(20, round(rnorm(N, 38, 7), 1))


# ---- Activity: declines with the index, and with age -------------------------

comm_step_offset <- rnorm(N_COMMUNITY, 0, 1200)[village_id]
ad_tot_step_count_0_24hr <- pmax(200, round(
  STEPS_MEAN + BETA_INDEX_STEP * (industrial_index - mean(index_by_c)) +
    comm_step_offset - 60 * (age_years - 40) + rnorm(N, 0, STEPS_SD * 0.6)))

comm_enmo_offset <- rnorm(N_COMMUNITY, 0, 2.5)[village_id]
ad_mean_enmo_mg_0_24hr <- pmax(3, round(
  ENMO_MEAN + (BETA_INDEX_STEP / 4500) * 10 * (industrial_index - mean(index_by_c)) +
    comm_enmo_offset - 0.15 * (age_years - 40) + rnorm(N, 0, ENMO_SD * 0.6), 2))


# ---- Outcomes ----------------------------------------------------------------

steps_5k  <- ad_tot_step_count_0_24hr / 5000
comm_sos_offset <- rnorm(N_COMMUNITY, 0, 45)[village_id]
tibia_sos <- round(
  SOS_MEAN + BETA_PA_SOS * (steps_5k - mean(steps_5k)) +
    BETA_AGE_SOS * pmax(0, age_years - 35) + comm_sos_offset +
    rt(N, df = 6) * SOS_SD * 0.55, 1)

osteocalcin_pg_ml <- round(exp(log(OSTEO_MED) - 0.004 * (age_years - 40) +
                                 rnorm(N, 0, OSTEO_LSD)))
ctx1_ng_ml        <- round(exp(log(CTX_MED) - 0.003 * (age_years - 40) +
                                 rnorm(N, 0, CTX_LSD)), 4)


# ---- Missingness, set so the analytic samples match the real design ----------
#
# Tibial SOS is measured in N_SOS_COMM communities and biomarkers in a smaller
# N_BIO_COMM subset of those, which is what gives the fitted models their
# community counts (and so the number of community indicators they carry).

sos_comms <- sort(sample(seq_len(N_COMMUNITY), N_SOS_COMM))
bio_comms <- sort(sample(sos_comms,            N_BIO_COMM))

drop_to <- function(eligible, target) {
  idx <- which(eligible)
  if (length(idx) <= target) return(!eligible)     # nothing to drop
  setdiff_idx <- sample(idx, length(idx) - target)
  out <- rep(FALSE, length(eligible)); out[setdiff_idx] <- TRUE; out
}

tibia_sos[drop_to(village_id %in% sos_comms, N_SOS)] <- NA
is_bio <- village_id %in% bio_comms
osteocalcin_pg_ml[drop_to(is_bio, N_BIO)] <- NA
ctx1_ng_ml[is.na(osteocalcin_pg_ml)]      <- NA     # assayed together
tibia_sos[!(village_id %in% sos_comms)]   <- NA

accel_missing <- drop_to(rep(TRUE, N), N_ACCEL)
ad_tot_step_count_0_24hr[accel_missing] <- NA
ad_mean_enmo_mg_0_24hr[accel_missing]   <- NA

# A little missingness on body composition, which is why the confounder-DAG
# specs fit on slightly fewer rows than their mediator-DAG counterparts.
bc_missing <- sample(c(TRUE, FALSE), N, replace = TRUE, prob = c(0.03, 0.97))
fat_mass_kg[bc_missing]      <- NA
fat_free_mass_kg[bc_missing] <- NA


# ---- Write -------------------------------------------------------------------

dat_example <- data.frame(
  age_years, sex, village_id, industrial_index,
  tibia_sos, ctx1_ng_ml, osteocalcin_pg_ml,
  ad_tot_step_count_0_24hr, ad_mean_enmo_mg_0_24hr,
  functional_status_n_y_0_1, pregnant_or_breastfeeding_n_y_0_1,
  smoking_binary_n_y_0_1, alcohol_binary_n_y_0_1,
  fat_mass_kg, fat_free_mass_kg
)

out_dir <- here::here("data", "processed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out_dir, "orang-asli-pa-bone-data.csv")
write_guarded(dat_example, out, "analytical dataset")
cat("wrote", out, sprintf("(%d rows, %d communities)\n",
                          nrow(dat_example), length(unique(village_id))))


# ---- The walking-experiment dataset, same treatment --------------------------

n_walk <- 10
walk <- data.frame(
  participant_id = rep(seq_len(n_walk), each = 3),
  time_point     = rep(c("pre", "t0", "t4"), times = n_walk)
)
walk_base <- rep(exp(rnorm(n_walk, log(22), 0.3)), each = 3)
walk$osteocalcin_ng_ml <- round(
  walk_base * exp(rep(c(0, 0.02, 0.01), times = n_walk) + rnorm(nrow(walk), 0, 0.08)), 2)

out_walk <- file.path(out_dir, "walking-experiment-data.csv")
write_guarded(walk, out_walk, "walking-experiment dataset")
cat("wrote", out_walk, sprintf("(%d participants x 3 timepoints)\n", n_walk))

# Written last, and only on success: the marker is what lets a re-run overwrite
# its own output while still refusing to touch anything else.
writeLines(WRITTEN, MARKER)

cat("\nThese files are SYNTHETIC. Results computed from them are not the study's.\n")
cat("Recorded in", MARKER, "-- delete it and the generator will treat these as real.\n")
