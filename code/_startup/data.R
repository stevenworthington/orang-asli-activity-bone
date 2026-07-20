###############################################################################
# Read + prep the project's working dataset (`dat`). The OA HeLP physical-
# activity / bone CSV from Ian; the 2026-06-01 revision (n = 1,007, unique
# tid/rid/vid, single canonical `industrial_index` column — supersedes the
# prior two-file `dat` / `dat_full` topology) drops the 75 Sarok (village_id
# 43) participants, excluded because the bone-density ultrasound system was
# not properly calibrated there (new operator). Derives rescaled columns so
# the shared priors on Intercept / sigma are on the right scale: tibial SOS
# → /1000, osteocalcin → /10000, daily step count → /1000. CTX-1 and ENMO
# are already on small scales and left raw. Factors `sex` and `village_id`.
###############################################################################


read_csv(here("data", "processed", "Orang-Asli-pa-vs-bone-60126.csv"),
         show_col_types = FALSE) |>
  janitor::clean_names(case = "snake") |>
  mutate(
    sex                   = factor(sex, levels = c("female", "male")),
    village_id            = factor(village_id),
    # Outcome rescaling (priors are calibrated to these scaled outcomes;
    # `outcome_scale_factor` in specifications.R reverses the rescaling for
    # downstream pred-draws / slope-draws / figures).
    tibia_sos_1k          = tibia_sos / 1000,
    osteocalcin_pg_ml_10k = osteocalcin_pg_ml / 10000,
    # Exposure rescaling (canonical k-step pattern).
    ad_steps_1k           = ad_tot_step_count_0_24hr / 1000,
    # Body composition (z-scored). Used by the confounder-DAG variants of the
    # PA -> bone analyses; the canonical 6-var MSAS includes Fat mass & lean
    # body mass, operationalized as two separate regression covariates per the
    # 2026-05-19 review (issue 1). Z-scoring keeps the b ~ student_t(3, 0, 2.5)
    # prior interpretable as "per-SD effect on the outcome scale."
    fat_mass_kg_z         = as.numeric(scale(fat_mass_kg)),
    fat_free_mass_kg_z    = as.numeric(scale(fat_free_mass_kg))
  ) ->
dat
