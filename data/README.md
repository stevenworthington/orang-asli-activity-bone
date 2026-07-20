# Data

**No data files are included in this snapshot.** Place the two analytical CSVs
below into `data/processed/` to run the pipeline. They are individual-level
records from the Orang Asli cohort; obtain them per the manuscript's
data-availability statement and the OA HeLP access terms (see `../LICENSE-data.md`).

## Expected files

### `data/processed/Orang-Asli-pa-vs-bone-60126.csv`

The main analytical dataset (n = 1,007). Column names are snake-cased at load
(`janitor::clean_names`). The analysis code reads, among others:

| Column | Role |
|---|---|
| `age_years`, `sex`, `village_id` | age (tensor smooth), sex, village identity |
| `tibia_sos` | tibial speed-of-sound outcome (m/s; rescaled to `/1000` internally) |
| `ctx1_ng_ml` | CTX-1 resorption biomarker (ng/mL) |
| `osteocalcin_pg_ml` | osteocalcin formation biomarker (pg/mL; rescaled `/10000`) |
| `ad_tot_step_count_0_24hr` | daily step count exposure (rescaled `/1000`) |
| `ad_mean_enmo_mg_0_24hr` | mean daily ENMO exposure (mg) |
| `industrial_index` | community-level industrialization index (Set 2 exposure) |
| `functional_status_n_y_0_1`, `pregnant_or_breastfeeding_n_y_0_1`, `smoking_binary_n_y_0_1`, `alcohol_binary_n_y_0_1` | adjustment-set covariates |
| `fat_mass_kg`, `fat_free_mass_kg` | body composition (z-scored; confounder-DAG sensitivity specs) |

The 2026-06-01 revision excludes the 75 Sarok (village 43) participants, whose
bone-density ultrasound was flagged as uncalibrated.

### `data/processed/compiled-walking-experiment-18june2026.csv`

The acute walking-experiment dataset (n = 10; pre / t0 / t4 osteocalcin
measurements), used only by `code/_experiments/walking-experiment-rope.R` for
the Supplementary Material 11 reverse-ROPE analysis.
