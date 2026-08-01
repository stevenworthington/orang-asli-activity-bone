# Data

**The data are not included in this repository. They are available upon request
from the corresponding author, Ian J. Wallace (iwallace@unm.edu).**

The analytical data are individual-level records from the Orang Asli, an
Indigenous population of Peninsular Malaysia. Access is governed by the Orang
Asli Health and Lifeways Project (OA HeLP) and the associated ethics approvals,
so the files are shared on request rather than distributed publicly.

`data/processed/` is kept as an empty directory so the pipeline has the location
it expects. Once you have obtained the data, place the two CSVs described below
into `data/processed/` and the pipeline (see the top-level `README.md`) will
reproduce every fitted model, figure, and table.

## Expected files

### `data/processed/Orang-Asli-pa-vs-bone-60126.csv`

The main analytical dataset (n = 1,007). Column names are snake-cased at load
(`janitor::clean_names`). The analysis code reads, among others:

| Column | Role |
|---|---|
| `age_years`, `sex`, `village_id` | age (tensor smooth), sex, community identity |
| `tibia_sos` | tibial speed-of-sound outcome (m/s; rescaled to `/1000` internally) |
| `ctx1_ng_ml` | CTX-1 resorption biomarker (ng/mL) |
| `osteocalcin_pg_ml` | osteocalcin formation biomarker (pg/mL; rescaled `/10000`) |
| `ad_tot_step_count_0_24hr` | daily step count exposure (rescaled `/1000`) |
| `ad_mean_enmo_mg_0_24hr` | mean daily ENMO exposure (mg) |
| `industrial_index` | community-level industrialization index (industrialization exposure) |
| `functional_status_n_y_0_1`, `pregnant_or_breastfeeding_n_y_0_1`, `smoking_binary_n_y_0_1`, `alcohol_binary_n_y_0_1` | adjustment-set covariates |
| `fat_mass_kg`, `fat_free_mass_kg` | body composition (z-scored; confounder-DAG sensitivity specs) |

`village_id` identifies the participant's community. The column keeps the data
file's original name, so the code, the variable names and the output paths all
say "village" where the manuscript says "community"; the two refer to the same
unit throughout.

This file excludes the 75 participants of community 43, whose bone-density
ultrasound was measured on an uncalibrated system.

### `data/processed/compiled-walking-experiment-18june2026.csv`

The acute walking-experiment dataset (n = 10; pre / t0 / t4 osteocalcin
measurements), used only by `code/_experiments/walking-experiment-rope.R` for
the supplementary reverse-ROPE analysis.
