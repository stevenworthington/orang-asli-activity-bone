# Replication code for "Industrialization, physical activity, and skeletal health"

Replication code for:

> Sena CMT, Worthington S, Kraus VB, Huebner JL, Lim YAL, Tan BT, Tam KL, Heng ESZ, Kivell TL, Holowka NB, Fadzil KS, Nicholas C, Sayed IM, Christopher KL, Watowich MM, Venkataraman VV, Lea AJ, Kraft TS, Wallace IJ. *Industrialization, physical activity, and skeletal health: Evidence from a within-population study of the Indigenous peoples of Peninsular Malaysia.* (under review, *Journal of Bone and Mineral Research*).

The study is a causal analysis of physical activity and bone health in the Orang Asli Health and Lifeways Project (OA HeLP) cohort (n = 1,007), across two analysis sets:

- **Set 1 — physical activity → bone (mediator DAG, primary).** Daily step count and mean daily ENMO against tibial speed-of-sound (SOS), CTX-1, and osteocalcin; within-village (village fixed-effects) estimand. Adjustment: age (in the tensor smooth) + sex + functional status + pregnancy/lactation + smoking + alcohol.
- **Set 2 — industrialization → outcome.** A community-level industrialization index against tibial SOS, daily steps, and mean daily ENMO; across-community estimand via a two-stage cluster-honest village-level estimator. Adjustment: age + sex.

Estimands are average exposure–response functions (AERFs), average marginal-effect functions (AMEFs), and their curvature, reported as linear-projection slopes with 95% simultaneous highest-posterior-density intervals.

## Data availability

**The analytical data are not included in this repository snapshot.** `data/processed/` is present but empty except for a `README.md` describing the two CSVs the code expects. See `data/README.md` and the manuscript's data-availability statement for how to obtain them. The data are individual-level records from a vulnerable Indigenous cohort; access terms are set by the OA HeLP team.

Once the data files are in place under `data/processed/`, the pipeline below reproduces every fitted model, figure, and table.

## Software environment

- **R 4.6.0** with the exact package versions pinned in `renv.lock` (329 packages).
- **CmdStan** via `{cmdstanr}` is the Stan backend for `{brms}`. Install it after restoring the library (`cmdstanr::install_cmdstan()`); the manuscript used Stan 2.38.0.

Restore the environment from the project root:

```r
# from R, in the repository root
renv::restore()          # rebuilds the pinned library from renv.lock
cmdstanr::install_cmdstan()   # one-time: installs the Stan toolchain
```

The model fits are computationally heavy (each brms spec is 10 chains × several thousand iterations). Budget accordingly and run them on a machine with adequate cores/RAM.

## Reproduce

Run order, from the repository root. `just` recipes are the front door (`just` with no arguments lists them); the underlying `Rscript` calls are shown for reference.

1. **Restore the environment** — `renv::restore()` (see above).
2. **Place the data** — copy the two CSVs into `data/processed/` (see `data/README.md`).
3. **Fit the models** — `just fit-all` runs the `_targets.R` pipeline: 15 specifications (6 mediator-DAG + 6 confounder-DAG + 3 industrialization), each producing a fitted `brms` object plus AERF / AMEF / curvature posterior draws. Only out-of-date targets re-run. Individual specs can also be run directly (`just fit-script sos-steps`, or `Rscript code/sos-steps.R`); a fast smoke test is `just smoke sos-steps`.
4. **Fit the industrialization primary estimator** — `Rscript code/_experiments/industrialization-village-two-stage.R` (the two-stage cluster-honest village-level meta-GAM; writes `outputs/_experiments/industrialization-village-two-stage/`).
5. **Run the supplementary analyses** — once the fits from step 3–4 exist, these are independent readouts (each `Rscript code/_experiments/<name>.R`):
   - `calibration-cluster-honest.R`, `calibration-sos-bayesian.R` — calibration of the cluster-honest estimator (frequentist proxy + Bayesian twin).
   - `power-curves-pa-bone.R`, `power-curves-industrialization.R` — a-priori power sweeps; `power-curves-figure.R` renders them; `power-landmarks-table.R` extracts landmark power.
   - `effect-size-probability-fig.R` — effect-size-vs-posterior-probability ("reverse-ROPE") curves.
   - `pa-contrast-effects.R`, `pa-bone-village-fixed-effects.R` — locked-contrast per-person effects and 95% upper-bound magnitudes for the PA→bone nulls.
   - `conf-slope-summaries.R` — confounder-DAG sensitivity magnitudes; `industrialization-linear-summaries.R` — industrialization linear-projection slopes.
   - `age-subset-amef.R` — young/old age-subset AMEF sensitivity.
   - `age-anchors.R`, `outcome-noise-anchors.R` — reference-scale anchors (age-equivalent and outcome-variability yardsticks).
   - `prior-predictive-check.R` — prior-predictive validation for the GAM specs.
   - `walking-experiment-rope.R` — the acute-walking-experiment reverse-ROPE analysis (needs `compiled-walking-experiment-18june2026.csv`).
6. **Assemble figures and tables** — `code/_final/`: `figures.R` (main + supplementary AERF/AMEF figures), `figure-industrialization.R` (industrialization figure), `supp-slope-table.R` and `supp-industrialization-table.R` (supplementary tables). Outputs land in `outputs/figures/final/` and `outputs/tables/`.
7. **DAG figures** — the four causal diagrams are TikZ/`dagitty` sources in `dags/` (`.dag` = graph, `.tex` = figure). Build the `.tex` files with XeLaTeX to produce the mediator, confounder, and industrialization DAG figures.

`code/_experiments/run-village-analyses.sh` records the driver order for the village-level analyses and can be consulted as a worked example.

## Repository layout

```
.
├── _targets.R                  # targets pipeline manifest (15 specs × 4 targets)
├── justfile                    # command surface (front door to the pipeline)
├── renv.lock, renv/, .Rprofile # pinned R environment
├── code/
│   ├── _startup/               # session bootstrap: packages, data prep, spec registry, helpers
│   ├── *.R                     # the 9 reported Tier-2 analysis scripts
│   ├── _experiments/           # supplementary analyses backing the paper's SM figures/tables
│   └── _final/                 # figure and table assembly
├── dags/                       # causal DAG sources (.dag graphs + .tex figures)
└── data/processed/             # (empty) — place the analytical CSVs here; see data/README.md
```

The fitted model objects (`outputs/models/`, `_targets/` cache) are **not** included: they are large and fully regenerable by running the pipeline above.

## License

- **Code** — MIT License (see `LICENSE`).
- **Data** — Creative Commons Attribution 4.0 International (CC-BY-4.0), when included (see `LICENSE-data.md`), subject to the OA HeLP access terms.
