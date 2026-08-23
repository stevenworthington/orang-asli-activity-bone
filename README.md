# Replication code for "Reduced physical activity is unlikely to explain declining bone health during industrialization"

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22062619.svg)](https://doi.org/10.5281/zenodo.22062619)

Replication code for:

> Cecilia M. T. Sena, Steven Worthington, Thomas S. Kraft, Fabiano T. Amorim, Kumaresh Louis Christopher, Melissa Emery Thompson, Kamal Solhaimi bin Fadzil, Edwin Ser Ze Heng, Nicholas B. Holowka, Tan Bee Ting A/P Tan Boon Huat, Janet L. Huebner, Tracy L. Kivell, Virginia B. Kraus, Yvonne A. L. Lim, Colin Nicholas, Marissa Ramirez, Izandis bin Mohd Sayed, Kar Lye Tam, Marina M. Watowich, Vivek V. Venkataraman, Amanda J. Lea, Ian J. Wallace. *Reduced physical activity is unlikely to explain declining bone health during industrialization.*
>
> Cecilia M. T. Sena and Steven Worthington contributed equally. Corresponding author: Ian J. Wallace (iwallace@unm.edu).

The analysis code in this repository was written by Steven Worthington, and `CITATION.cff` records that. The byline above is the article's, not the code's — please cite the article for the study itself, which is what `CITATION.cff` names as the preferred citation.

The study is a causal analysis of physical activity and bone health in the Orang Asli Health and Lifeways Project (OA HeLP) cohort (n = 1,007), across two analysis sets:

- **Physical activity → bone.** Daily step count and mean daily ENMO against tibial speed-of-sound (SOS), CTX-1, and osteocalcin; within-community (community fixed-effects) estimand. Adjustment (mediator DAG): age (in a tensor smooth) + sex + functional status + pregnancy/lactation + smoking + alcohol.
- **Industrialization → outcome.** A community-level industrialization index against tibial SOS, daily steps, and mean daily ENMO; across-community estimand, from a one-stage hierarchical model with a community random intercept. Adjustment: age + sex.

The reported causal estimands are the average exposure–response function (AERF) and its first derivative, the average marginal-effect function (AMEF). Where the AERF is approximately linear it is summarized by its linear-projection slope, reported with a 95% highest-posterior-density interval and the posterior probability of the hypothesized sign. Curvature (the AERF's second derivative) is computed as a diagnostic of departure from linearity.

The two analysis sets use different model forms because their exposures sit at different levels. Activity varies within a community, so community identity enters as fixed indicators and the estimate comes from within-community contrasts only. The industrialization index is a community attribute — it takes one value per community, about 25 of them — so it carries a community random intercept instead, the smooths are additive rather than tensor because an age × index interaction is informed by only those ~25 distinct index values and is too weakly determined to resolve, and the uncertainty on it is governed primarily by how many communities were sampled rather than how many individuals. Individual observations still inform the outcome means and the variance components; what they do not do is supply independent realisations of the industrialization gradient.

## Data availability

**The data are not included in this repository. They are available upon request from the corresponding author, Ian J. Wallace (iwallace@unm.edu).** `data/processed/` is kept as an empty directory (with a `README.md` describing the two CSVs the code expects, and the columns it reads from each). The data are individual-level records from the Orang Asli, an Indigenous population of Peninsular Malaysia; access is governed by the Orang Asli Health and Lifeways Project (OA HeLP) and the associated ethics approvals, so the files are shared on request rather than distributed publicly.

Once the data files are in place under `data/processed/`, the pipeline below reproduces every fitted model, figure, and table.

**A note on "village" and "community".** The data file identifies a participant's community with a column called `village_id`, and the object names and output paths derived from it keep that word. Everything written for a reader — this README, the code comments, the figure labels, the manuscript — says *community* instead, because several of the sampled settlements are sizeable peri-urban areas rather than villages, so *community* is the more accurate word for the unit. Both refer to the same thing. The column is deliberately left as it is, so the code matches the data file you receive rather than one you would have to rename first. `data/README.md` lists it with the other columns the analysis reads.

## Running it without the data

The pipeline can be run end to end on a **synthetic** dataset, so the code is executable by anyone, not only by those granted access to the real data.

```sh
just example        # synthetic data + static gate + one spec, under a minute
just example-full   # all 15 fits and their draws, the model-derived figures and tables
```

`just example-full` takes roughly 10-15 minutes on a first run, most of it Stan compilation, and nearer 5 once the compile cache is warm. It builds Figures 3 and 4 and two supplementary figures, plus the slope tables; the power sweeps, the effect-size figure and the age-subset refits are separate scripts, left out because they cost far more than the fits do.

**Read `code/_experiments/make-example-data.R` before trusting anything these produce.** Nothing in that generator is derived from the real data: no model was fitted to the real participants, nothing was resampled or perturbed, and every value comes from a distribution written out in full in the script. What it borrows from the study is only its published design — the number of participants and communities, which communities carry which measurement, and the approximate scale of each variable. The effects it plants are deliberately *not* the study's: the physical-activity effect on bone is set clearly positive there, where the paper reports a null. **Results computed from it describe the generator, not the Orang Asli.**

Both example recipes sample in smoke mode — 4 chains × 1,500 iterations, 500 warm-up, thinned by 5, so **800 posterior draws** against the 10,000 the paper reports. A simultaneous band widens its pointwise mass until whole curves are contained, and on a posterior this thin that calibration runs close to the sample's extremes, so bands drawn from 800 draws look visibly rougher than the published figures. That is the thin posterior showing, not a fault.

The generated CSVs land in `data/processed/` under the same filenames the real data use, so the pipeline needs no switch to find them. To keep that from becoming a hazard, the generator records what it wrote in `data/processed/.synthetic-example` and **refuses to overwrite any file it did not write itself** — so running the example on a machine that already holds the real data stops with an error instead of destroying them. That marker is also how you tell the two apart: if it lists a file, that file is synthetic. `data/processed/` is gitignored, so neither can be committed.

## Software environment

- **R 4.6.0** with the exact package versions pinned in `renv.lock`.
- **CmdStan** via `{cmdstanr}` is the Stan backend for `{brms}`; install it after restoring the library. The manuscript used CmdStan **2.38.0**; pin the install to that version for a faithful reproduction.

Restore the environment from the project root:

```r
# from R, in the repository root
renv::restore()          # rebuilds the pinned library from renv.lock
cmdstanr::install_cmdstan(version = "2.38.0")   # one-time: the version the manuscript used
```

The model fits are computationally heavy: each of the 15 specifications runs 10 chains of 7,000 iterations (2,000 warmup, thinned by 5) at `adapt_delta = 0.999`. Budget accordingly and run them on a machine with adequate cores and RAM.

## Reproduce

`code/_experiments/run-all-analyses.sh` runs steps 3–6 below in dependency order, unattended, and is the quickest way to regenerate the model-derived results. It stops if the fits fail, reports which steps failed, and exits nonzero. It does **not** restore the environment (step 1), place the data (step 2), build the DAG figures (step 7) or export PNGs — those need decisions or a toolchain beyond R. What follows is the full sequence, explained. `just` recipes are the front door (`just` with no arguments lists them); the underlying `Rscript` calls are shown for reference.

1. **Restore the environment** — `renv::restore()` (see above).
2. **Place the data** — copy the two CSVs into `data/processed/` (see `data/README.md`).
3. **Check the source** — `just check` runs `code/_checks/verify-startup.R`, which fits nothing and takes seconds. It asserts that the 15 specifications carry the families and constants they are supposed to, that the prediction grids are built from the right values, that `build_priors()` reproduces a hand computation, and that the documented switches actually switch something. Run it before committing to step 4.
4. **Fit the models** — `just fit-all` runs the `_targets.R` pipeline: 15 specifications (6 mediator-DAG + 6 confounder-DAG + 3 industrialization), each producing a fitted `brms` object plus AERF / AMEF / curvature posterior draws. Only out-of-date targets re-run. Individual specs can also be run directly (`just fit-script sos-steps`, or `Rscript code/sos-steps.R`); a fast smoke test is `just smoke sos-steps`.
5. **Assemble figures and tables** — `code/_final/`: `figures.R` (Figure 4 plus the confounder-DAG and age-conditional supplementary grids), `figure-industrialization.R` (Figure 3), `supp-slope-table.R` and `supp-industrialization-table.R` (the two supplementary slope tables). Outputs land in `outputs/figures/final/` and `outputs/tables/`. `just figures-png` re-exports the PDFs to 300-dpi PNG (needs `pdftoppm`, from poppler).
6. **Run the supplementary analyses** — each is `Rscript code/_experiments/<name>.R`. Those in the first group read the fits from step 4, so they must run after it; the rest need only the data.

   Reading the fits:
   - `pa-contrast-effects.R` — the main-text effect-magnitude bounds for the six physical-activity → bone nulls: the linear-projection slope over the reported contrast, and U95, the 95% upper bound on the magnitude. `conf-slope-summaries.R` is the confounder-DAG counterpart. Both read the reference scales in `_reference-scales.R`, so they quote the same denominators.
   - `industrialization-linear-summaries.R` — linear-projection slopes and linearity thresholds for the three industrialization analyses.
   - `age-subset-amef.R` — young/old age-subset sensitivity (refits each spec below and above age 35).
   - `effect-size-probability-fig.R` — effect-size-vs-posterior-probability ("reverse-ROPE") curves for all nine analyses.

   Self-contained:
   - `age-anchors.R`, `outcome-noise-anchors.R` — the reference scales the effect-size bounds are read against: an age-equivalent yardstick, and the assay/measurement resolution floor per outcome.
   - `power-curves-pa-bone.R`, `power-curves-industrialization-onestage.R` — a-priori power sweeps, one per analysis set, each using a frequentist twin matched to the production likelihood and estimator. `power-curves-figure.R` renders them; `power-landmarks-table.R` extracts landmark power.
   - `calibration-pa-bone.R`, `calibration-onestage-production.R` — false-positive rate, coverage and power for the two reported estimators, simulated on the real study design at the production priors, on the linear-projection estimand the paper reports. These are what license the paper's magnitude bounds: a bound is only as good as the coverage of the interval it comes from. Each invocation calibrates **one** arm — select it with `CAL_SPEC`/`CAL_COND` and `CAL_OUTCOME` — so the full set needs the loop the driver runs; the bare defaults cover one arm each.
   - `walking-experiment-rope.R` — the acute walking-experiment reverse-ROPE analysis (needs `walking-experiment-data.csv`).

7. **DAG figures** — the four causal diagrams are TikZ/`dagitty` sources in `dags/` (`.dag` = graph, `.tex` = figure). Building them requires a LaTeX distribution with XeLaTeX (e.g. TeX Live or MacTeX). Compile the `.tex` files to produce the main-text three-hypothesis composite and the mediator, confounder, and industrialization DAG figures.

### Band construction

Bands drawn over a *curve* — the nested grey ribbons on every AERF and AMEF panel, and the red envelope on the linear projection of the AERF — are **simultaneous** credible bands: the stated level is the probability that an entire posterior curve lies inside the band, not the per-x probability.

One band is not, and deliberately so. The horizontal red band in each AMEF panel shows the interval for the *constant slope* of that linear projection. A constant slope is a single scalar rather than a curve, so there is no simultaneity to calibrate and it is a plain HPDI. The figure captions draw the same distinction: SHPDI for the projection over the AERF, HPDI for the constant slope beneath it.

Two constructions of the simultaneous bands are implemented, and the reported figures use the first:

- **HPDI** (the default) — nested pointwise highest-density intervals, whose common mass is raised until whole curves are contained at the stated level. "Pointwise" describes the interval at each grid point; the calibration target is simultaneous coverage of the curve. Asymmetric wherever the posterior is skewed, which is what makes it a highest-density region.
- **symmetric** — mean ± c\* × pointwise SD, with c\* calibrated on the sup-norm. Symmetric about the mean however the posterior is shaped, so not a highest-density region. Kept so the two can be compared directly.

Setting `BAND_TYPE=symmetric` in the environment rebuilds every simultaneous band the second way; `just figures-band-comparison` builds both sets side by side. Each band carries its construction in an `interval_type` column, so a saved band always says how it was made.

The per-script diagnostic figures under `outputs/figures/working/` are outside all of this: they use *pointwise* intervals at three levels and are labelled as such. They exist to check a fit, not to report one.

## Repository layout

```
.
├── _targets.R                  # targets pipeline manifest (15 specs × 4 targets)
├── justfile                    # command surface (front door to the pipeline)
├── renv.lock, renv/, .Rprofile # pinned R environment
├── code/
│   ├── _startup/               # session bootstrap: packages, data prep, spec registry, helpers
│   ├── _checks/                # static verification of the startup set (`just check`)
│   ├── *.R                     # the 9 reported analysis scripts
│   ├── _experiments/           # supplementary analyses backing the paper's figures and tables
│   └── _final/                 # figure and table assembly
├── dags/                       # causal DAG sources (.dag graphs + .tex figures)
└── data/processed/             # (empty) — place the analytical CSVs here; see data/README.md
```

The fitted model objects (`outputs/models/`, `_targets/` cache) are **not** included: they are large and fully regenerable by running the pipeline above.

## License

- **Code** — MIT License (see `LICENSE`). This repository contains code only.
- **Data** — not distributed here; available on request from the corresponding author under the Orang Asli Health and Lifeways Project (OA HeLP) data-governance terms (see `data/README.md`).
