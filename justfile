# orang-asli-activity-bone -- project command surface
#
# `just` is the project's front door; the heavier `targets` pipeline runner
# lives underneath. Run `just` with no args to see all recipes.

# Show the recipe list when invoked with no args.
default:
    @just --list


# ---- targets pipeline (15 specs: 6 PA-bone + 6 confounder-DAG variants + 3 industrialization) ----

# Run the full targets pipeline (15 fits + their AERF/AMEF/curvature draws; 60 targets).
# Skips any target whose inputs haven't changed since the last run.
fit-all:
    Rscript -e 'targets::tar_make()'

# Run a single target and its dependencies.
# Usage: just fit-one fit_sos_steps
fit-one TARGET:
    Rscript -e 'targets::tar_make({{TARGET}})'

# Run just the 15 brms fits (skip the draw extraction stages).
fit-models-only:
    Rscript -e 'targets::tar_make(matches("^fit_"))'

# Show which targets are out-of-date vs cached.
status:
    Rscript -e 'targets::tar_outdated()'

# Show the dependency graph (opens an interactive HTML widget).
viz:
    Rscript -e 'targets::tar_visnetwork()'

# Show the static manifest (target list + their dependencies).
manifest:
    Rscript -e 'print(targets::tar_manifest(), n = Inf)'

# Pull a saved target back into an interactive R session.
# Usage: just load pred_draws_sos_steps
load TARGET:
    Rscript -e 'print(targets::tar_read({{TARGET}}))'

# Destroy the targets cache. Forces all targets to re-run next tar_make().
# Confirms before nuking.
clean-targets:
    Rscript -e 'targets::tar_destroy(ask = TRUE)'


# ---- Per-spec Rscript wrappers (parallel path to targets) ----------------

# Fit a single spec via the named Rscript at code/<stem>.R. Each script also
# saves model.Rdata + pred-draws.Rdata + slope-draws.Rdata + curvature-draws.Rdata
# + a working aerf-amef.pdf to outputs/{models,figures/working}/<stem>/.
# Usage: just fit-script sos-steps
fit-script STEM:
    Rscript code/{{STEM}}.R

# Smoke-test a single script with fast 4-chain sampling.
# Usage: just smoke sos-steps
smoke STEM:
    SAMPLING_MODE=smoke Rscript code/{{STEM}}.R


# ---- Runnable example (synthetic data) -----------------------------------
#
# The real data are not distributed here, so these recipes let the pipeline be
# run end to end on a synthetic dataset instead. Read
# code/_experiments/make-example-data.R before trusting anything they produce:
# nothing in it is derived from the real data, and the effects it plants are
# deliberately NOT the study's. The numbers, figures and tables that come out
# describe the generator, not the Orang Asli.
#
# Sampling is smoke mode throughout -- 4 chains x 1,500 iterations, 500 warm-up,
# thinned by 5, so 800 draws against the 10,000 the paper reports. Bands drawn
# from 800 draws look visibly rougher than the published ones; that is the thin
# posterior showing, not a fault.

# Write the synthetic dataset into data/processed/ (gitignored).
example-data:
    Rscript code/_experiments/make-example-data.R

# Synthetic data, static gate, then one spec fitted in smoke mode. Under a
# minute once Stan has compiled that model.
#
# Smallest end-to-end check that the pipeline runs (synthetic data).
example: example-data
    Rscript code/_checks/verify-startup.R
    SAMPLING_MODE=smoke Rscript code/sos-steps.R

# The whole pipeline on synthetic data. Roughly 10-15 minutes on a first run,
# most of it Stan compilation; nearer 5 once the compile cache is warm.
#
# SAMPLING_MODE is part of each target's command, so switching between this and
# `just fit-all` correctly invalidates the cache and refits rather than mixing
# smoke and production fits.
#
# All 15 fits and their draws, the model-derived figures and slope tables.
# The power sweeps and age-subset refits are separate scripts and are not run
# here; they cost far more than the fits do.
example-full: example-data
    Rscript code/_checks/verify-startup.R
    SAMPLING_MODE=smoke TAR_SEQUENTIAL=1 Rscript -e 'targets::tar_make()'
    Rscript code/_final/figures.R
    Rscript code/_final/figure-industrialization.R
    Rscript code/_final/supp-slope-table.R
    Rscript code/_final/supp-industrialization-table.R


# ---- Diagnostics ---------------------------------------------------------

# Show renv consistency state. Should print "No issues found".
renv-status:
    Rscript -e 'renv::status()'


# ---- Figure export -------------------------------------------------------

# The R generators write PDF only; PNGs are what get inserted into the
# manuscript. Run this after any figure regeneration, so the two cannot fall
# out of step.
#
# supp-fig-5-power-curves is NOT listed: power-curves-figure.R writes its own
# PNG directly (at 200 dpi). The four DAG figures are built by the dags/
# pipeline and are not model-derived.
#
# Re-export the model-derived final figures to 300-dpi PNG.
figures-png:
    #!/usr/bin/env bash
    set -euo pipefail
    cd outputs/figures/final
    for f in fig-3-urb fig-4-pa-bone supp-fig-4-effect-size \
             supp-fig-6-pa-bone-conf supp-fig-7-age-conditional-bands \
             supp-fig-8-age-subset; do
      pdftoppm -png -r 300 -singlefile "$f.pdf" "$f"
      echo "  exported $f.png"
    done


# The "symmetric" set (mean +/- c* x pointwise SD) is built first and moved
# aside; the default HPDI set is built second, so the canonical
# outputs/figures/final/ ends up holding the HPDI version -- the one the paper
# reports. Deliberately not parameterised by an output suffix: only one of these
# is the manuscript's figure, and a suffix would make it easy to ship the other.
#
# Build the model-derived figures under both band constructions.
figures-band-comparison:
    #!/usr/bin/env bash
    set -euo pipefail
    cmp_dir=outputs/_experiments/band-type-comparison
    mkdir -p "$cmp_dir/symmetric" "$cmp_dir/HPDI"
    for bt in symmetric HPDI; do
      echo "=== building figures with BAND_TYPE=$bt ==="
      BAND_TYPE=$bt Rscript code/_final/figures.R
      BAND_TYPE=$bt Rscript code/_final/figure-industrialization.R
      for f in fig-3-urb fig-4-pa-bone supp-fig-6-pa-bone-conf supp-fig-7-age-conditional-bands; do
        cp "outputs/figures/final/$f.pdf" "$cmp_dir/$bt/$f.pdf"
      done
    done
    echo "both sets in $cmp_dir; outputs/figures/final/ holds the HPDI version"


# Fits nothing; takes seconds. It asserts that every call into an attached
# package is namespaced, that the 15 specs carry the families and constants the
# specification names, that the postestimation grid is the right size and drawn
# from the right values, that build_priors() reproduces a hand computation, and
# that the documented switches (interval_type, contrast_units, rug_data) switch
# something. Exits non-zero on the first failure.
#
# `CHECK_FIT=1 just check` adds a smoke fit of one community-level spec, the
# only way to confirm re_formula survives the trip through marginaleffects.
#
# Verify the startup set statically. Gate for any long run.
check:
    Rscript code/_checks/verify-startup.R
