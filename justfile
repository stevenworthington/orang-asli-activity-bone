# wallace-bone-turnover -- project command surface
#
# General pattern (per ~/.config/agents/agentic-ai-tool-chain.md §15a):
# `just` is the project's front door; the heavier `targets` pipeline runner
# (§15b) lives underneath. Run `just` with no args to see all recipes.

# Show the recipe list when invoked with no args.
default:
    @just --list


# ---- targets pipeline (9 specs: 6 PA-bone + 3 industrialization) ---------

# Run the full targets pipeline (9 fits + their AERF/AMEF/curvature draws).
# Skips any target whose inputs haven't changed since the last run.
fit-all:
    Rscript -e 'targets::tar_make()'

# Run a single target and its dependencies.
# Usage: just fit-one fit_sos_steps
fit-one TARGET:
    Rscript -e 'targets::tar_make({{TARGET}})'

# Run just the 9 brms fits (skip the draw extraction stages).
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


# ---- Diagnostics ---------------------------------------------------------

# Show renv consistency state. Should print "No issues found".
renv-status:
    Rscript -e 'renv::status()'
