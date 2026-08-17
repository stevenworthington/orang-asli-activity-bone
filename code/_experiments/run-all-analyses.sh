#!/usr/bin/env bash
###############################################################################
# Unattended end-to-end regeneration of every reported analysis, in dependency
# order. It is both a driver and the executable statement of what depends on
# what: the figures and tables read the targets cache, so step 1 must finish
# before anything else runs.
#
# Steps continue on error, so one failure never blocks the rest; per-step exit
# status is appended to the STATUS file and each step's output goes to its own
# log. It overwrites outputs/ in place.
#
# Continuing past a failure is deliberate, but the run must not then report
# success: failed step names are collected and the script prints FAILED and
# exits nonzero if any step failed. Step 1 is treated as a hard dependency --
# everything downstream reads the cache it writes, so if it fails the rest would
# assemble figures and tables from whatever fits happen to be there already.
#
# It does NOT build the DAG figures or export PNGs; those are steps 7 and the
# figures-png recipe in the README, and both need a toolchain beyond R.
#
# The heavy step by far is step 1 (15 brms fits). Launch the whole thing
# detached so it survives the shell, e.g.:
#   nohup bash code/_experiments/run-all-analyses.sh &
###############################################################################
set -u
# Resolve the project root from this script's own location (code/_experiments/),
# so the script runs from any checkout and under a job runner that starts in an
# arbitrary working directory.
cd "$(dirname "$0")/../.." || exit 1
# Robust env for a background daemon, which may not inherit the interactive PATH.
export PATH="/usr/local/bin:/opt/homebrew/bin:/Library/TeX/texbin:$PATH"
# Point CMDSTAN at your cmdstan installation if it is not on the default path.
export CMDSTAN="${CMDSTAN:-$HOME/.cmdstan/cmdstan-2.38.0}"
# Run targets in the main R process rather than through crew workers. Slower,
# but it keeps the whole run in one process, which is what a detached job wants.
export TAR_SEQUENTIAL=1

LOG="${TMPDIR:-/tmp}/oa-analyses"; mkdir -p "$LOG"
S="$LOG/STATUS.txt"; : > "$S"
say() { echo "[$(date '+%m-%d %H:%M')] $*" | tee -a "$S"; }

# A plain string, not an array: macOS still ships bash 3.2, where expanding an
# empty array under `set -u` is an unbound-variable error.
FAILED_STEPS=""

run() {                       # run <name> <Rscript-or-cmd...>
  local name="$1"; shift
  say "START $name"
  "$@" < /dev/null > "$LOG/$name.log" 2>&1
  local rc=$?
  say "END   $name  rc=$rc  (log: $LOG/$name.log)"
  # `say` pipes through tee and would otherwise supply this function's exit
  # status, which is always 0 -- the reason an earlier version of this script
  # could report COMPLETE after a failed fit.
  [ "$rc" -eq 0 ] || FAILED_STEPS="$FAILED_STEPS $name"
  return "$rc"
}

say "=== analysis run begin ==="

# 0. Static gate. Fits nothing, takes seconds, and fails loudly on a source
#    inconsistency that would otherwise surface hours into step 1.
run check Rscript code/_checks/verify-startup.R

# 1. The 15 fits and their AERF / AMEF / curvature draws. targets persists all
#    60 objects in _targets/objects/, and skips any that is already up to date.
if ! run tar-make Rscript -e 'targets::tar_make(callr_function = NULL)'; then
  say "FATAL: the fits failed; everything below reads their output. Stopping."
  say "=== analysis run FAILED (tar-make) ==="
  exit 1
fi
say "  fits present: $(ls _targets/objects/fit_* 2>/dev/null | wc -l | tr -d ' ')"

# 2. PA -> bone figures (Figure 4, the confounder-DAG and age-conditional
#    supplementary grids) + spec-summary.csv.
run figures Rscript code/_final/figures.R

# 3. Industrialization figure and its slope numbers.
run fig-industrialization Rscript code/_final/figure-industrialization.R
run ind-table             Rscript code/_final/supp-industrialization-table.R
run ind-linear            Rscript code/_experiments/industrialization-linear-summaries.R

# 4. Supplementary slope table for PA -> bone (slopes + shape diagnostics).
run slope-table Rscript code/_final/supp-slope-table.R

# 5. Age-subset sensitivity (refits age<35 / >=35 from the registry).
run age-subset Rscript code/_experiments/age-subset-amef.R

# 6. Effect-magnitude reporting: the bounds the Results quote, the reverse-ROPE
#    figure, and the reference scales both are read against.
run pa-contrast      Rscript code/_experiments/pa-contrast-effects.R
run conf-slopes      Rscript code/_experiments/conf-slope-summaries.R
run effsize-fig      Rscript code/_experiments/effect-size-probability-fig.R
run noise-anchors    Rscript code/_experiments/outcome-noise-anchors.R
run age-anchors      Rscript code/_experiments/age-anchors.R

# 7. A-priori power. The two sweeps write the grids; the figure and the
#    landmark table read them, so they run after.
run power-pa    Rscript code/_experiments/power-curves-pa-bone.R
run power-ind   Rscript code/_experiments/power-curves-industrialization-onestage.R
run power-fig   Rscript code/_experiments/power-curves-figure.R
run power-table Rscript code/_experiments/power-landmarks-table.R

# 8. Calibration of the two reported estimators. Independent of everything
#    above -- these simulate their own data on the real study design -- but
#    slow, so they run last. Each script calibrates ONE arm per invocation, so
#    every arm is enumerated: three PA specs x two conditions, and three
#    industrialization outcomes at the production prior cell. Running either
#    script bare covers only its own defaults.
for spec in sos-steps ctx-steps osteo-steps; do
  for cond in null effect; do
    run "calib-pa-$spec-$cond" \
      env CAL_SPEC="$spec" CAL_COND="$cond" Rscript code/_experiments/calibration-pa-bone.R
  done
done
for outcome in sos steps enmo; do
  run "calib-ind-$outcome" \
    env CAL_OUTCOME="$outcome" Rscript code/_experiments/calibration-onestage-production.R
done

# 9. The acute walking experiment, a separate dataset and a standalone readout.
run walking Rscript code/_experiments/walking-experiment-rope.R

if [ -z "$FAILED_STEPS" ]; then
  say "=== analysis run COMPLETE ==="
else
  say "=== analysis run FAILED:$FAILED_STEPS ==="
  say "    (per-step logs in $LOG)"
  exit 1
fi
