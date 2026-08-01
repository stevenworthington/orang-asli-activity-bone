#!/usr/bin/env bash
###############################################################################
# Unattended end-to-end regeneration of the village-level analyses.
#
# Runs the full pipeline end to end. Every R step is
# run under a pseudo-terminal (script -q /dev/null ... < /dev/null) to
# avoid the no-controlling-terminal 100% CPU spin; targets runs SEQUENTIALLY
# (TAR_SEQUENTIAL=1 -> no crew, no ttyless worker spin). Steps continue on
# error so one failure never blocks the rest; per-step status is appended to
# the STATUS file. It overwrites outputs/ in place.
#
# Launch it detached so it survives the shell, e.g.:
#   nohup bash code/_experiments/run-village-analyses.sh &
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
export TAR_SEQUENTIAL=1

LOG="${TMPDIR:-/tmp}/village-analyses"; mkdir -p "$LOG"
S="$LOG/STATUS.txt"; : > "$S"
say() { echo "[$(date '+%m-%d %H:%M')] $*" | tee -a "$S"; }

run() {                       # run <name> <Rscript-or-cmd...>
  local name="$1"; shift
  say "START $name"
  script -q /dev/null "$@" < /dev/null > "$LOG/$name.log" 2>&1
  local rc=$?
  say "END   $name  rc=$rc  (log: $LOG/$name.log)"
}

say "=== village-level analyses run begin ==="

# 1. Full canonical re-fit (FE for the 12 PA->bone specs; industrialization
#    specs re-fit identically since specs_file changed). targets persists every
#    fit + derivative (pred/slope/curvature) in _targets/objects/.
run tar-make Rscript -e 'targets::tar_make(callr_function = NULL)'
say "  fits present: $(ls _targets/objects/fit_* 2>/dev/null | wc -l | tr -d ' ')"

# 2. Canonical figures (Fig 4 PA->bone, Fig 3 industrialization, conf +
#    age-conditional supp) + spec-summary.csv. Reads FE draws via tar_read_raw.
run figures Rscript code/_final/figures.R

# 3. Supplementary Table 1 (slopes + shape diagnostics). Reads spec-summary.csv.
run slope-table Rscript code/_final/supp-slope-table.R

# 4. Age-subset sensitivity (refits age<35 / >=35; inherits FE from the registry).
run age-subset Rscript code/_experiments/age-subset-amef.R

# 5. Effect-size summaries on the FE fits (contrast table + reverse-ROPE figure;
#    both read the now-FE pred_draws from _targets/objects/).
run pa-contrast    Rscript code/_experiments/pa-contrast-effects.R
run effsize-fig    Rscript code/_experiments/effect-size-probability-fig.R

# 6. Industrialization village-level two-stage estimator (PRIMARY for the industrialization analyses;
#    saved to outputs/_experiments/industrialization-village-two-stage/).
run village-two-stage Rscript code/_experiments/industrialization-village-two-stage.R

say "=== village-level analyses run COMPLETE ==="
