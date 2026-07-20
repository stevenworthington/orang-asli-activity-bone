#!/usr/bin/env bash
###############################################################################
# Overnight unattended run: village-fixed-effects canonical regeneration.
#
# Runs the full pipeline end-to-end after the 2026-06-15 village-FE registry
# change. Every R step is pty-wrapped (script -q /dev/null ... < /dev/null) to
# avoid the no-controlling-terminal 100% CPU spin; targets runs SEQUENTIALLY
# (TAR_SEQUENTIAL=1 -> no crew, no ttyless worker spin). Steps continue on
# error so one failure never blocks the rest; per-step status is appended to
# the STATUS file. Pre-change results were archived to
# .graveyard/2026-06-15-pre-village-fe/ before this runs.
#
# Launch (survives the session) via pueue, sandbox off:
#   pueue add -- bash code/_experiments/overnight-village-fe-run.sh
###############################################################################
set -u
cd "/Users/sworthin/Library/CloudStorage/Dropbox/workspace/papers/wallace-bone-turnover" || exit 1
# Robust env for pueue's daemon (which may not inherit the interactive PATH).
export PATH="/usr/local/bin:/opt/homebrew/bin:/Library/TeX/texbin:$PATH"
export CMDSTAN="$HOME/.cmdstan/cmdstan-2.38.0"
export TAR_SEQUENTIAL=1

LOG=/tmp/claude/overnight; mkdir -p "$LOG"
S="$LOG/STATUS.txt"; : > "$S"
say() { echo "[$(date '+%m-%d %H:%M')] $*" | tee -a "$S"; }

run() {                       # run <name> <Rscript-or-cmd...>
  local name="$1"; shift
  say "START $name"
  script -q /dev/null "$@" < /dev/null > "$LOG/$name.log" 2>&1
  local rc=$?
  say "END   $name  rc=$rc  (log: $LOG/$name.log)"
}

say "=== overnight village-FE run begin ==="

# 1. Full canonical re-fit (FE for the 12 PA->bone specs; industrialization
#    specs re-fit identically since specs_file changed). targets persists every
#    fit + derivative (pred/slope/curvature) in _targets/objects/.
run tar-make Rscript -e 'targets::tar_make(callr_function = NULL)'
say "  fits present: $(ls _targets/objects/fit_* 2>/dev/null | wc -l | tr -d ' ')"

# 2. Canonical figures (Fig 4 PA->bone, Fig 3 urb, conf + age-conditional supp)
#    + spec-summary.csv. Reads FE draws via tar_read_raw.
run figures Rscript code/_final/figures.R

# 3. Supplementary Table 1 (slopes + shape diagnostics). Reads spec-summary.csv.
run slope-table Rscript code/_final/supp-slope-table.R

# 4. Age-subset sensitivity (refits age<35 / >=35; inherits FE from the registry).
run age-subset Rscript code/_experiments/age-subset-amef.R

# 5. Effect-size summaries on the FE fits (contrast table + reverse-ROPE figure;
#    both read the now-FE pred_draws from _targets/objects/).
run pa-contrast    Rscript code/_experiments/pa-contrast-effects.R
run effsize-fig    Rscript code/_experiments/effect-size-probability-fig.R

# 6. Industrialization cluster-honest analyses (saved to outputs/_experiments/).
run option-h Rscript code/_experiments/set2-option-h-meta-gam.R
run option-g Rscript code/_experiments/set2-option-g-cluster-bootstrap.R
run option-b Rscript code/_experiments/set2-option-b-ranef-diagnostic.R

say "=== overnight village-FE run COMPLETE ==="
