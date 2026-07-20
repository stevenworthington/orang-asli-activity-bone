###############################################################################
# Session-start orchestrator. Sourced by every Tier-2 analysis script. Loads
# the six per-role _startup/ files in the canonical order required by their
# dependencies (functions first so packages.R can use load_pkgs(); data.R uses
# tidyverse from packages.R; specifications.R references column names from
# `dat`; pipeline-helpers.R references `model_templates` from specifications.R).
# See `~/.config/agents/CONVENTIONS.md` § "_startup/ conventions".
###############################################################################


library(here)
source(here("code", "_startup", "functions.R"))
source(here("code", "_startup", "packages.R"))
source(here("code", "_startup", "options.R"))
source(here("code", "_startup", "data.R"))
source(here("code", "_startup", "specifications.R"))
source(here("code", "_startup", "pipeline-helpers.R"))
