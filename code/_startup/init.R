###############################################################################
# Session-start orchestrator. Sourced by every analysis script. Loads the six
# per-role _startup/ files in dependency order:
#
#   functions.R          first: packages.R calls load_pkgs(), data.R calls
#                        prep_dat(), both defined here
#   packages.R           attaches the core set
#   options.R            session options, including the band-construction switch
#   data.R               builds `dat`
#   specifications.R     the model registry and its constants
#   pipeline-helpers.R   last: its functions reference `model_templates` and the
#                        MCMC constants from specifications.R
#
# specifications.R does not require `dat` to exist. Its formulas name columns as
# symbols inside brms::bf(), which are not evaluated until a fit, and the data
# is referenced by the string "dat". `_targets.R` relies on this: it sources
# functions.R, specifications.R and pipeline-helpers.R with no data.R at all.
###############################################################################


source(here::here("code", "_startup", "functions.R"))
source(here::here("code", "_startup", "packages.R"))
source(here::here("code", "_startup", "options.R"))
source(here::here("code", "_startup", "data.R"))
source(here::here("code", "_startup", "specifications.R"))
source(here::here("code", "_startup", "pipeline-helpers.R"))
