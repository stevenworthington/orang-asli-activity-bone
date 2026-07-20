###############################################################################
# wallace-bone-turnover -- targets pipeline manifest
#
# Use `just fit-all` (or `targets::tar_make()` directly) to materialize.
# Use `just status` to see what's stale; `just viz` for the dependency graph.
#
# Pipeline shape: for each of the 9 specs in model_templates (6 causal
# PA → bone + 3 industrialization), four targets are generated dynamically
# via `tarchetypes::tar_map`:
#
#   <spec>_fit              -- fitted brms object (the long-running step)
#   <spec>_pred_draws       -- AERF posterior draws (level)
#   <spec>_slope_draws      -- AMEF posterior draws (first derivative)
#   <spec>_curvature_draws  -- finite-difference curvature draws (second derivative)
#
# 36 targets total. The pipeline invalidates downstream targets when any of
# (data file, spec registry, helper functions) changes, so refits happen
# automatically when the spec is edited but never when an unrelated script is
# touched.
#
# Design notes on startup cost:
# This file is parsed on EVERY `tar_*()` call (status, manifest, outdated,
# load, make...). Anything sourced or attached at the top runs every time.
# To keep the front door fast, only the packages needed to BUILD the manifest
# are loaded here (targets, tarchetypes, brms for bf() / set_prior() inside
# specifications.R). Heavier per-worker packages (cmdstanr, marginaleffects,
# etc.) are attached just-in-time via `tar_option_set(packages = ...)`.
###############################################################################


suppressPackageStartupMessages({
  library(targets)
  library(tarchetypes)
  library(here)
  library(brms)     # needed at manifest build time for bf() / set_prior() in specifications.R
  library(dplyr)
  library(tidyr)
})


# Source the spec registry + pipeline helpers. None of these load heavy
# packages or read data at sourcing time -- they're pure code definitions.
source(here("code", "_startup", "functions.R"))
source(here("code", "_startup", "specifications.R"))
source(here("code", "_startup", "pipeline-helpers.R"))


# Agent-environment escape hatch: crew workers can spin at 100% CPU without ever
# making progress when spawned from the automation/sandbox context (the mirai
# worker<->controller connection never forms). Set TAR_SEQUENTIAL=1 to run all
# targets in the main R process (no crew). Interactive runs keep the crew default.
.tar_deploy <- if (nzchar(Sys.getenv("TAR_SEQUENTIAL"))) "main" else "worker"

tar_option_set(
  # Packages workers attach when running targets that need them.
  packages = c("brms", "cmdstanr", "tidyverse", "tidybayes",
               "marginaleffects", "mgcv", "splines",
               "HDInterval", "posterior", "ggdist", "janitor"),
  format             = "qs",
  memory             = "transient",
  garbage_collection = TRUE,
  storage            = .tar_deploy,
  retrieval          = .tar_deploy,
  deployment         = .tar_deploy
)


# ---- Pipeline ---------------------------------------------------------------

# Per-spec target template, expanded via tar_map() to produce four targets
# (fit, pred_draws, slope_draws, curvature_draws) for each of the 9 specs.
spec_targets <- tar_map(
  values = tibble::tibble(spec_key = names(model_templates)),
  names  = spec_key,

  tar_target(
    fit,
    {
      file.exists(specs_file)  # invalidate on spec change
      fit_one_spec(spec_key, dat, mode = "full")
    }
  ),

  tar_target(
    pred_draws,
    aerf_draws(fit, model_templates[[spec_key]], dat)
  ),

  tar_target(
    slope_draws,
    amef_draws(fit, model_templates[[spec_key]], dat)
  ),

  tar_target(
    curvature_draws,
    curvature_draws_from_amef(slope_draws, model_templates[[spec_key]])
  )
)


list(

  # -- Inputs ----------------------------------------------------------------

  tar_target(
    data_file,
    here("data", "processed", "Orang-Asli-pa-vs-bone-60126.csv"),
    format = "file"
  ),
  tar_target(
    specs_file,
    here("code", "_startup", "specifications.R"),
    format = "file"
  ),

  # Load + prep the dataset as a target so downstream fits get invalidated
  # when the CSV changes (file-format dependency on `data_file`).
  tar_target(
    dat,
    {
      readr::read_csv(data_file, show_col_types = FALSE) |>
        janitor::clean_names(case = "snake") |>
        dplyr::mutate(
          sex                   = factor(sex, levels = c("female", "male")),
          village_id            = factor(village_id),
          tibia_sos_1k          = tibia_sos / 1000,
          osteocalcin_pg_ml_10k = osteocalcin_pg_ml / 10000,
          ad_steps_1k           = ad_tot_step_count_0_24hr / 1000,
          # Body composition (z-scored). Added 2026-05-26 for the confounder-
          # DAG variant adjustment set; MUST mirror data.R since _targets.R
          # has its own data-prep block (convention duplication).
          fat_mass_kg_z         = as.numeric(scale(fat_mass_kg)),
          fat_free_mass_kg_z    = as.numeric(scale(fat_free_mass_kg))
        )
    }
  ),


  # -- Per-spec fits + draws (9 specs x 4 targets = 36 targets) -------------

  spec_targets

)
