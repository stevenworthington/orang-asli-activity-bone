###############################################################################
# Supplementary slope table for the physical-activity -> bone analyses
# (referenced from the main-text Figure 4 and Supplementary Figure 6 captions).
# 12 rows: the 6 mediator-DAG specs (Figure 4 panels A-F) and the 6 confounder-DAG
# *-conf specs (Supplementary Figure 6 panels A-F), reporting:
#   - linear-projection slope of the AERF + 95% HPDI + P(beta < 0)
#   - linearity threshold: 95th-percentile posterior upper bound on the AERF's
#     maximum deviation from its best linear projection, as % of the AERF range
#   - flatness threshold: 95th-percentile posterior upper bound on the maximum
#     absolute local slope of the AMEF, in outcome units per exposure unit
#
# Both thresholds are defined in the supplement's effect-magnitude section, and
# the flatness column must match the identically-named column in the
# industrialization slope table -- so both are computed the same way, from the
# `slope_draws` target.
#
# All quantities computed from the cached pred_draws / slope_draws, and reported
# in each outcome's native unit -- osteocalcin pg/mL, CTX-1 ng/mL, SOS m/s. These
# are the units the registry's `outcome_label` gives the figure axes and that
# `_experiments/_reference-scales.R` gives the magnitude bounds, so this table
# can be compared against Figure 4 without converting anything.
#
# Outputs:
#   outputs/tables/supp-slope-table.csv   machine-readable
#   stdout                                GitHub-flavored markdown table to paste
###############################################################################


source(here::here("code", "_startup", "init.R"))


# ---- Spec ordering + reporting metadata (Figure 4 panel order A-F) ----

rows <- tibble::tribble(
  ~dag,         ~panel, ~spec_key,          ~outcome_lab,  ~unit,   ~report_div, ~outcome_col,
  "Mediator",   "A",    "osteo-steps",      "Osteocalcin", "pg/mL",           1,  "osteocalcin_pg_ml",
  "Mediator",   "B",    "ctx-steps",        "CTX-1",       "ng/mL",           1,  "ctx1_ng_ml",
  "Mediator",   "C",    "sos-steps",        "Tibial SOS",  "m/s",             1,  "tibia_sos",
  "Mediator",   "D",    "osteo-enmo",       "Osteocalcin", "pg/mL",           1,  "osteocalcin_pg_ml",
  "Mediator",   "E",    "ctx-enmo",         "CTX-1",       "ng/mL",           1,  "ctx1_ng_ml",
  "Mediator",   "F",    "sos-enmo",         "Tibial SOS",  "m/s",             1,  "tibia_sos",
  "Confounder", "A",    "osteo-steps-conf", "Osteocalcin", "pg/mL",           1,  "osteocalcin_pg_ml",
  "Confounder", "B",    "ctx-steps-conf",   "CTX-1",       "ng/mL",           1,  "ctx1_ng_ml",
  "Confounder", "C",    "sos-steps-conf",   "Tibial SOS",  "m/s",             1,  "tibia_sos",
  "Confounder", "D",    "osteo-enmo-conf",  "Osteocalcin", "pg/mL",           1,  "osteocalcin_pg_ml",
  "Confounder", "E",    "ctx-enmo-conf",    "CTX-1",       "ng/mL",           1,  "ctx1_ng_ml",
  "Confounder", "F",    "sos-enmo-conf",    "Tibial SOS",  "m/s",             1,  "tibia_sos"
)

# The reported contrast comes from the registry, not from a hard-coded string
# or a guess at the exposure name. A label that does not track the scaling the
# fits were built on mislabels every row by a factor of 5 or 10, with nothing
# raised.
exposure_unit <- function(spec) spec$contrast_label


# ---- Per-spec computation ----

spec_summary <- readr::read_csv(here::here("outputs", "tables", "spec-summary.csv"),
                                show_col_types = FALSE)

one_row <- function(dag, panel, spec_key, outcome_lab, unit, report_div, outcome_col) {
  spec <- model_templates[[spec_key]]
  ex   <- spec$exposure

  pred_draws  <- spec_draws(spec_key, "pred")
  slope_draws <- spec_draws(spec_key, "slope")

  # --- linear-projection slope, 95% HPDI only (natural units; no point estimate) ---
  lp   <- linear_projection(pred_draws, exposure = !!rlang::sym(ex), level = 0.95)
  slope_lo  <- lp$beta_hpdi$lo / report_div
  slope_hi  <- lp$beta_hpdi$hi / report_div
  p_neg     <- mean(lp$beta_draws$beta < 0)          # scale-invariant

  # --- linearity threshold: q95 of per-draw max |resid from linear fit|,
  #     as % of the outcome's observed range -- the stable,
  #     comparable-across-analyses scale ---
  max_dev <- pred_draws |>
    dplyr::filter(!is.na(draw)) |>
    dplyr::group_by(drawid) |>
    # max |deviation from the closed-form OLS projection| per draw (lm-free posterior functional)
    dplyr::summarize(md = {
                       x <- .data[[ex]]; xc <- x - mean(x)
                       b <- sum(xc * (draw - mean(draw))) / sum(xc^2)
                       a <- mean(draw) - b * mean(x); max(abs(draw - (a + b * x)))
                     },
                     .groups = "drop") |>
    dplyr::pull(md)
  outcome_range <- diff(range(dat[[outcome_col]], na.rm = TRUE))
  lin_thr_pct   <- unname(quantile(max_dev, 0.95)) / outcome_range * 100

  # --- flatness threshold: q95 of per-draw max |local slope|, outcome units ---
  max_abs <- slope_draws |>
    dplyr::filter(!is.na(draw)) |>
    dplyr::group_by(drawid) |>
    dplyr::summarize(ma = max(abs(draw)), .groups = "drop") |>
    dplyr::pull(ma)
  flat_thr <- unname(quantile(max_abs, 0.95)) / report_div

  n_a <- spec_summary$n_analytic[spec_summary$spec_key == spec_key]

  tibble::tibble(
    dag, panel, spec_key, outcome = outcome_lab,
    exposure   = ifelse(grepl("steps", ex), "Daily step count", "Mean daily ENMO"),
    n          = n_a,
    slope_unit = paste(unit, exposure_unit(spec)),
    hpdi_lo    = slope_lo,
    hpdi_hi    = slope_hi,
    p_beta_neg = p_neg,
    linearity_pct_of_outcome_range = lin_thr_pct,
    flatness_unit = paste(unit, exposure_unit(spec)),
    flatness   = flat_thr
  )
}

tab <- do.call(rbind, Map(one_row, rows$dag, rows$panel, rows$spec_key, rows$outcome_lab,
                          rows$unit, rows$report_div, rows$outcome_col))


# ---- Write CSV ----

out_dir <- here::here("outputs", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(out_dir, "supp-slope-table.csv")
readr::write_csv(tab, out_csv)
cat("Wrote", out_csv, "\n\n")


# ---- Emit a GitHub-flavored markdown table ----

# significant-figure formatting that matches the main-text prose precision
fmt <- function(x, sos_msstyle = FALSE) {
  ifelse(abs(x) >= 100, sprintf("%.0f", x),
  ifelse(abs(x) >= 1,   sprintf("%.1f", x),
  ifelse(abs(x) >= 0.01, sprintf("%.3f", x), sprintf("%.4f", x))))
}
sgn <- function(x) ifelse(x >= 0, paste0("+", fmt(x)), fmt(x))

# Compact layout: units stated once per row (Outcome carries the outcome unit,
# Exposure carries the per-exposure unit); Slope and Flatness are bare numbers
# in [outcome unit] per [exposure unit].
# Derived from `rows` rather than restated: a second copy of the unit mapping is
# how the markdown table came to disagree with the CSV's own `slope_unit`.
unit_of <- setNames(rows$unit, rows$outcome_lab)
md <- c(
  "| DAG | Panel | Outcome | Exposure | *n* | Slope (95% HPDI) | P(β < 0) | Linearity threshold | Flatness threshold |",
  "|:---:|:---:|---|---|---:|---|---:|---:|---:|"
)
for (i in seq_len(nrow(tab))) {
  r <- tab[i, ]
  md <- c(md, sprintf(
    "| %s | %s | %s (%s) | %s | %d | %s to %s | %.2f | %.0f%% | %s |",
    r$dag, r$panel, r$outcome, unit_of[[r$outcome]],
    exposure_unit(model_templates[[r$spec_key]]),
    r$n,
    sgn(r$hpdi_lo), sgn(r$hpdi_hi),
    r$p_beta_neg,
    round(r$linearity_pct_of_outcome_range),
    fmt(r$flatness)
  ))
}
cat(paste(md, collapse = "\n"), "\n")
