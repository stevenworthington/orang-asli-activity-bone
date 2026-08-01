###############################################################################
# Linear-projection slopes for the plug-in / cut / joint industrialization
# estimators, so the propagation comparison is made on the quantity the
# manuscript actually reports.
#
# industrialization-village-two-stage-propagated.R compares the three estimators on the
# smooth's across-gradient ENDPOINT CONTRAST. The Results report the
# LINEAR-PROJECTION SLOPE of the AERF, per 10 index units, led by the posterior
# probability of decline. Those are different summaries of the same posterior and
# can disagree about whether an interval covers zero, so the reported quantity has
# to be recomputed for each method before any claim is made about whether full
# propagation changes a result.
#
# Read-only: consumes the saved AERF draws, prints numbers, edits nothing.
# Mirrors industrialization-linear-summaries.R exactly (same linear_projection helper).
#
# Run: Rscript code/_experiments/propagation-linear-summaries.R
###############################################################################

library(here)
source(here("code", "_startup", "init.R"))
suppressMessages({ library(dplyr) })

prop_dir <- here("outputs", "_experiments", "industrialization-village-two-stage-propagated")
out_rows <- list()

slope_summary <- function(ad, grid, ex, method, spec_key) {
  nd <- nrow(ad); ng <- ncol(ad)
  pred <- tibble::tibble(drawid = rep(seq_len(nd), times = ng), draw = as.vector(ad))
  pred[[ex]] <- rep(grid, each = nd)
  lp <- linear_projection(pred, exposure = !!rlang::sym(ex), level = 0.95)
  data.frame(spec = spec_key, method = method,
             slope10_lo = lp$beta_hpdi$lo * 10,          # per 10 index units
             slope10_hi = lp$beta_hpdi$hi * 10,
             p_decline  = mean(lp$beta_draws$beta < 0),
             excludes_0 = (lp$beta_hpdi$lo > 0) | (lp$beta_hpdi$hi < 0))
}

for (key in c("steps-urb", "enmo-urb", "sos-urb")) {
  spec <- model_templates[[key]]
  ex   <- spec$exposure
  p    <- readRDS(file.path(prop_dir, paste0("propagated_", key, ".rds")))
  s2   <- readRDS(here("outputs", "_experiments", "industrialization-village-two-stage",
                       paste0("stage2_", key, ".rds")))

  rows <- rbind(slope_summary(s2$aerf_draws, s2$grid, ex, "plugin", key),
                slope_summary(p$ep_cut,      p$grid,  ex, "cut",    key),
                slope_summary(p$ep_joint,    p$grid,  ex, "joint",  key))

  cat(sprintf("\n=== %s  (%s) ===\n", key, spec$outcome_label))
  cat(sprintf("  %-8s %14s %14s %12s %11s\n", "method", "slope/10 lo", "slope/10 hi",
              "P(decline)", "excl. zero"))
  for (i in seq_len(nrow(rows)))
    cat(sprintf("  %-8s %14.4g %14.4g %12.3f %11s\n", rows$method[i],
                rows$slope10_lo[i], rows$slope10_hi[i], rows$p_decline[i],
                ifelse(rows$excludes_0[i], "yes", "no")))
  out_rows[[key]] <- rows
}

tab <- do.call(rbind, out_rows)
write.csv(tab, file.path(prop_dir, "propagation-linear-slopes.csv"), row.names = FALSE)
cat("\nsaved", file.path(prop_dir, "propagation-linear-slopes.csv"), "\n")
