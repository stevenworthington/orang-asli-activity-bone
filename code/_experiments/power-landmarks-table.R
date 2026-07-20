###############################################################################
# Power-landmarks table — companion to Supp Fig 8 (the a priori power curves).
# Reads the two power-sweep CSVs that drive the figure and extracts, per analysis,
# the power at the two landmarks marked on each curve: the measurement-resolution
# floor and a small standardized effect of 0.2 SD. No re-simulation — pure readout:
# the two sweep scripts now include the landmark effect sizes as explicit grid points,
# so each landmark power is read directly (the interp() below is an exact lookup that
# only falls back to linear interpolation if a landmark were ever off-grid).
#
# Outputs:
#   outputs/tables/power-landmarks.csv  (machine-readable, one row per analysis)
#   stdout: a GitHub-markdown table for pasting into Supplementary Material 10.
#
# Run: script -q /dev/null Rscript code/_experiments/power-landmarks-table.R < /dev/null
###############################################################################

suppressMessages(library(here))

pb <- read.csv(here("outputs", "_experiments", "power-curves", "pa-bone-power.csv"))
ub <- read.csv(here("outputs", "_experiments", "power-curves", "industrialization-power.csv"))

# village counts (level of inference for the village-level analyses), from the
# Option H stage-2 summaries (same real design as the power simulation).
n_vill <- c("sos-urb" = 25L, "steps-urb" = 29L, "enmo-urb" = 29L)

# interpolate y at x0 from a monotone-x grid (NA if x0 is NA / out of range edges held flat)
interp <- function(x, y, x0) {
  if (is.na(x0)) return(NA_real_)
  o <- order(x); x <- x[o]; y <- y[o]
  as.numeric(approx(x, y, xout = x0, rule = 2)$y)   # rule=2: clamp at the ends
}


# ---- physical activity -> bone (individual level) ----------------------------

pb_lab <- c("sos-steps" = "Steps → tibial SOS",   "sos-enmo" = "ENMO → tibial SOS",
            "ctx-steps" = "Steps → CTX-1",         "ctx-enmo" = "ENMO → CTX-1",
            "osteo-steps" = "Steps → osteocalcin", "osteo-enmo" = "ENMO → osteocalcin")

pa <- do.call(rbind, lapply(names(pb_lab), function(k) {
  d <- pb[pb$key == k, ]
  res_sd <- d$resolution_sd[1]
  data.frame(
    key = k, label = pb_lab[[k]], level = "Physical activity → bone (individual level)",
    n = d$n[1], n_villages = d$n_villages[1], resolution_sd = res_sd,
    power_floor  = interp(d$effect_sd, d$power, res_sd),
    power_0.2sd  = interp(d$effect_sd, d$power, 0.2),
    stringsAsFactors = FALSE)
}))


# ---- industrialization -> outcome (village level) ----------------------------

ub_lab <- c("sos-urb" = "Industrialization → tibial SOS",
            "steps-urb" = "Industrialization → daily steps",
            "enmo-urb" = "Industrialization → mean ENMO")

ur <- do.call(rbind, lapply(names(ub_lab), function(k) {
  d <- ub[ub$key == k, ]
  sd_y      <- d$sd0.2[1] / 0.2                 # outcome between-person SD (natural units)
  eff_0.2sd <- d$sd0.2[1]                       # natural-unit effect equal to 0.2 SD
  res_nat   <- d$resolution[1]                  # NA for the two activity outcomes
  data.frame(
    key = k, label = ub_lab[[k]], level = "Industrialization → outcome (village level)",
    n = NA_integer_, n_villages = n_vill[[k]],
    resolution_sd = if (is.na(res_nat)) NA_real_ else res_nat / sd_y,
    power_floor   = interp(d$effect, d$power, res_nat),
    power_0.2sd   = interp(d$effect, d$power, eff_0.2sd),
    stringsAsFactors = FALSE)
}))

tab <- rbind(pa, ur)


# ---- write machine-readable CSV ----------------------------------------------

out_csv <- here("outputs", "tables", "power-landmarks.csv")
write.csv(tab, out_csv, row.names = FALSE)
cat("wrote", out_csv, "\n\n")


# ---- emit a markdown table for the supplement --------------------------------

fmtp  <- function(p) ifelse(is.na(p), "—", sprintf("%.2f", p))
fmtsd <- function(s) ifelse(is.na(s), "—", sprintf("%.2f SD", s))
fmtN  <- function(r) if (grepl("individual", r$level)) sprintf("%d", r$n) else sprintf("%d villages", r$n_villages)

cat("| Analysis | N | Resolution floor | Power at floor | Power at 0.2 SD |\n")
cat("|---|---|---|---|---|\n")
for (lv in unique(tab$level)) {
  sub <- tab[tab$level == lv, ]
  cat(sprintf("| **%s** | | | | |\n", lv))
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    cat(sprintf("| %s | %s | %s | %s | %s |\n",
                r$label, fmtN(r), fmtsd(r$resolution_sd), fmtp(r$power_floor), fmtp(r$power_0.2sd)))
  }
}
