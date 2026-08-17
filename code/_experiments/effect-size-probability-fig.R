# Supplement figure: effect-size-vs-posterior-probability curves (the "reverse-ROPE").
# Per-panel effect-size metric + bound chosen by shape and result:
#   - DETECTED monotonic:   effect = |across-gradient change|;        report T95 (>=).
#   - NULL (PA -> bone):    effect = |proj slope x LOCKED contrast|;  report U95 (<=, rule out larger).
# PA-null contrasts are LOCKED: daily steps = 5,000 steps, mean daily ENMO = 10 mg.
# One unit of `ad_steps_5k` IS 5,000 steps and one unit of `enmo_10` IS 10 mg, so
# the contrast weight w is 1. The weights must move with the column names, or the
# reported effects are wrong by a factor of 5 and 10 with nothing raising an error.
# Each panel also marks the reconciled measurement-noise floor for its OUTCOME.
#
# ALL nine panels read the production targets cache, so there is one source for
# every panel and no panel can be rebuilt from a different set of fits than its
# neighbour.
# x = candidate effect size (outcome units); y = P(|effect| > x). Base R + qs2,
# plus the shared reference scales -- which are plain constants, so this script
# still runs without the startup set.
suppressMessages(library(qs2))
source(here::here("code", "_experiments", "_reference-scales.R"))

read_obj <- function(path) {
  for (fn in c("qs_read", "qd_read"))
    if (exists(fn, where = asNamespace("qs2"), inherits = FALSE)) {
      o <- tryCatch(getExportedValue("qs2", fn)(path), error = function(e) NULL)
      if (!is.null(o)) return(o)
    }
  stop("cannot read ", path)
}

# Measurement-resolution floor per OUTCOME, natural units, from
# `_reference-scales.R` (NA = no clean instrument floor, i.e. the ENMO outcome).
floor_of <- NOISE_FLOOR

specs <- list(
  list(key="sos.urb",    ex="industrial_index",       lab="SOS ~ industrialization",   unit="m/s",   type="detected", out="sos",   tag="community-level"),
  list(key="enmo.urb",   ex="industrial_index",       lab="ENMO ~ industrialization",  unit="mg",    type="detected", out="enmo",  tag="community-level"),
  list(key="steps.urb",  ex="industrial_index",       lab="Steps ~ industrialization", unit="steps", type="detected", out="steps", tag="community-level"),
  list(key="sos.steps",  ex="ad_steps_5k",            lab="SOS ~ daily steps",         unit="m/s",   type="null",     out="sos",   w=1,  wlab="per 5,000 steps"),
  list(key="ctx.steps",  ex="ad_steps_5k",            lab="CTX-1 ~ daily steps",       unit="ng/mL", type="null",     out="ctx",   w=1,  wlab="per 5,000 steps"),
  list(key="osteo.steps",ex="ad_steps_5k",            lab="Osteocalcin ~ daily steps", unit="pg/mL", type="null",     out="osteo", w=1,  wlab="per 5,000 steps"),
  list(key="sos.enmo",   ex="enmo_10",                lab="SOS ~ ENMO",                unit="m/s",   type="null",     out="sos",   w=1,  wlab="per 10 mg ENMO"),
  list(key="ctx.enmo",   ex="enmo_10",                lab="CTX-1 ~ ENMO",              unit="ng/mL", type="null",     out="ctx",   w=1,  wlab="per 10 mg ENMO"),
  list(key="osteo.enmo", ex="enmo_10",                lab="Osteocalcin ~ ENMO",        unit="pg/mL", type="null",     out="osteo", w=1,  wlab="per 10 mg ENMO")
)

green <- "#3F7E4E"; amber <- "#C77B30"; floorc <- "#2563EB"; dark <- "#1E293B"
fmtn <- function(x) if (abs(x) >= 100) sprintf("%.0f", x) else if (abs(x) >= 1) sprintf("%.1f", x) else sprintf("%.3g", x)

out_dir <- here::here("outputs", "figures", "final")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out_dir, "supp-fig-4-effect-size.pdf")
# Sized for a journal page. Margins are tight and thousands ticks are
# k-formatted so nine panels fit at a pointsize that stays readable in print.
pdf(out, width = 6.3, height = 6.72, pointsize = 11)
par(mfrow = c(3, 3), oma = c(0, 2.0, 2.6, 0), mgp = c(2.5, 0.6, 0))
for (i in seq_along(specs)) {
  s   <- specs[[i]]
  col <- (i - 1) %% 3 + 1                                          # 1 = left column (only one with y-axis numbers)
  par(mar = c(4.4, 1.8, 2.9, 0.5))                                 # uniform: every panel the same plot size
  # One source for every panel: the production targets cache.
  pred <- as.data.frame(read_obj(here::here("_targets", "objects", paste0("pred_draws_", s$key))))
  pred <- pred[is.finite(pred$draw), ]
  pred <- pred[order(pred$drawid, pred[[s$ex]]), ]
  x  <- sort(unique(pred[[s$ex]])); ng <- length(x); nd <- nrow(pred) / ng
  M  <- matrix(pred$draw, nrow = nd, ncol = ng, byrow = TRUE)

  if (s$type == "null") {
    xc <- x - mean(x)                                              # OLS proj slope per draw
    eff <- abs(as.numeric((M %*% xc) / sum(xc^2)) * s$w)           # |slope x locked contrast|
    xlab_metric <- paste0("Absolute effect (", s$unit, ")\n", s$wlab)
  } else {
    eff <- abs(M[, ng] - M[, 1])                                  # detected: |across-gradient change|
    xlab_metric <- paste0("Absolute across-gradient\nchange (", s$unit, ")")
  }

  ts_hi <- as.numeric(quantile(eff, 0.995))
  if (!is.na(floor_of[[s$out]])) ts_hi <- max(ts_hi, floor_of[[s$out]] * 1.1)  # keep the noise floor in frame
  ts  <- seq(0, ts_hi, length.out = 250)
  p   <- vapply(ts, function(t) mean(eff > t), numeric(1))

  plot(ts, p, type = "l", lwd = 2.0, col = dark, ylim = c(0, 1), xaxt = "n", yaxt = "n",
       xlab = xlab_metric, ylab = "", main = "",
       cex.lab = 0.9)
  at <- pretty(ts)
  axis(1, at = at, labels = ifelse(at >= 1000, paste0(at / 1000, "k"), format(at, trim = TRUE)),
       cex.axis = 0.85)
  if (col == 1) axis(2, cex.axis = 0.85)                          # left column only; mid/right share it
  title(main = s$lab, line = 1.5, cex.main = 0.95)                # raised to leave a gap above the subtitle tag
  grid(col = "grey92")
  abline(h = c(0.05, 0.95), lty = 2, lwd = 1.4, col = "grey70")

  fl <- floor_of[[s$out]]                                         # measurement-noise floor line
  if (!is.na(fl) && fl <= max(ts)) {
    abline(v = fl, lty = 1, col = floorc, lwd = 1.3)
    lab_adj <- if (fl > 0.5 * max(ts)) 1.05 else -0.05            # label on the clear side of the line
    y_floor <- if (mean(eff > fl) > 0.5) 0.15 else 0.55          # and clear of the curve vertically
    if (s$key == "ctx.enmo") y_floor <- y_floor - 0.04           # small nudge down, clear of the curve here
    text(fl, y_floor, sprintf("noise floor %s", fmtn(fl)), adj = lab_adj, cex = 0.65, col = floorc)
  }

  if (s$type == "null") {
    U95 <- as.numeric(quantile(eff, 0.95))
    abline(v = U95, lty = 1, col = amber, lwd = 2.0)
    mtext("null", side = 3, line = 0.1, adj = 0.02, cex = 0.65, col = amber, font = 3)
    mtext(sprintf("U95 <= %s %s", fmtn(U95), s$unit), side = 3, line = -2.4, adj = 0.97, cex = 0.68, col = amber)
  } else {
    T95 <- as.numeric(quantile(eff, 0.05))
    abline(v = T95, lty = 1, col = green, lwd = 2.0)
    tag <- if (!is.null(s$tag)) s$tag else "detected"
    mtext(tag, side = 3, line = 0.1, adj = 0.02, cex = 0.65, col = green, font = 3)
    mtext(sprintf("T95 >= %s %s", fmtn(T95), s$unit), side = 3, line = -2.4, adj = 0.97, cex = 0.68, col = green)
  }
}
mtext("Effect size vs posterior probability", outer = TRUE, cex = 1.0, font = 2, line = 1.2)
mtext("Green T95: 95% lower bound (industrialization)    Amber U95: 95% upper bound (PA->bone null)    Blue: noise floor",
      outer = TRUE, cex = 0.62, line = 0.1)
mtext("P( | effect | > x)", side = 2, outer = TRUE, line = 0.4, cex = 0.9)
dev.off()
cat("saved", out, "\n")
