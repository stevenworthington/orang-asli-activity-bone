# Per-person contrast effects for the 6 PA -> bone nulls, on LOCKED contrasts:
#   daily steps: 5,000 steps/day (= 5 units of ad_steps_1k)
#   mean daily ENMO: 10 mg       (= 10 units of ad_mean_enmo_mg_0_24hr)
# Effect = best-linear-projection slope of the AERF (OLS of AERF on exposure, per draw,
# matching the figure caption) x contrast width. Reports signed 95% HPDI, P(direction),
# and U95 = 95th pctile of |effect| (the "rule out effects larger than" bound), then
# compares U95 to the reconciled measurement-noise floor and the cohort SD. Base R + qs2.
suppressMessages(library(qs2))
read_obj <- function(p) { for (fn in c("qs_read","qd_read")) if (exists(fn, where=asNamespace("qs2"), inherits=FALSE)) { o<-tryCatch(getExportedValue("qs2",fn)(p), error=function(e) NULL); if(!is.null(o)) return(o) }; stop("cannot read ",p) }
hpdi <- function(x, m=0.95) { x<-sort(x); n<-length(x); k<-floor(m*n); i<-which.min(x[(k+1):n]-x[1:(n-k)]); c(x[i], x[i+k]) }

# reconciled floors (2026-06-14) and cohort SDs, in the outcome's natural units
floor_sd <- list(
  sos   = list(unit="m/s",   floor=31.6,  sd=182.2,    name="Tibial SOS"),
  ctx   = list(unit="ng/mL", floor=0.008, sd=0.116,    name="CTX-1"),
  osteo = list(unit="pg/mL", floor=990,   sd=19082.9,  name="Osteocalcin")
)
specs <- list(
  list(key="sos.steps",  out="sos",   ex="ad_steps_1k",            w=5,  wlab="5,000 steps", hyp="+"),
  list(key="ctx.steps",  out="ctx",   ex="ad_steps_1k",            w=5,  wlab="5,000 steps", hyp="-"),
  list(key="osteo.steps",out="osteo", ex="ad_steps_1k",            w=5,  wlab="5,000 steps", hyp="+"),
  list(key="sos.enmo",   out="sos",   ex="ad_mean_enmo_mg_0_24hr", w=10, wlab="10 mg ENMO",  hyp="+"),
  list(key="ctx.enmo",   out="ctx",   ex="ad_mean_enmo_mg_0_24hr", w=10, wlab="10 mg ENMO",  hyp="-"),
  list(key="osteo.enmo", out="osteo", ex="ad_mean_enmo_mg_0_24hr", w=10, wlab="10 mg ENMO",  hyp="+")
)

cat(sprintf("%-13s %-12s | AERF mean | proj slope/unit (HPDI) | contrast HPDI | P(hyp) | U95 |x| | vs floor | vs SD\n", "spec","contrast"))
for (s in specs) {
  pred <- as.data.frame(read_obj(file.path("_targets","objects",paste0("pred_draws_",s$key))))
  pred <- pred[is.finite(pred$draw), ]
  pred <- pred[order(pred$drawid, pred[[s$ex]]), ]
  x  <- sort(unique(pred[[s$ex]])); ng <- length(x); nd <- nrow(pred)/ng
  M  <- matrix(pred$draw, nrow=nd, ncol=ng, byrow=TRUE)        # draws x grid (AERF)
  xc <- x - mean(x); sxx <- sum(xc^2)
  slope <- as.numeric((M %*% xc) / sxx)                        # OLS proj slope per draw
  eff   <- slope * s$w                                         # contrast effect per draw
  fs <- floor_sd[[s$out]]
  sl_h <- hpdi(slope); ef_h <- hpdi(eff)
  p_hyp <- if (s$hyp=="+") mean(eff>0) else mean(eff<0)
  u95   <- as.numeric(quantile(abs(eff), 0.95))
  cat(sprintf("%-13s %-12s | %9.3g | %.4g [%.4g, %.4g] | [%.4g, %.4g] | %.2f | %.4g %s | %.2fx | %.0f%%SD\n",
      s$key, s$wlab, mean(M), median(slope), sl_h[1], sl_h[2], ef_h[1], ef_h[2],
      p_hyp, u95, fs$unit, u95/fs$floor, 100*u95/fs$sd))
}
