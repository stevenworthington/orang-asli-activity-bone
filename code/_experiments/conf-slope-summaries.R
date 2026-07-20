###############################################################################
# Confounder-DAG (SM4) linear-projection summaries for the six physical-activity
# -> bone *-conf specs. Uses the SAME functionals as the main-text reporting:
#   - per-draw OLS slope of the AERF on exposure (matches linear_projection /
#     code/_final/supp-slope-table.R) -> P(slope < 0)
#   - per-person contrast effect = slope * contrast width (5 = 5,000 steps,
#     10 = 10 mg ENMO); U95 = 95th percentile of |effect| (matches
#     code/_experiments/fe-set1-experiment.R, which produced the main-text
#     "95% upper bound" magnitudes), reported as % of the outcome SD and as a
#     multiple of the assay-noise floor.
# Posterior probability of the HYPOTHESIZED direction:
#   - osteocalcin (formation) and tibial SOS (density): hypothesis positive,
#     P(hyp) = 1 - P(slope < 0)
#   - CTX-1 (resorption): hypothesis negative, P(hyp) = P(slope < 0)
###############################################################################


library(here)
source(here("code", "_startup", "init.R"))

hpdi <- function(x, m = 0.95) { x <- sort(x); n <- length(x); k <- floor(m * n)
  i <- which.min(x[(k + 1):n] - x[1:(n - k)]); c(x[i], x[i + k]) }

# reconciled measurement-noise floors + outcome SDs (natural units), from
# fe-set1-experiment.R (the source of the main-text magnitude numbers)
sds    <- c(sos = 182.15, ctx = 0.116, osteo = 19082.9)
floors <- c(sos = 31.6,   ctx = 0.008, osteo = 990)
units  <- c(sos = "m/s",  ctx = "ng/mL", osteo = "pg/mL")
cwidth <- function(ex) if (identical(ex, "ad_steps_1k")) 5 else 10

summ <- function(pred, ex, w) {
  pred <- pred[is.finite(pred$draw), ]
  pred <- pred[order(pred$drawid, pred[[ex]]), ]
  x <- sort(unique(pred[[ex]])); ng <- length(x); nd <- nrow(pred) / ng
  M <- matrix(pred$draw, nrow = nd, ncol = ng, byrow = TRUE)
  xc <- x - mean(x); slope <- as.numeric((M %*% xc) / sum(xc^2)); eff <- slope * w
  list(sl = hpdi(slope), pneg = mean(slope < 0), ch = hpdi(eff),
       u95 = as.numeric(quantile(abs(eff), 0.95)))
}

specs <- tibble::tribble(
  ~key,               ~hyp_pos,
  "osteo-steps-conf", TRUE,
  "ctx-steps-conf",   FALSE,
  "sos-steps-conf",   TRUE,
  "osteo-enmo-conf",  TRUE,
  "ctx-enmo-conf",    FALSE,
  "sos-enmo-conf",    TRUE
)

for (i in seq_len(nrow(specs))) {
  key <- specs$key[i]; hyp_pos <- specs$hyp_pos[i]
  spec <- model_templates[[key]]; ex <- spec$exposure; w <- cwidth(ex)
  oc <- sub("-.*", "", key); sd <- sds[[oc]]; fl <- floors[[oc]]; u <- units[[oc]]
  pred <- as.data.frame(qs2::qs_read(
    file.path("_targets", "objects", paste0("pred_draws_", gsub("-", ".", key)))))
  s <- summ(pred, ex, w)
  p_hyp <- if (hyp_pos) 1 - s$pneg else s$pneg
  cat(sprintf("\n=== %-16s (%s, contrast %s) ===\n", key, u,
              if (w == 5) "5,000 steps" else "10 mg ENMO"))
  cat(sprintf("  P(hypothesized direction) = %.3f   [P(slope<0)=%.3f]\n", p_hyp, s$pneg))
  cat(sprintf("  contrast effect 95%% HPDI  = [%.4g, %.4g] %s\n", s$ch[1], s$ch[2], u))
  cat(sprintf("  U95 (95th pct |effect|)   = %.4g %s = %.0f%% of SD = %.1fx assay-noise floor\n",
              s$u95, u, 100 * s$u95 / sd, s$u95 / fl))
}
