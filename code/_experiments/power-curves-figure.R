###############################################################################
# Assemble the nine a priori power curves into one supplement figure. Reads the
# two power CSVs (PA->bone in SD units; industrialization converted to SD units),
# plots power vs true effect size per analysis, and marks the landmarks:
#   - dashed vertical  = 0.2 SD (small standardized effect)
#   - solid vertical   = measurement-resolution floor (where defined)
#   - dotted horizontal = 0.80 power
# PA->bone panels (green) are well-powered; industrialization panels (red, village
# level) are not. Output: outputs/figures/final/supp-fig-5-power-curves.{pdf,png}
#
# Run: script -q /dev/null Rscript code/_experiments/power-curves-figure.R < /dev/null
###############################################################################

suppressMessages({ library(here); library(ggplot2); library(dplyr) })

pb <- read.csv(here("outputs", "_experiments", "power-curves", "pa-bone-power.csv"))
ub <- read.csv(here("outputs", "_experiments", "power-curves", "industrialization-power.csv"))

pa <- pb |>
  transmute(key, group = "Physical activity → bone (individual level)",
            effect_sd, power, resolution_sd)

ur <- ub |>
  mutate(sd = sd0.2 / 0.2,
         effect_sd     = effect / sd,
         resolution_sd = ifelse(is.na(resolution), NA_real_, resolution / sd)) |>
  transmute(key, group = "Industrialization → outcome (village level)",
            effect_sd, power, resolution_sd)

lab <- c("sos-steps" = "Steps → tibial SOS",     "sos-enmo" = "ENMO → tibial SOS",
         "ctx-steps" = "Steps → CTX-1",          "ctx-enmo" = "ENMO → CTX-1",
         "osteo-steps" = "Steps → osteocalcin",  "osteo-enmo" = "ENMO → osteocalcin",
         "sos-urb" = "Industrialization → tibial SOS",
         "steps-urb" = "Industrialization → steps",
         "enmo-urb" = "Industrialization → ENMO")

dat <- bind_rows(pa, ur) |>
  mutate(analysis = factor(unname(lab[key]), levels = unname(lab)))

lm_res <- dat |> group_by(analysis, group) |> summarize(res = resolution_sd[1], .groups = "drop")

cols <- c("Physical activity → bone (individual level)"    = "#3F7E4E",
          "Industrialization → outcome (village level)"     = "#A85D52")

p <- ggplot(dat, aes(effect_sd, power, color = group)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "grey40", linewidth = 0.8) +
  geom_vline(xintercept = 0.2, linetype = "22", color = "grey40", linewidth = 0.8) +
  geom_vline(data = subset(lm_res, !is.na(res)), aes(xintercept = res),
             color = "grey20", linewidth = 0.9) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.1) +
  facet_wrap(~ analysis, scales = "free_x", ncol = 3) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_color_manual(values = cols) +
  labs(x = "True effect size (outcome SD units, over the meaningful contrast: 5,000 steps / 10 mg ENMO, or across the industrialization gradient)",
       y = "Power", color = NULL,
       title = "A priori power curves (frequentist proxy, validated against the Bayesian estimator; N = 2,000 per point)",
       subtitle = "Solid vertical = measurement-resolution floor   |   dashed vertical = 0.2 SD   |   dotted horizontal = 0.80 power") +
  theme_bw(base_size = 9) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 7.5, color = "grey30"),
        strip.text = element_text(size = 7.5))

ggsave(here("outputs", "figures", "final", "supp-fig-5-power-curves.pdf"), p, width = 8.2, height = 7)
ggsave(here("outputs", "figures", "final", "supp-fig-5-power-curves.png"), p, width = 8.2, height = 7, dpi = 200)
cat("wrote outputs/figures/final/supp-fig-5-power-curves.{pdf,png}\n")
