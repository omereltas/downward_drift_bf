# =============================================================================
# Downward Prior Drift — Figure and Table Generator
# Dosya : 02_analysis.R
# Calistir : source("02_analysis.R")
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ---- 0. Dizinler -------------------------------------------------------------
for (d in c("results", "figures", "tables")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ---- 1. -----------------------------------------------------------
find_rds <- function(filename) {
  candidates <- c(file.path("results", filename), filename)
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) stop(sprintf("'%s' bulunamadi.", filename))
  found[1]
}

er  <- readRDS(find_rds("error_rates.rds"))
sim <- readRDS(find_rds("all_simulations.rds"))
message(sprintf("Veri yuklendi: %d kosul, %d satir.", nrow(er), nrow(sim)))

# ---- 2. ------------------------------------------------
COLS <- c(
  "none"    = "#2C7BB6",
  "gradual" = "#D7191C",
  "sudden"  = "#FC8D59"
)

er <- er %>%
  mutate(
    drift_label = case_when(
      drift_magnitude == 0 ~ "none",
      TRUE                 ~ drift_type
    ),
    drift_label = factor(drift_label, levels = c("none", "gradual", "sudden")),
    prior_label = sprintf("r\u2080 = %.3f", prior_r0),
    delta_label = sprintf("\u03b4 = %.1f", true_delta),
    speed_label = ifelse(is.infinite(drift_speed), "speed = \u221e",
                         sprintf("speed = %g", drift_speed))
  )

# ---- 3. Figure 1: ----------------------------
fig1_data <- er %>%
  filter(true_delta == 0, !is.na(type_I_error)) %>%
  filter(drift_magnitude == 0 | drift_type %in% c("gradual", "sudden"))

p1 <- ggplot(fig1_data,
             aes(x        = drift_magnitude,
                 y        = type_I_error,
                 color    = drift_label,
                 linetype = factor(prior_r0),
                 group    = interaction(drift_label, prior_r0))) +
  geom_hline(yintercept = 0.05, linetype = "dashed",
             color = "gray50", linewidth = 0.6) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~ speed_label) +
  scale_color_manual(values = COLS,
                     name   = "Drift type",
                     labels = c("None (baseline)", "Gradual", "Sudden")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"),
                        name   = "Prior r\u2080") +
  scale_y_continuous(labels = percent_format(accuracy = 0.1),
                     limits = c(0, NA)) +
  scale_x_continuous(breaks = c(0, 0.1, 0.3, 0.6)) +
  labs(
    title    = "Type I Error Rate \u00d7 Downward Prior Drift Magnitude",
    subtitle = expression(H[0]~"true ("*delta*" = 0) | BF thresholds: 10 / 0.1 | "*italic(n)[max]*" = 500"),
    x        = expression("Drift magnitude ("*Delta*italic(r)*", downward)"),
    y        = "Type I error rate",
    caption  = sprintf("K = %d replications per condition", max(er$K))
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.box       = "vertical",
    strip.background = element_rect(fill = "gray95"),
    panel.grid.minor = element_blank()
  )

ggsave("figures/fig1_typeI_downward.pdf", p1, width = 8, height = 5.5)
ggsave("figures/fig1_typeI_downward.png", p1, width = 8, height = 5.5, dpi = 300)
message("Figure 1 kaydedildi.")

# ---- 4. Figure 2: -------------------------------------
fig2_data <- er %>%
  filter(true_delta > 0, !is.na(power)) %>%
  filter(drift_magnitude == 0 | drift_type %in% c("gradual", "sudden"))

p2 <- ggplot(fig2_data,
             aes(x        = drift_magnitude,
                 y        = power,
                 color    = drift_label,
                 linetype = factor(prior_r0),
                 group    = interaction(drift_label, prior_r0))) +
  geom_hline(yintercept = 0.80, linetype = "dashed",
             color = "gray50", linewidth = 0.6) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_grid(delta_label ~ speed_label) +
  scale_color_manual(values = COLS,
                     name   = "Drift type",
                     labels = c("None (baseline)", "Gradual", "Sudden")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"),
                        name   = "Prior r\u2080") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(
    title    = "Statistical Power \u00d7 Downward Prior Drift Magnitude",
    subtitle = expression(H[1]~"true | BF thresholds: 10 / 0.1 | "*italic(n)[max]*" = 500"),
    x        = expression("Drift magnitude ("*Delta*italic(r)*", downward)"),
    y        = "Power (1 \u2212 Type II error)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "gray95"),
    panel.grid.minor = element_blank()
  )

ggsave("figures/fig2_power_downward.pdf", p2, width = 9, height = 7)
ggsave("figures/fig2_power_downward.png", p2, width = 9, height = 7, dpi = 300)
message("Figure 2 kaydedildi.")


# ---- 6. Tablo 1: Hata oranlari (r0 = 0.707) ---------------------------------
tbl1 <- er %>%
  filter(prior_r0 == 0.707) %>%
  mutate(across(c(type_I_error, type_II_error, power, inconclusive),
                ~ round(.x, 4))) %>%
  select(condition_id, true_delta, drift_type, drift_magnitude,
         drift_speed, type_I_error, type_II_error, power,
         inconclusive, median_n, mean_n) %>%
  arrange(true_delta, drift_magnitude, drift_speed)

write.csv(tbl1, "tables/table1_error_rates_JZS_downward.csv", row.names = FALSE)
message("Tablo 1 -> tables/table1_error_rates_JZS_downward.csv")

# ---- 7. Tablo 2: Tum prior'lar ----------------------------------------------
tbl2 <- er %>%
  mutate(across(c(type_I_error, type_II_error, power, inconclusive),
                ~ round(.x, 4))) %>%
  select(condition_id, true_delta, prior_r0, drift_type,
         drift_magnitude, drift_speed,
         type_I_error, type_II_error, power, inconclusive,
         median_n, mean_n) %>%
  arrange(prior_r0, true_delta, drift_magnitude, drift_speed)

write.csv(tbl2, "tables/table2_all_error_rates_downward.csv", row.names = FALSE)
message("Tablo 2 -> tables/table2_all_error_rates_downward.csv")

message("=== 02_analysis.R tamamlandi ===")
