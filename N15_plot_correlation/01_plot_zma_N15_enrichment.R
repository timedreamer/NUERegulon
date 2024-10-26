# This script plots the N15 enrichment for the maize time-course shoot and root data.

# Author: Ji Huang

# Output: zma_N15_enrichment_plot.pdf


# 0. Prep -----------------------------------------------------------------

library(here)
library(tidyverse)
library(readxl)
library(ggbeeswarm)
library(cowplot)
library(ggpubr)

# 1. Load data ------------------------------------------------------------

n15 <- read_excel(here("N15_result.xlsx"),
    skip = 1,
    col_names = c("sample_id", "N15", "sample_name", "time", "tissue")) %>%
    mutate(time_linear = str_remove(time, "m")) %>%
    mutate(time_linear = as.numeric(time_linear)) %>%
    mutate(time = factor(time,
        levels = c("0m", "5m", "10m", "15m", "20m", "30m", "45m", "60m", "90m", "120m")))

n15_shoot <- n15 %>%
    filter(tissue == "shoot")

n15_root <- n15 %>%
    filter(tissue == "root")

# 2. Plot -----------------------------------------------------------------

# Include a linear model
p3_lm <- ggplot(n15_shoot, aes(x = time_linear, y = N15)) +
    geom_beeswarm(color = "#009E73") +
    stat_smooth(method = "lm", se = FALSE, color = "#009E73") +
    scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 30, 45, 60, 90, 120)) +
    cowplot::theme_minimal_hgrid() +
    labs(title = expression("Maize Shoot "^15 * "N enrichment"),
        x = "Time", y = expression(""^15 * "N at-% (atomic percent)")) +
    stat_regline_equation(aes(label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~")), color = "black", size = 5)

p4_lm <- ggplot(n15_root, aes(x = time_linear, y = N15)) +
    geom_beeswarm(color = "#E69F00") +
    stat_smooth(method = "lm", se = FALSE, color = "#E69F00") +
    scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 30, 45, 60, 90, 120)) +
    cowplot::theme_minimal_hgrid() +
    labs(title = expression("Maize Root "^15 * "N enrichment"),
        x = "Time", y = expression(""^15 * "N at-% (atomic percent)")) +
    stat_regline_equation(aes(label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~")), color = "black", size = 5)

# 3. Save -----------------------------------------------------------------

combined_plot4  <- plot_grid(p3_lm, p4_lm, ncol = 1)


pdf(here("result", "N15", "zma_N15_enrichment_plot.pdf"),
    width = 5, height = 8)

print(combined_plot4)
dev.off()

