library(nimble)
library(tidyverse)
library(MCMCpack)
library(dplyr)

set.seed(11)
bc <- readRDS("../BadgersSG2024.rds")
names(bc)
bc <- filter(bc, !is.na(Cult_Sum), .preserve = T)

#------------------------------------------------------------------------------#
#
# Occasion
#
#------------------------------------------------------------------------------#

# Reduce to one row per badger per year
bc_unique <- bc %>%
  group_by(tattoo, occ) %>%
  summarise(
    culture_pos = as.integer(any(Cult_Sum > 0)),
    .groups = "drop"
  )

# Summarise by occ
yearly_summary <- bc_unique %>%
  group_by(occ) %>%
  summarise(
    n_tested = n(),
    n_culture_pos = sum(culture_pos),
    prop_culture_pos = n_culture_pos / n_tested,
    .groups = "drop"
  ) %>%
  mutate(
    sensitivity = 0.11,
    estimated_prevalence = pmin(prop_culture_pos / sensitivity, 1)
  )

# Create long-format plot data
plot_data <- yearly_summary %>%
  dplyr::select(occ, prop_culture_pos, estimated_prevalence, n_tested) %>%
  pivot_longer(
    cols = c(prop_culture_pos, estimated_prevalence),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
                    prop_culture_pos = "Culture Positivity Rate",
                    estimated_prevalence = "Estimated Prevalence (11% sensitivity)"
    )
  )

# Plot with points and sample size labels only for the culture line
ggplot(plot_data, aes(x = occ, y = value, color = metric, linetype = metric)) +
  geom_line(size = 1.2) +
  geom_point(
    data = filter(plot_data, metric == "Culture Positivity Rate"),
    size = 2) +
  geom_text(
    data = filter(plot_data, metric == "Culture Positivity Rate"),
    aes(label = n_tested),
    vjust = -1,
    size = 3) +
    labs(
    x = "Trapping occasion (1/4 years)",
    y = "Proportion / Estimated Prevalence",
    title = "Badger TB: Culture Positivity and Estimated Prevalence Over Time",
    color = "Metric",
    linetype = "Metric",
    caption = "Numbers above points = number of badgers tested") +
  theme_bw() +
  theme(legend.position = "bottom")

#------------------------------------------------------------------------------#

# Reduce to one row per badger per year
bc_unique <- bc %>%
  group_by(tattoo, capture_yr) %>%
  summarise(
    culture_pos = as.integer(any(Cult_Sum > 0)),
    .groups = "drop"
  )

# Summarise by occ
yearly_summary <- bc_unique %>%
  group_by(capture_yr) %>%
  summarise(
    n_tested = n(),
    n_culture_pos = sum(culture_pos),
    prop_culture_pos = n_culture_pos / n_tested,
    .groups = "drop"
  ) %>%
  mutate(
    sensitivity = 0.11,
    estimated_prevalence = pmin(prop_culture_pos / sensitivity, 1)
  )

# Create long-format plot data
plot_data <- yearly_summary %>%
  dplyr::select(capture_yr, prop_culture_pos, estimated_prevalence, n_tested) %>%
  pivot_longer(
    cols = c(prop_culture_pos, estimated_prevalence),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
                    prop_culture_pos = "Culture Positivity Rate",
                    estimated_prevalence = "Estimated Prevalence (11% sensitivity)"
    )
  )

# Plot with points and sample size labels only for the culture line
ggplot(plot_data, aes(x = capture_yr, y = value, color = metric, linetype = metric)) +
  geom_line(size = 1.2) +
  geom_point(
    data = filter(plot_data, metric == "Culture Positivity Rate"),
    size = 2) +
  geom_text(
    data = filter(plot_data, metric == "Culture Positivity Rate"),
    aes(label = n_tested),
    vjust = -1,
    size = 3) +
  labs(
    x = "Year",
    y = "Proportion / Estimated Prevalence",
    title = "Badger TB: Culture Positivity and Estimated Prevalence Over Time",
    color = "Metric",
    linetype = "Metric",
    caption = "Numbers above points = number of badgers tested") +
  theme_bw() +
  theme(legend.position = "bottom")

