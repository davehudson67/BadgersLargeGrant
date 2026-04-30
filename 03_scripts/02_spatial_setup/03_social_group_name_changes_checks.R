library(nimble)
library(tidyverse)
library(MCMCpack)

bc <- readRDS("../BadgersSG2024.rds")
bc$socg <- as.factor(bc$socg)
levels(bc$socg)

bc <- bc %>%
  group_by(tattoo) %>%
  filter(pm != "Yes") %>%
  filter(!is.na(socg)) %>%
  mutate(count = n()) %>%
  filter(count > 2) %>%
  ungroup() %>%
  droplevels() %>%
  mutate(state = 1)

# Get unique individual IDs and time occasions
ids <- unique(bc$tattoo)
occs <- seq(1:(max(bc$occ)))

# Function to count SG changes for each individual
count_sg_changes <- function(sg) {
  return(sum(sg[-1] != sg[-length(sg)]))
}

## SG changes
# Summarize the DataFrame to count SG changes and captures for each individual
individual_summary <- bc %>%
  group_by(tattoo) %>%
  summarize(
    Captures = n(),
    SG_Changes = count_sg_changes(socg)
  )

# Calculate the rate of SG changes per capture for each individual
individual_summary <- individual_summary %>%
  mutate(SG_Change_Rate = SG_Changes / Captures)

# Merge the individual summary back to the original data to associate each individual with their rate
merged_df <- bc %>%
  left_join(individual_summary, by = "tattoo")

# Summarize the rate of SG changes for each Social Group
sg_summary <- merged_df %>%
  group_by(socg) %>%
  summarize(Average_SG_Change_Rate = mean(SG_Change_Rate, na.rm = TRUE))

# Create a bar plot to visualize the average rate of SG changes for each Social Group
ggplot(sg_summary, aes(x = socg, y = Average_SG_Change_Rate)) +
  geom_bar(stat = "identity") +
  labs(title = "Average Rate of SG Changes for Each Social Group",
       x = "Social Group",
       y = "Average Rate of SG Changes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
