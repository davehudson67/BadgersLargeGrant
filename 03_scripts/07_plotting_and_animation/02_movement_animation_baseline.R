library(tidyverse)
library(gganimate)
library(viridis)  # For color scales
library(magick)

run <- readRDS("EG_Badgers_NoCV_200indSummary.rds")
load("spatialInfo.RData")

samples <- as.data.frame(t(as.matrix(run$all.chains)))

samplesS <- samples %>%
  select(starts_with("S"))

# Assuming your df is called `samples` and column names follow the pattern S[i,1,k] or S[i,2,k]
# Let's start by extracting individual (i), coordinate (1/2), and time (k) from the column names

# Step 1: Extract indices from column names
col_names <- colnames(samplesS)
pattern <- "S\\[(\\d+), ?(\\d+), ?(\\d+)\\]"
matches <- stringr::str_match(col_names, pattern)

# Convert matches to a data frame
index_df <- data.frame(
  individual = as.numeric(matches[,2]),
  coordinate = as.numeric(matches[,3]),
  time = as.numeric(matches[,4]),
  col_name = col_names
)

# Melt the samples data frame to long format
samples_long <- samplesS %>%
  pivot_longer(cols = everything(), names_to = "col_name", values_to = "value")

# Merge the melted data with the index information
samples_reformatted <- left_join(samples_long, index_df, by = "col_name")

# Calculate the mean for each individual, coordinate, and time step
mean_activity_centers <- samples_reformatted %>%
  group_by(individual, coordinate, time) %>%
  summarize(mean_value = mean(value), .groups = 'drop')

mean_activity_centers <- mean_activity_centers %>%
  pivot_wider(names_from = coordinate, values_from = mean_value, names_prefix = "S_") %>%
  rename(S_x = S_1, S_y = S_2)

mean_activity_centers <- mean_activity_centers %>%
  filter(!is.na(S_x) & !is.na(S_y))

#p <- ggplot(mean_activity_centers, aes(x = S_x, y = S_y, group = individual)) +
#  geom_point(size = 2) +
#  geom_path(alpha = 0.7) +  # Draw the movement path
#  scale_color_viridis_d() +  # Use a colorblind-friendly color scale
#  labs(x = "X Coordinate", y = "Y Coordinate") +
#  theme_minimal() +
#  ggtitle("Movement of Individuals Over Time")

# Add animation over time
#p_anim <- p +
#  transition_time(time) +
#  labs(title = "Time Step: {frame_time}")  # Dynamic title showing the time step

# Render the animation
#animation <- animate(p_anim, nframes = 100, fps = 10, width = 800, height = 800)

## adding study area:
p <- ggplot(mean_activity_centers, aes(x = S_x, y = S_y, group = individual)) +
  geom_tile(data = grid_df, aes(x = col, y = row, fill = as.factor(value)), inherit.aes = FALSE) +
  scale_fill_manual(values = c("0" = "white", "1" = "lightgreen"), guide = "none") +  # Custom colors for the study area
  geom_point(size = 2) +
  geom_path(alpha = 0.7) +  # Draw the movement path
  scale_color_viridis_d() +  # Use a colorblind-friendly color scale
  labs(x = "X Coordinate", y = "Y Coordinate", color = "Individual") +
  theme_minimal() +
  ggtitle("Movement of Individuals Over Time")

# Add animation over time
p_anim <- p +
  transition_time(time) +
  labs(title = "Time Step: {frame_time}")  # Dynamic title showing the time step

# Render the animation
animation <- animate(p_anim, nframes = 50, fps = 10, width = 800, height = 800)

# Save the animation as a GIF
image_write(animation, path = "movement.gif")
