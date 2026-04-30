library(tidyverse)
library(gganimate)
library(viridis)  # For color scales
library(magick)
library(coda)

run <- readRDS("EG_Badgers_InfOut_adjS.rds")
load("spatialInfo.RData")

samples <- as.data.frame(as.matrix(run$samples2)) %>%
  select(starts_with("S"))

# Extract indices from column names
col_names <- colnames(samples)
pattern <- "S\\[(\\d+), ?(\\d+), ?(\\d+)\\]"
matches <- stringr::str_match(col_names, pattern)

# Convert matches to a data frame
index_df <- data.frame(
  individual = as.numeric(matches[,2]),
  coordinate = as.numeric(matches[,3]),
  time = as.numeric(matches[,4]),
  col_name = col_names
)

# Adjust time
index_df$time <- index_df$time + 1981

# Melt the samples data frame to long format
samples_long <- samples %>%
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


# Create a data frame of infection status mapped to individuals
inf_df <- data.frame(
  individual = 1:length(inf),
  infection_status = ifelse(inf == 1, "positive as cub", "never positive")  # Convert 1/0 to "infected"/"not_infected"
)

# Ensure that the infection status is added to the mean_activity_centers
mean_activity_centers <- mean_activity_centers %>%
  left_join(inf_df, by = "individual")

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

# Create the plot with individuals colored by infection status
p <- ggplot(mean_activity_centers, aes(x = S_x, y = S_y, group = individual, color = infection_status)) +
  geom_tile(data = grid_df, aes(x = col, y = row, fill = as.factor(value)), inherit.aes = FALSE) +
  scale_fill_manual(values = c("0" = "white", "1" = "lightgreen"), guide = "none") +  # Custom colors for the study area
  geom_point(size = 2) +
  geom_path(alpha = 0.7) +
  scale_color_manual(values = c("positive as cub" = "red", "never positive" = "blue")) +
  labs(x = "X Coordinate", y = "Y Coordinate", color = "Infection Status") +
  theme_minimal() +
  ggtitle("Movement of Individuals Over Time")

# Add animation over time
p_anim <- p +
  transition_time(time) +
  labs(title = "Time Step: {frame_time}")  # Dynamic title showing the time step

# Render the animation
animation <- animate(p_anim, nframes = 50, fps = 10, width = 800, height = 800)

# Save the animation as a GIF
image_write(animation, path = "movement_infSadj.gif")
