# Load necessary libraries
library(sf)
library(tidyverse)

# Read the shapefiles
grid <- st_read("Grid/GridShapeFile.shp")
studyArea <- st_read("StudyArea/StudyAreaShapeFile.shp")
#settLocations <- st_read("SettLocations/SettLocationsShapeFile.shp")
settGrid <- st_read("SettGrid/SettGridLocationsR.shp")

# Calculate the number of rows and columns
cols <- length(min(grid$col_index):max(grid$col_index))
rows <- length(min(grid$row_index):max(grid$row_index))

# Initialize a matrix to store the grid values
grid_matrix <- matrix(0, nrow = rows, ncol = cols)

studyArea$row_index <- studyArea$row_index + 1
studyArea$col_index <- studyArea$col_index + 1

# Mark the study area grid squares in the matrix
for (i in 1:nrow(studyArea)) {
  grid_matrix[studyArea$row_index[i], studyArea$col_index[i]] <- 1
}

# Adjust sett location indices to 1-based indexing for R
settGrid$row_index <- settGrid$row_index + 1
settGrid$col_index <- settGrid$col_index + 1

# Mark the sett grid squares in the matrix
for (i in 1:nrow(settGrid)) {
  grid_matrix[settGrid$row_index[i], settGrid$col_index[i]] <- 2
}

# Convert the matrix to a data frame for plotting
grid_df <- as.data.frame(as.table(grid_matrix))
colnames(grid_df) <- c("row", "col", "value")

# Ensure the row and col columns are numeric
grid_df$row <- as.numeric(grid_df$row)
grid_df$col <- as.numeric(grid_df$col)

grid_df$row <- 34 - grid_df$row

# Plot the grid matrix with study area and sett locations
ggplot(grid_df, aes(x = col, y = row, fill = factor(value))) +
  geom_tile() +
  scale_fill_manual(values = c("white", "grey", "red"), 
                    labels = c("Background", "Study Area", "Sett Locations"), 
                    name = "Legend") +
#  scale_y_reverse() +  # Reverse the y-axis
  coord_equal() +
  labs(title = "Study Area with Sett Locations", x = "X", y = "Y") +
  theme_minimal() +
  theme(legend.position = "right")


# Setup data for analysis
settGrid <- settGrid %>%
  select(Sett, row_index, col_index) %>%
  mutate(Sett = toupper(Sett)) %>%
  mutate(row_index = 34 - row_index)

studyArea <- studyArea %>%
  select(row_index, col_index) %>%
  mutate(row_index = 34 - row_index) %>%
  mutate(studyArea = 1)

rm(grid)

save.image("spatialInfo.RData")
