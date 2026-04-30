library(sf)
library(sp)
library(raster)
library(tidyverse)

# Read the full grid and selected grid shapefiles
grid <- st_read("Grid/GridShapeFile.shp")
studyArea <- st_read("StudyArea/StudyAreaShapeFile.shp")
settLocations <- st_read("SettLocations/SettLocationsShapeFile.shp")

# Extract centroids for both grids
grid_centroids <- st_centroid(grid)
studyArea_centroids <- st_centroid(studyArea)

grid_coords <- st_coordinates(grid_centroids)
studyArea_coords <- st_coordinates(studyArea_centroids)

# Sort sett locations
sett_coords <- st_coordinates(settLocations)
settLocations <- settLocations %>%
  mutate(longitude = sett_coords[,1], latitude = sett_coords[,2]) %>%
  mutate(SettName = toupper(SettName)) %>%
  dplyr::select(SettName, longitude, latitude)

# Define the extent and grid size
x_min <- min(grid_coords[,1])
x_max <- max(grid_coords[,1])
y_min <- min(grid_coords[,2])
y_max <- max(grid_coords[,2])
grid_size <- 100

# Calculate the number of rows and columns
cols <- ceiling((x_max - x_min) / grid_size)
rows <- ceiling((y_max - y_min) / grid_size)

# Initialize a matrix to store the grid values
grid_matrix <- matrix(0, nrow = rows, ncol = cols)

# Function to convert coordinates to matrix indices
coord_to_index <- function(x, y, x_min, y_min, grid_size) {
  col <- ceiling((x - x_min) / grid_size)
  row <- ceiling((y - y_min) / grid_size)
  return(c(row, col))
}

# Mark the studyArea grid squares
for (i in 1:nrow(studyArea_coords)) {
  indices <- coord_to_index(studyArea_coords[i, 1], studyArea_coords[i, 2], x_min, y_min, grid_size)
  grid_matrix[indices[1], indices[2]] <- 1
}

# Convert the matrix to a data frame for plotting
grid_df <- as.data.frame(as.table(grid_matrix))

# Plot the grid matrix
ggplot(grid_df, aes(Var2, Var1, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  coord_equal() +
  labs(title = "Selected Grid Squares", x = "X", y = "Y") +
  theme_minimal()

save.image("SpatialInfo.RData")
