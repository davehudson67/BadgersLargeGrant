library(sf)
library(dplyr)
library(ggplot2)

grid <- st_read("Habitat_Final.shp")

names(grid)
sum(is.na(grid$SG_id))
sum(is.na(grid$habitat))
st_crs(grid)

# Check lakes
ggplot(grid) +
  geom_sf(aes(fill = factor(habitat)), colour = NA) +
  coord_sf() +
  theme_minimal()

# Check SG boundaries
ggplot(grid) +
  geom_sf(aes(fill = factor(SG_id)), colour = NA) +
  coord_sf() +
  theme_minimal()

# Use cell size used in QGIS
cell_size <- 25  # meters

# Get centroids of the grid squares
grid_cent <- st_centroid(grid)
coords <- st_coordinates(grid_cent)

grid <- grid %>%
  mutate(
    x = coords[, 1],
    y = coords[, 2]
  )

# Get grid bounding box
bb <- st_bbox(grid)

xmin <- as.numeric(bb["xmin"])
xmax <- as.numeric(bb["xmax"])
ymin <- as.numeric(bb["ymin"])
ymax <- as.numeric(bb["ymax"])

# Convert coordinates to matrix indices
grid <- grid %>%
  mutate(
    col_id = as.integer(floor((x - xmin) / cell_size)) + 1,
    row_id = as.integer(floor((ymax - y) / cell_size)) + 1
  )

head(grid[, c("id", "x", "y", "col_id", "row_id", "SG_id", "habitat")])
n_rows <- max(grid$row_id)
n_cols <- max(grid$col_id)

n_rows
n_cols
n_rows * n_cols
nrow(grid)

habitat_mat <- matrix(NA_integer_, nrow = n_rows, ncol = n_cols)
SG_mat      <- matrix(NA_integer_, nrow = n_rows, ncol = n_cols)

habitat_mat[cbind(grid$row_id, grid$col_id)] <- as.integer(grid$habitat)
SG_mat[cbind(grid$row_id, grid$col_id)]      <- as.integer(grid$SG_id)

any(is.na(habitat_mat))
any(is.na(SG_mat))

all(habitat_mat[cbind(grid$row_id, grid$col_id)] == grid$habitat)
all(SG_mat[cbind(grid$row_id, grid$col_id)] == grid$SG_id)

landscape_inputs <- list(
  habitat_mat = habitat_mat,
  SG_mat = SG_mat,
  grid = grid,
  xmin = xmin,
  xmax = xmax,
  ymin = ymin,
  ymax = ymax,
  cell_size = cell_size,
  n_rows = n_rows,
  n_cols = n_cols
)

saveRDS(landscape_inputs, "landscape_inputs.rds")
#landscape_inputs <- readRDS("landscape_inputs.rds")
#habitat_mat <- landscape_inputs$habitat_mat
#SG_mat      <- landscape_inputs$SG_mat
#xmin        <- landscape_inputs$xmin
#xmax        <- landscape_inputs$xmax
#ymin        <- landscape_inputs$ymin
#ymax        <- landscape_inputs$ymax
#cell_size   <- landscape_inputs$cell_size
#n_rows      <- landscape_inputs$n_rows
#n_cols      <- landscape_inputs$n_cols

# Quick check
lookup_cell <- function(x, y, xmin, ymax, cell_size, n_rows, n_cols) {
  
  col_id <- floor((x - xmin) / cell_size) + 1
  row_id <- floor((ymax - y) / cell_size) + 1
  
  col_id <- pmax(1, pmin(n_cols, col_id))
  row_id <- pmax(1, pmin(n_rows, row_id))
  
  data.frame(
    x = x,
    y = y,
    row_id = row_id,
    col_id = col_id
  )
}

test <- grid[1:10, ]

test_lookup <- lookup_cell(
  x = test$x,
  y = test$y,
  xmin = xmin,
  ymax = ymax,
  cell_size = cell_size,
  n_rows = n_rows,
  n_cols = n_cols
)

test_lookup

habitat_mat[cbind(test_lookup$row_id, test_lookup$col_id)]
test$habitat

SG_mat[cbind(test_lookup$row_id, test_lookup$col_id)]
test$SG_id
