# ============================================================================
# Badger spatial CMR model: infection-history groups, adjusted first location
# MASTER CANDIDATE SCRIPT
# ============================================================================
#
# Purpose
# -------
# Fits an Ergon & Gardner-style spatial CMR model to badger capture data.
#
# Biological comparison:
#
#   group 1 = never positive
#   group 2 = test-positive as cub
#
# Main movement parameter:
#
#   dmean[1] = mean annual activity-centre movement for never-positive badgers
#   dmean[2] = mean annual activity-centre movement for test-positive-as-cub badgers
#
# Main contrast:
#
#   dmean[2] - dmean[1]
#
# Notes
# -----
# 1. First activity centre is fixed to the first observed social-group/sett
#    location for each individual.
#
# 2. The script currently retains the original `primary = as.numeric(primary) + 1L`
#    indexing. This means the first real year is usually primary index 2, and
#    primary index 1 is effectively unused. Do not change this until the model
#    is otherwise stable.
#
# 3. NIMBLE may still emit warnings such as:
#
#      warning: value in right-hand-side-only variable is NA or NaN, in variable: G
#
#    If the post-run checks show no NA/NaN/Inf values in the main monitored
#    parameters, this is probably coming from undefined unused array elements,
#    rather than from the core posterior samples. Keep documenting it, but do
#    not necessarily treat it as fatal if all validation checks pass.
#
# ============================================================================


# ---- 0. User options --------------------------------------------------------
TEST_RUN <- FALSE          # TRUE = short MCMC for debugging
SAMPLE_N <- NA           # use NA_integer_ for full dataset
RANDOM_SEED <- 42

# Monitor latent S/z?
#
# TRUE  = monitor active-window S/z nodes in run$samples2 for animation later.
# FALSE = only monitor main parameters. Use this if latent monitoring keeps
#         producing nuisance NA values while debugging the core model.
MONITOR_LATENT <- TRUE

# MCMC settings
MCMC_NITER_FULL <- 50000
MCMC_NBURN_FULL <- 12000
MCMC_NCHAINS_FULL <- 2

MCMC_NITER_TEST <- 1000
MCMC_NBURN_TEST <- 200
MCMC_NCHAINS_TEST <- 2
# ---- 1. Packages ------------------------------------------------------------
library(sf)
library(nimble)
library(lubridate)
library(tidyverse)
library(coda)
#library(mcmcplots)


# ---- 2. Helper functions ----------------------------------------------------
find_project_root <- function(start_dir = getwd()) {
  current <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  
  repeat {
    looks_like_root <- file.exists(file.path(current, "README_project_summary.md")) &&
      dir.exists(file.path(current, "02_data")) &&
      dir.exists(file.path(current, "03_scripts"))
    
    if (looks_like_root) return(current)
    
    parent <- dirname(current)
    
    if (identical(parent, current)) {
      stop(
        "Could not find the NewBadgers_organised project root.\n",
        "Open the NewBadgers_organised folder as your RStudio project, or set ",
        "`project_root` manually near the top of this script.",
        call. = FALSE
      )
    }
    
    current <- parent
  }
}


stop_if_missing_cols <- function(dat, cols, object_name = deparse(substitute(dat))) {
  missing <- setdiff(cols, names(dat))
  
  if (length(missing) > 0) {
    stop(
      "The object `", object_name, "` is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}


replace_na_if_exists <- function(dat, cols, value = 0) {
  # Replaces NAs in diagnostic columns that exist; creates absent columns as 0.
  for (nm in cols) {
    if (!nm %in% names(dat)) {
      dat[[nm]] <- value
    } else {
      dat[[nm]][is.na(dat[[nm]])] <- value
    }
  }
  
  dat
}


check_range <- function(mat, pattern, lower = -Inf, upper = Inf, label = pattern) {
  cols <- grep(pattern, colnames(mat), value = TRUE)
  
  if (length(cols) == 0) {
    warning("No columns found for ", label)
    return(invisible(NULL))
  }
  
  vals <- mat[, cols, drop = FALSE]
  
  n_low <- sum(vals < lower, na.rm = TRUE)
  n_high <- sum(vals > upper, na.rm = TRUE)
  
  cat(label, "\n")
  cat("  columns: ", length(cols), "\n")
  cat("  min:     ", min(vals, na.rm = TRUE), "\n")
  cat("  max:     ", max(vals, na.rm = TRUE), "\n")
  cat("  below ", lower, ": ", n_low, "\n", sep = "")
  cat("  above ", upper, ": ", n_high, "\n", sep = "")
  
  if (n_low == 0 && n_high == 0) {
    message("  PASS")
  } else {
    warning("  FAIL: values outside expected range for ", label)
  }
  
  invisible(NULL)
}


boundary_check <- function(mat, parameter, lower, upper, tolerance) {
  cols <- grep(paste0("^", parameter, "\\["), colnames(mat), value = TRUE)
  
  if (length(cols) == 0) {
    warning("No columns found for ", parameter)
    return(invisible(NULL))
  }
  
  vals <- mat[, cols, drop = FALSE]
  
  prop_near_lower <- colMeans(vals <= lower + tolerance, na.rm = TRUE)
  prop_near_upper <- colMeans(vals >= upper - tolerance, na.rm = TRUE)
  
  out <- tibble(
    parameter = cols,
    prop_near_lower = prop_near_lower,
    prop_near_upper = prop_near_upper
  )
  
  print(out)
  
  if (any(out$prop_near_lower > 0.05 | out$prop_near_upper > 0.05)) {
    warning(
      parameter,
      " has samples close to prior boundaries. This may indicate weak identifiability or overly tight priors."
    )
  } else {
    message("PASS: ", parameter, " not obviously stuck at prior boundaries.")
  }
  
  invisible(out)
}


make_spatial_initial_values <- function(bc, ids, nind, n.prim, first, K) {
  
  bc_for_init <- bc %>%
    mutate(.tattoo_order = match(tattoo, ids)) %>%
    filter(!is.na(.tattoo_order)) %>%
    arrange(.tattoo_order, primary, trap_season, date)
  
  S_init <- array(NA_real_, dim = c(nind, 2, n.prim))
  
  for (i in seq_len(nind)) {
    individual_data <- bc_for_init %>%
      filter(tattoo == ids[i]) %>%
      arrange(primary, trap_season, date)
    
    if (nrow(individual_data) == 0) {
      stop(
        "No capture rows found while creating initial S for individual index ",
        i, " (tattoo ", ids[i], ").",
        call. = FALSE
      )
    }
    
    first_x <- individual_data$col_index[1]
    first_y <- individual_data$row_index[1]
    
    if (is.na(first_x) || is.na(first_y)) {
      stop(
        "First observed spatial location is NA for individual ",
        ids[i], ". Check the settGrid join and SG names.",
        call. = FALSE
      )
    }
    
    # Fill all years with a sensible initial location.
    # This does NOT make the animal alive before first capture; it just avoids
    # NIMBLE receiving NA initial coordinates in unused array elements.
    last_x <- first_x
    last_y <- first_y
    
    for (k in seq_len(n.prim)) {
      matching_rows <- individual_data %>%
        filter(primary == k)
      
      if (nrow(matching_rows) > 0) {
        last_x <- matching_rows$col_index[1]
        last_y <- matching_rows$row_index[1]
      }
      
      S_init[i, 1, k] <- last_x
      S_init[i, 2, k] <- last_y
    }
  }
  
  if (anyNA(S_init)) {
    bad <- which(is.na(S_init), arr.ind = TRUE)
    
    bad_df <- tibble(
      individual_index = bad[, 1],
      coordinate = ifelse(bad[, 2] == 1, "x_col_index", "y_row_index"),
      primary_index = bad[, 3],
      tattoo = ids[bad[, 1]]
    ) %>%
      distinct() %>%
      arrange(individual_index, primary_index)
    
    print(head(bad_df, 30))
    
    stop(
      "Initial S array still contains NA values. The first 30 missing positions are printed above.",
      call. = FALSE
    )
  }
  
  # Latent alive/dead initial values.
  z_init <- matrix(0L, nrow = nind, ncol = n.prim)
  
  for (i in seq_len(nind)) {
    if (!is.na(first[i]) && !is.na(K[i]) && K[i] >= first[i]) {
      z_init[i, first[i]:K[i]] <- 1L
    }
  }
  
  d_init <- matrix(1, nrow = nind, ncol = max(1, n.prim - 1))
  theta_init <- matrix(0, nrow = nind, ncol = max(1, n.prim - 1))
  
  list(
    PL = c(0.5, 0.5),
    kappa = c(2, 2),
    sigma = c(8, 8),
    phi = matrix(runif((n.prim - 1) * 2, 0.7, 0.95), nrow = 2),
    dmean = c(8 + runif(1, -1, 1), 8 + runif(1, -1, 1)),
    z = z_init,
    d = d_init,
    theta = theta_init,
    S = S_init
  )
}


# ---- 3. Paths and data loading ---------------------------------------------

project_root <- find_project_root()

capture_data_path <- file.path(
  project_root,
  "02_data/01_processed_badger_data/badger_capture_diagnostic_cleaned_2024.rds"
)

spatial_object_path <- file.path(
  project_root,
  "02_data/04_saved_spatial_objects/spatial_grid_and_sett_objects.RData"
)

out_dir <- file.path(project_root, "04_model_outputs/01_saved_model_runs")
plot_dir <- file.path(project_root, "05_figures_and_animations/model_diagnostics")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

bc <- readRDS(capture_data_path) %>%
  droplevels()

load(spatial_object_path)


# ---- 4. Infection-history variables ----------------------------------------

stop_if_missing_cols(
  bc,
  c("tattoo", "date", "sett", "pm", "socg", "sex", "age", "trap_season"),
  "bc"
)

# Some projects may not contain every diagnostic column, so absent columns are
# created as zero rather than failing later.
bc <- replace_na_if_exists(bc, c("brock", "GAMMA", "statpak"), value = 0)

cult_cols <- grep("^Cult_", names(bc), value = TRUE)

if (length(cult_cols) == 0) {
  warning("No columns beginning with `Cult_` were found. Culture positivity will be treated as zero.")
  bc$culture_positive_any <- 0L
} else {
  bc <- bc %>%
    mutate(
      culture_positive_any = as.integer(
        rowSums(across(all_of(cult_cols), ~ replace_na(as.numeric(.x), 0) > 0)) > 0
      )
    )
}

bc <- bc %>%
  mutate(
    test_positive = as.integer(
      culture_positive_any == 1 |
        replace_na(as.numeric(brock), 0) > 0 |
        replace_na(as.numeric(statpak), 0) > 0 |
        replace_na(as.numeric(GAMMA), 0) > 0
    )
  )

individual_summary <- bc %>%
  group_by(tattoo) %>%
  summarise(
    positive_as_cub = as.integer(any(test_positive == 1 & age == 0, na.rm = TRUE)),
    ever_positive   = as.integer(any(test_positive == 1, na.rm = TRUE)),
    never_positive  = as.integer(all(test_positive == 0 | is.na(test_positive))),
    .groups = "drop"
  )

bc <- bc %>%
  select(-any_of(c("positive_as_cub", "ever_positive", "never_positive", "cub_positive"))) %>%
  left_join(individual_summary, by = "tattoo")

message("Overall infection-history summary before filtering:")
print(
  individual_summary %>%
    summarise(
      n_individuals = n(),
      n_positive_as_cub = sum(positive_as_cub == 1),
      n_ever_positive = sum(ever_positive == 1),
      n_never_positive = sum(never_positive == 1)
    )
)

# Clean two-group comparison:
#   1 = never positive
#   2 = positive as cub
#
# Individuals positive later in life but not as cubs are removed from this
# specific model because they do not belong cleanly to either group.
bc <- bc %>%
  filter(positive_as_cub == 1 | never_positive == 1) %>%
  mutate(
    infection_group = case_when(
      never_positive == 1 ~ 1L,
      positive_as_cub == 1 ~ 2L,
      TRUE ~ NA_integer_
    ),
    infection_group_label = factor(
      infection_group,
      levels = c(1, 2),
      labels = c("Never positive", "Test-positive as cub")
    )
  )

if (anyNA(bc$infection_group)) {
  stop("Some rows have missing infection_group after filtering. Check infection grouping logic.")
}

message("Model comparison groups after infection filtering:")
print(
  bc %>%
    distinct(tattoo, infection_group_label) %>%
    count(infection_group_label, name = "n_individuals")
)


# ---- 4b. Optional test subset ----------------------------------------------

if (!is.na(SAMPLE_N)) {
  set.seed(RANDOM_SEED)
  
  eligible_ids <- bc %>%
    distinct(tattoo, infection_group)
  
  if (SAMPLE_N >= nrow(eligible_ids)) {
    sampled_ids <- eligible_ids$tattoo
  } else {
    sample_per_group <- max(1, floor(SAMPLE_N / 2))
    
    ids_by_group <- split(eligible_ids$tattoo, eligible_ids$infection_group)
    
    sampled_ids <- unlist(
      lapply(ids_by_group, function(ids) {
        sample(ids, size = min(sample_per_group, length(ids)), replace = FALSE)
      }),
      use.names = FALSE
    )
    
    if (length(sampled_ids) < SAMPLE_N) {
      remaining_ids <- setdiff(eligible_ids$tattoo, sampled_ids)
      
      if (length(remaining_ids) > 0) {
        sampled_ids <- c(
          sampled_ids,
          sample(
            remaining_ids,
            size = min(SAMPLE_N - length(sampled_ids), length(remaining_ids)),
            replace = FALSE
          )
        )
      }
    }
  }
  
  bc <- bc %>%
    filter(tattoo %in% sampled_ids)
  
  message("Using test subset of ", length(unique(bc$tattoo)), " individuals.")
  
  print(
    bc %>%
      distinct(tattoo, infection_group_label) %>%
      count(infection_group_label, name = "n_individuals")
  )
}

# ---- 5. Spatial data preparation -------------------------------------------

load(spatial_object_path)

# ---- 5. Spatial data preparation using NEW habitat/SG grid ------------------
#
# Important distinction:
#
#   capture_loc_id / capture_loc_num = observation-model capture location
#                                     = cleaned socg, because settGrid only has
#                                       the coarser 23-ish spatial locations.
#
#   sett_id = cleaned exact sett/location name from bc$sett.
#             Kept as a diagnostic/descriptive column, but NOT used for X/H.
#
#   SG_id / SG_mat = landscape social-group territory ID from the new GIS layer.
#
#   habitat_mat = 1 land, 0 lake.
#
# The old settGrid geometry is reused, but its old row_index/col_index values
# are discarded and rebuilt from the new habitat/SG grid.


# ---- 5a. Preserve original capture-location spatial object ------------------

if (!exists("settGrid")) {
  stop("Object `settGrid` was not found after loading spatial_object_path.")
}

if (!inherits(settGrid, "sf")) {
  stop(
    "`settGrid` is not an sf object. ",
    "You need the original capture-location geometry to rebuild indices.",
    call. = FALSE
  )
}

settGrid_old_sf <- settGrid


# ---- 5b. Clean IDs in the capture data --------------------------------------

bc <- bc %>%
  mutate(
    # This is the spatial capture-location ID used by the model.
    # It should match the old settGrid$Sett values.
    capture_loc_id = iconv(socg, from = "latin1", to = "UTF-8", sub = ""),
    capture_loc_id = gsub(" ", "", capture_loc_id),
    
    # Exact sett/location name. Keep for diagnostics only.
    sett_id = iconv(sett, from = "latin1", to = "UTF-8", sub = ""),
    sett_id = gsub(" ", "", sett_id)
  )


# ---- 5c. Read the NEW final habitat/SG grid ---------------------------------

new_grid_path <- file.path(
  project_root,
  "02_data/03_spatial_inputs/Hetro_habitat_SGs/Habitat_Final.shp"
)

if (!file.exists(new_grid_path)) {
  stop(
    "Cannot find new grid shapefile at:\n",
    new_grid_path,
    "\nCheck the path and file name.",
    call. = FALSE
  )
}

new_grid <- sf::st_read(new_grid_path, quiet = FALSE)


# ---- 5d. Standardise/check new grid columns ---------------------------------

# The final grid must contain:
#   SG_id   = integer territory/buffer ID
#   habitat = 1 land, 0 lake

if (!"SG_id" %in% names(new_grid) && "SG_ID" %in% names(new_grid)) {
  new_grid <- new_grid %>% rename(SG_id = SG_ID)
}

if (!"SG_id" %in% names(new_grid) && "sg_id" %in% names(new_grid)) {
  new_grid <- new_grid %>% rename(SG_id = sg_id)
}

if (!"habitat" %in% names(new_grid) && "Habitat" %in% names(new_grid)) {
  new_grid <- new_grid %>% rename(habitat = Habitat)
}

required_grid_cols <- c("SG_id", "habitat")
missing_grid_cols <- setdiff(required_grid_cols, names(new_grid))

if (length(missing_grid_cols) > 0) {
  stop(
    "The new grid is missing required columns: ",
    paste(missing_grid_cols, collapse = ", "),
    call. = FALSE
  )
}

new_grid <- new_grid %>%
  mutate(
    SG_id = as.integer(SG_id),
    habitat = as.integer(habitat)
  )

if (anyNA(new_grid$SG_id)) {
  stop("new_grid$SG_id contains NA values. Fix this before modelling.")
}

if (anyNA(new_grid$habitat)) {
  stop("new_grid$habitat contains NA values. Fix this before modelling.")
}

if (!all(sort(unique(new_grid$habitat)) %in% c(0L, 1L))) {
  stop("new_grid$habitat must contain only 0 = lake and 1 = land.")
}


# ---- 5e. Check CRS consistency ----------------------------------------------

if (is.na(sf::st_crs(settGrid_old_sf))) {
  stop("settGrid has no CRS. Assign/check CRS before continuing.")
}

if (is.na(sf::st_crs(new_grid))) {
  stop("new_grid has no CRS. Assign/check CRS before continuing.")
}

if (sf::st_crs(settGrid_old_sf) != sf::st_crs(new_grid)) {
  message("Transforming settGrid to match new_grid CRS.")
  settGrid_old_sf <- sf::st_transform(settGrid_old_sf, sf::st_crs(new_grid))
}


# ---- 5f. Recreate row/column indices for the NEW grid -----------------------

# IMPORTANT: this must match the cell size used in QGIS.
cell_size <- 25L

new_grid_cent <- sf::st_centroid(new_grid)
grid_coords <- sf::st_coordinates(new_grid_cent)

bb <- sf::st_bbox(new_grid)

xmin <- as.numeric(bb["xmin"])
xmax <- as.numeric(bb["xmax"])
ymin <- as.numeric(bb["ymin"])
ymax <- as.numeric(bb["ymax"])

new_grid <- new_grid %>%
  mutate(
    grid_x = grid_coords[, 1],
    grid_y = grid_coords[, 2],
    col_index = as.integer(floor((grid_x - xmin) / cell_size)) + 1L,
    row_index = as.integer(floor((ymax - grid_y) / cell_size)) + 1L
  )

n_rows <- max(new_grid$row_index, na.rm = TRUE)
n_cols <- max(new_grid$col_index, na.rm = TRUE)

cat("\nNew grid dimensions:\n")
cat("  n_rows:           ", n_rows, "\n")
cat("  n_cols:           ", n_cols, "\n")
cat("  n cells expected: ", n_rows * n_cols, "\n")
cat("  n rows in grid:   ", nrow(new_grid), "\n")

if (n_rows * n_cols != nrow(new_grid)) {
  warning(
    "n_rows * n_cols does not equal nrow(new_grid). ",
    "This suggests the grid may not be a complete rectangle. ",
    "For this matrix lookup approach, a complete rectangular grid is strongly preferred."
  )
}


# ---- 5g. Build habitat and SG matrices from the NEW grid --------------------

habitat_mat <- matrix(NA_integer_, nrow = n_rows, ncol = n_cols)
SG_mat      <- matrix(NA_integer_, nrow = n_rows, ncol = n_cols)

habitat_mat[cbind(new_grid$row_index, new_grid$col_index)] <-
  as.integer(new_grid$habitat)

SG_mat[cbind(new_grid$row_index, new_grid$col_index)] <-
  as.integer(new_grid$SG_id)

if (anyNA(habitat_mat)) {
  stop(
    "habitat_mat contains NA values. ",
    "This usually means the exported grid is not a complete rectangle, ",
    "or row/column indexing is inconsistent.",
    call. = FALSE
  )
}

if (anyNA(SG_mat)) {
  stop(
    "SG_mat contains NA values. ",
    "Every grid cell must have an SG_id, including the buffer.",
    call. = FALSE
  )
}

cat("\nHabitat matrix check:\n")
print(table(as.vector(habitat_mat), useNA = "ifany"))

cat("\nSG matrix check:\n")
print(table(as.vector(SG_mat), useNA = "ifany"))


# ---- 5h. Standardise old settGrid capture-location IDs ----------------------

settGrid_old_tab <- settGrid_old_sf %>%
  sf::st_drop_geometry()

if (!"Sett" %in% names(settGrid_old_tab) && !"capture_loc_id" %in% names(settGrid_old_tab)) {
  stop(
    "settGrid does not contain a `Sett` or `capture_loc_id` column. ",
    "Need a capture-location ID to match cleaned bc$socg.",
    call. = FALSE
  )
}

if ("Sett" %in% names(settGrid_old_sf) && !"capture_loc_id" %in% names(settGrid_old_sf)) {
  settGrid_old_sf <- settGrid_old_sf %>%
    rename(capture_loc_id = Sett)
}

settGrid_old_sf <- settGrid_old_sf %>%
  mutate(
    capture_loc_id = iconv(capture_loc_id, from = "latin1", to = "UTF-8", sub = ""),
    capture_loc_id = gsub(" ", "", capture_loc_id)
  )


# ---- 5i. Check capture-location matching ------------------------------------

bc_capture_ids <- sort(unique(bc$capture_loc_id))

grid_capture_ids <- settGrid_old_sf %>%
  sf::st_drop_geometry() %>%
  distinct(capture_loc_id) %>%
  pull(capture_loc_id) %>%
  sort()

cat("\nCapture-location ID match check:\n")
cat("  unique bc$capture_loc_id:       ", length(bc_capture_ids), "\n")
cat("  unique settGrid$capture_loc_id: ", length(grid_capture_ids), "\n")
cat("  bc capture IDs matching settGrid: ",
    sum(bc_capture_ids %in% grid_capture_ids), "\n")

unmatched_bc_capture_locs <- setdiff(bc_capture_ids, grid_capture_ids)

if (length(unmatched_bc_capture_locs) > 0) {
  warning(
    length(unmatched_bc_capture_locs),
    " capture_loc_id values are not present in settGrid. First examples:\n",
    paste(head(unmatched_bc_capture_locs, 20), collapse = ", ")
  )
}

# Optional diagnostic showing why exact sett names are not used for spatial X/H.
bc_sett_ids <- sort(unique(bc$sett_id))

cat("\nExact sett-name diagnostic only:\n")
cat("  unique bc$sett_id values: ", length(bc_sett_ids), "\n")
cat("  These are not used as model detection locations unless you have a full sett-point layer.\n")


# ---- 5j. Assign every capture location to the NEW grid ----------------------
#
# Use nearest grid cell rather than st_within because points can sometimes sit
# exactly on cell boundaries.

nearest_cell <- sf::st_nearest_feature(settGrid_old_sf, new_grid)

settGrid_new_sf <- settGrid_old_sf %>%
  mutate(
    row_index = new_grid$row_index[nearest_cell],
    col_index = new_grid$col_index[nearest_cell],
    SG_id_landscape = new_grid$SG_id[nearest_cell],
    habitat = new_grid$habitat[nearest_cell]
  )

sett_distance_to_cell <- as.numeric(
  sf::st_distance(
    sf::st_geometry(settGrid_old_sf),
    sf::st_geometry(new_grid[nearest_cell, ]),
    by_element = TRUE
  )
)

cat("\nDistance from capture-location points to assigned grid cells:\n")
print(summary(sett_distance_to_cell))

if (max(sett_distance_to_cell, na.rm = TRUE) > cell_size) {
  warning(
    "Some capture-location points are more than one cell size away from their assigned grid cell. ",
    "This suggests some points may fall outside the new state-space."
  )
}

settGrid <- settGrid_new_sf %>%
  sf::st_drop_geometry() %>%
  mutate(
    row_index = as.integer(row_index),
    col_index = as.integer(col_index),
    SG_id_landscape = as.integer(SG_id_landscape),
    habitat = as.integer(habitat)
  )

dup_check <- settGrid %>%
  count(capture_loc_id) %>%
  filter(n > 1)

if (nrow(dup_check) > 0) {
  coord_dup_check <- settGrid %>%
    group_by(capture_loc_id) %>%
    summarise(
      n_rows = n(),
      n_unique_cells = n_distinct(paste(row_index, col_index)),
      .groups = "drop"
    ) %>%
    filter(n_unique_cells > 1)
  
  if (nrow(coord_dup_check) > 0) {
    print(coord_dup_check)
    stop(
      "Some capture_loc_id values occur multiple times with different grid cells. ",
      "Resolve duplicate capture-location IDs before modelling.",
      call. = FALSE
    )
  }
  
  settGrid <- settGrid %>%
    distinct(capture_loc_id, .keep_all = TRUE)
}

stop_if_missing_cols(
  settGrid,
  c("capture_loc_id", "row_index", "col_index"),
  "settGrid"
)


# ---- 5k. Check rebuilt capture-location indices -----------------------------

if (anyNA(settGrid$row_index)) {
  stop("Some settGrid$row_index values are NA after assigning to the new grid.")
}

if (anyNA(settGrid$col_index)) {
  stop("Some settGrid$col_index values are NA after assigning to the new grid.")
}

if (any(settGrid$row_index < 1L | settGrid$row_index > n_rows)) {
  stop("Some settGrid$row_index values fall outside the new matrix.")
}

if (any(settGrid$col_index < 1L | settGrid$col_index > n_cols)) {
  stop("Some settGrid$col_index values fall outside the new matrix.")
}

capture_loc_habitat_check <- habitat_mat[cbind(settGrid$row_index, settGrid$col_index)]
capture_loc_SG_check      <- SG_mat[cbind(settGrid$row_index, settGrid$col_index)]

cat("\nCapture-location habitat check, from matrix:\n")
print(table(capture_loc_habitat_check, useNA = "ifany"))

cat("\nCapture-location SG territory check, from matrix:\n")
print(table(capture_loc_SG_check, useNA = "ifany"))

if (any(capture_loc_habitat_check == 0L, na.rm = TRUE)) {
  warning(
    "Some capture locations fall in cells classified as lake. ",
    "Check lake polygons, grid resolution, and nearest-cell assignment."
  )
}


# ---- 5l. Save consistent spatial inputs -------------------------------------

consistent_spatial_inputs_path <- file.path(
  project_root,
  "02_data/04_saved_spatial_objects/consistent_new_grid_spatial_inputs.rds"
)

saveRDS(
  list(
    new_grid = new_grid,
    settGrid = settGrid,
    habitat_mat = habitat_mat,
    SG_mat = SG_mat,
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    cell_size = cell_size,
    n_rows = n_rows,
    n_cols = n_cols
  ),
  consistent_spatial_inputs_path
)

message("Saved consistent spatial inputs to: ", consistent_spatial_inputs_path)


# ---- 5m. Join updated capture-location indices onto capture data ------------

capture_locs_all <- settGrid$capture_loc_id

bc <- bc %>%
  filter(capture_loc_id %in% capture_locs_all) %>%
  left_join(
    settGrid %>%
      dplyr::select(capture_loc_id, row_index, col_index, SG_id_landscape, habitat),
    by = "capture_loc_id"
  ) %>%
  mutate(
    date = lubridate::ymd(date),
    primary_year = lubridate::year(date)
  ) %>%
  filter(primary_year > 1981) %>%
  arrange(date) %>%
  mutate(
    primary = factor(primary_year, levels = as.character(1982:2020))
  ) %>%
  dplyr::select(
    tattoo, date,
    sett, sett_id,
    socg, capture_loc_id,
    pm, sex, age,
    positive_as_cub, ever_positive, never_positive,
    infection_group, infection_group_label,
    primary, trap_season,
    row_index, col_index,
    SG_id_landscape, habitat
  ) %>%
  ungroup() %>%
  filter(!is.na(primary))

bc <- bc %>%
  group_by(tattoo) %>%
  filter(!all(pm == "Yes")) %>%
  ungroup() %>%
  arrange(tattoo, date)

if (nrow(bc) == 0) {
  stop("No observations remain after spatial/date/pm filtering.")
}

bc <- bc %>%
  mutate(
    primary = as.numeric(primary) + 1L,
    trap_season = as.integer(trap_season)
  )

if (anyNA(bc$trap_season)) {
  stop("trap_season contains NA after conversion to integer.")
}

if (anyNA(bc$row_index) || anyNA(bc$col_index)) {
  stop("bc has NA row_index/col_index after joining updated settGrid.")
}


# ---- 6. Capture-history arrays and death constraints ------------------------

n.prim <- max(bc$primary)
n.sec <- max(bc$trap_season)
dt <- rep(1, times = n.prim - 1)
nind <- length(unique(bc$tattoo))

bc <- bc %>%
  group_by(tattoo) %>%
  mutate(
    minPrimary = min(primary),
    maxPrimary = max(primary)
  ) %>%
  group_by(tattoo, primary) %>%
  mutate(lastSecondary = max(trap_season)) %>%
  ungroup() %>%
  group_by(tattoo) %>%
  mutate(
    firstSecondary = ifelse(primary == minPrimary, min(trap_season), NA_integer_)
  ) %>%
  fill(firstSecondary, .direction = "downup") %>%
  mutate(
    lastSecondary = ifelse(primary == maxPrimary, lastSecondary, NA_integer_)
  ) %>%
  fill(lastSecondary, .direction = "downup") %>%
  mutate(
    death.occasion = ifelse(pm == "Yes", primary, n.prim + 1L),
    death.occasion = min(death.occasion, na.rm = TRUE),
    death.secondary = ifelse(pm == "Yes", trap_season, 0L)
  ) %>%
  mutate(
    lastSecondary = ifelse(primary == maxPrimary & pm == "Yes", lastSecondary - 1L, lastSecondary),
    maxPrimary = ifelse(pm == "Yes" & lastSecondary == 0L, maxPrimary - 1L, maxPrimary),
    lastSecondary = ifelse(lastSecondary == 0L & pm == "Yes", 4L, lastSecondary)
  ) %>%
  arrange(tattoo, date) %>%
  mutate(maxPrimary = min(maxPrimary)) %>%
  ungroup()

bc <- bc %>%
  group_by(tattoo) %>%
  mutate(maxPrimary = ifelse(maxPrimary == death.occasion, maxPrimary - 1L, maxPrimary)) %>%
  ungroup()

bc <- bc %>%
  group_by(tattoo) %>%
  filter(min(maxPrimary, na.rm = TRUE) >= min(minPrimary, na.rm = TRUE)) %>%
  ungroup()

nind <- length(unique(bc$tattoo))

if (nind == 0) {
  stop("No individuals remain after death-constraint filtering.")
}

# Capture locations used by the observation model.
capture_locs <- levels(as.factor(bc$capture_loc_id))

settGrid <- settGrid %>%
  filter(capture_loc_id %in% capture_locs)

ids <- unique(bc$tattoo)

first <- bc %>%
  distinct(tattoo, .keep_all = TRUE) %>%
  pull(minPrimary)

J <- matrix(n.sec, nind, n.prim)

death.occasion <- bc %>%
  distinct(tattoo, .keep_all = TRUE) %>%
  pull(death.occasion)

death.primary <- death.occasion
death.primary[death.primary > n.prim] <- NA

death.secondary <- bc %>%
  group_by(tattoo) %>%
  mutate(death.secondary = max(death.secondary)) %>%
  distinct(tattoo, .keep_all = TRUE) %>%
  pull(death.secondary)

death.secondary[death.secondary == 0] <- NA

K <- rep(n.prim, nind)
K[death.occasion < n.prim] <- death.occasion[death.occasion < n.prim]

for (i in seq_along(death.primary)) {
  if (!is.na(death.primary[i])) {
    K[i] <- death.primary[i] - 1L
    
    if (!is.na(death.secondary[i])) {
      J[i, death.primary[i]] <- death.secondary[i]
    }
  }
}

capture_loc_map <- setNames(seq_along(capture_locs), capture_locs)

bc$capture_loc_num <- as.integer(capture_loc_map[bc$capture_loc_id])

if (anyNA(bc$capture_loc_num)) {
  stop("Some capture_loc_id values could not be mapped to capture-location indices.")
}

H <- array(1L, dim = c(nind, n.sec, n.prim), dimnames = list(ids, NULL, NULL))

for (row_i in seq_len(nrow(bc))) {
  p <- bc$primary[row_i]
  s <- bc$trap_season[row_i]
  ind <- bc$tattoo[row_i]
  
  H[ind, s, p] <- bc$capture_loc_num[row_i] + 1L
}

R <- length(unique(bc$capture_loc_num))

settGrid$capture_loc_id <- factor(settGrid$capture_loc_id, levels = capture_locs)
settGrid <- settGrid[order(settGrid$capture_loc_id), ]
rownames(settGrid) <- NULL

X <- settGrid %>%
  dplyr::select(col_index, row_index)

last.prim <- bc %>%
  distinct(tattoo, .keep_all = TRUE) %>%
  pull(maxPrimary)

z_data <- matrix(NA, nrow = nind, ncol = n.prim)

death.occNA <- death.occasion
death.occNA[death.occNA > n.prim] <- NA

for (i in seq_len(nind)) {
  z_data[i, first[i]:last.prim[i]] <- 1L
  
  if (!is.na(death.occNA[i])) {
    z_data[i, death.occNA[i]:n.prim] <- 0L
  }
}

group <- bc %>%
  distinct(tattoo, .keep_all = TRUE) %>%
  pull(infection_group)

first_location <- bc %>%
  arrange(date) %>%
  distinct(tattoo, .keep_all = TRUE) %>%
  pull(capture_loc_num)

ord <- order(K - first)

K <- K[ord]
J <- J[ord, ]
first <- first[ord]
H <- H[ord, , ]
death.occasion <- death.occasion[ord]
z_data <- z_data[ord, ]
group <- group[ord]
first_location <- first_location[ord]
ids <- ids[ord]
last.prim <- last.prim[ord]

if (!all(group %in% c(1L, 2L))) {
  stop("`group` must only contain 1 and 2.")
}

if (any(K < first)) {
  stop("Some individuals have K < first after sorting; check death constraints.")
}

N <- c(sum((K - first) == 0), length(first))

message("N[1] individuals with no later movement/survival interval: ", N[1])
message("N[2] total individuals: ", N[2])

if (N[1] < 1) {
  stop(
    "N[1] is zero. The current NIMBLE model has a for-loop `for (i in 1:N[1])`, ",
    "which requires at least one individual in this group. You can fix this later ",
    "by rewriting the model loops, but for now this dataset/subset needs N[1] > 0.",
    call. = FALSE
  )
}

if (N[1] >= N[2]) {
  stop(
    "All individuals have K == first, so there are no individuals with a movement/survival interval. ",
    "Increase SAMPLE_N or use the full dataset.",
    call. = FALSE
  )
}


# ---- 6b. Input validity checks ---------------------------------------------

cat("\nInput validity checks before nimbleModel()\n")
cat("Any NA in X: ", anyNA(as.matrix(X)), "\n")
cat("Any NA in H: ", anyNA(H), "\n")
cat("Any NA in first_location: ", anyNA(first_location), "\n")
cat("Any NA in group: ", anyNA(group), "\n")
cat("Any NA in first: ", anyNA(first), "\n")
cat("Any NA in K: ", anyNA(K), "\n")
cat("Any K < first: ", any(K < first), "\n")
cat("Any NA in habitat_mat: ", anyNA(habitat_mat), "\n")
cat("Any NA in SG_mat: ", anyNA(SG_mat), "\n")

if (anyNA(as.matrix(X))) stop("X contains NA values.")
if (anyNA(H)) stop("H contains NA values.")
if (anyNA(first_location)) stop("first_location contains NA values.")
if (anyNA(group)) stop("group contains NA values.")
if (any(K < first)) stop("Some individuals have K < first.")
if (anyNA(habitat_mat)) stop("habitat_mat contains NA values.")
if (anyNA(SG_mat)) stop("SG_mat contains NA values.")

if (!all(sort(unique(as.vector(habitat_mat))) %in% c(0L, 1L))) {
  stop("habitat_mat must only contain 0 = lake and 1 = land.")
}

if (!identical(dim(habitat_mat), dim(SG_mat))) {
  stop("habitat_mat and SG_mat have different dimensions.")
}

n_rows <- nrow(habitat_mat)
n_cols <- ncol(habitat_mat)


# ---- 6b. Input validity checks ---------------------------------------------

cat("\nInput validity checks before nimbleModel()\n")
cat("Any NA in X: ", anyNA(as.matrix(X)), "\n")
cat("Any NA in H: ", anyNA(H), "\n")
cat("Any NA in first_location: ", anyNA(first_location), "\n")
cat("Any NA in group: ", anyNA(group), "\n")
cat("Any NA in first: ", anyNA(first), "\n")
cat("Any NA in K: ", anyNA(K), "\n")
cat("Any K < first: ", any(K < first), "\n")
cat("Any NA in habitat_mat: ", anyNA(habitat_mat), "\n")
cat("Any NA in SG_mat: ", anyNA(SG_mat), "\n")

if (anyNA(as.matrix(X))) stop("X contains NA values.")
if (anyNA(H)) stop("H contains NA values.")
if (anyNA(first_location)) stop("first_location contains NA values.")
if (anyNA(group)) stop("group contains NA values.")
if (any(K < first)) stop("Some individuals have K < first.")
if (anyNA(habitat_mat)) stop("habitat_mat contains NA values.")
if (anyNA(SG_mat)) stop("SG_mat contains NA values.")

if (!all(sort(unique(as.vector(habitat_mat))) %in% c(0L, 1L))) {
  stop("habitat_mat must only contain 0 = lake and 1 = land.")
}

if (!identical(dim(habitat_mat), dim(SG_mat))) {
  stop("habitat_mat and SG_mat have different dimensions.")
}

n_rows <- nrow(habitat_mat)
n_cols <- ncol(habitat_mat)


# ---- 7. NIMBLE model --------------------------------------------------------

# ---- 7. NIMBLE model --------------------------------------------------------

code <- nimbleCode({
  
  ## PRIORS AND CONSTRAINTS
  
  for (grp in 1:2) {
    
    kappa[grp] ~ dunif(0.25, 10)
    sigma[grp] ~ dunif(0.25, 30)
    
    PL[grp] ~ dunif(0.01, 0.99)
    log.lambda0[grp] <- log(-log(1 - PL[grp]))
    lambda[grp] <- exp(log.lambda0[grp])
    
    for (k in 1:(n.prim - 1)) {
      phi[grp, k] ~ dunif(0.001, 0.999)
    }
    
    dmean[grp] ~ dunif(0.25, 100)
    dlambda[grp] <- 1 / dmean[grp]
  }
  
  # Shared landscape penalties for first implementation.
  lake_penalty ~ dexp(1)
  boundary_penalty ~ dexp(1)
  
  
  ## MODEL
  
  # Individuals seen only in their final/censored primary session.
  for (i in 1:N[1]) {
    
    z[i, first[i]] ~ dbern(1)
    
    S[i, 1, first[i]] <- X[first_location[i], 1]
    S[i, 2, first[i]] <- X[first_location[i], 2]
    
    col_S[i, first[i]] <- max(1, min(n_cols, trunc(S[i, 1, first[i]])))
    row_S[i, first[i]] <- max(1, min(n_rows, trunc(S[i, 2, first[i]])))
    
    habitat_here[i, first[i]] <- habitat_mat[row_S[i, first[i]], col_S[i, first[i]]]
    SG_here[i, first[i]]      <- SG_mat[row_S[i, first[i]], col_S[i, first[i]]]
    
    p_land[i, first[i]] <- habitat_here[i, first[i]] +
      (1 - habitat_here[i, first[i]]) * exp(-lake_penalty)
    
    land_ok[i, first[i]] ~ dbern(p_land[i, first[i]])
    
    g[i, first[i], 1] <- 0
    
    for (r in 1:R) {
      D[i, r, first[i]] <- sqrt(
        pow(S[i, 1, first[i]] - X[r, 1], 2) +
          pow(S[i, 2, first[i]] - X[r, 2], 2)
      )
      
      g[i, first[i], r + 1] <-
        exp(-pow(D[i, r, first[i]] / sigma[group[i]], kappa[group[i]]))
    }
    
    G[i, first[i]] <- sum(g[i, first[i], 1:(R + 1)])
    
    for (j in 1:J[i, first[i]]) {
      
      P[i, j, first[i]] <- 1 - exp(-lambda[group[i]] * G[i, first[i]])
      
      captureProb[i, first[i], j] <-
        step(H[i, j, first[i]] - 2) *
        (g[i, first[i], H[i, j, first[i]]] / (G[i, first[i]] + 0.000000001)) *
        P[i, j, first[i]] +
        (1 - step(H[i, j, first[i]] - 2)) *
        (1 - P[i, j, first[i]])
      
      Ones[i, j, first[i]] ~ dbern(captureProb[i, first[i], j])
    }
  }
  
  
  # Individuals with at least one later primary session.
  for (i in (N[1] + 1):N[2]) {
    
    z[i, first[i]] ~ dbern(1)
    
    S[i, 1, first[i]] <- X[first_location[i], 1]
    S[i, 2, first[i]] <- X[first_location[i], 2]
    
    col_S[i, first[i]] <- max(1, min(n_cols, trunc(S[i, 1, first[i]])))
    row_S[i, first[i]] <- max(1, min(n_rows, trunc(S[i, 2, first[i]])))
    
    habitat_here[i, first[i]] <- habitat_mat[row_S[i, first[i]], col_S[i, first[i]]]
    SG_here[i, first[i]]      <- SG_mat[row_S[i, first[i]], col_S[i, first[i]]]
    
    p_land[i, first[i]] <- habitat_here[i, first[i]] +
      (1 - habitat_here[i, first[i]]) * exp(-lake_penalty)
    
    land_ok[i, first[i]] ~ dbern(p_land[i, first[i]])
    
    g[i, first[i], 1] <- 0
    
    for (r in 1:R) {
      D[i, r, first[i]] <- sqrt(
        pow(S[i, 1, first[i]] - X[r, 1], 2) +
          pow(S[i, 2, first[i]] - X[r, 2], 2)
      )
      
      g[i, first[i], r + 1] <-
        exp(-pow(D[i, r, first[i]] / sigma[group[i]], kappa[group[i]]))
    }
    
    G[i, first[i]] <- sum(g[i, first[i], 1:(R + 1)])
    
    for (j in 1:J[i, first[i]]) {
      
      P[i, j, first[i]] <- 1 - exp(-lambda[group[i]] * G[i, first[i]])
      
      captureProb[i, first[i], j] <-
        step(H[i, j, first[i]] - 2) *
        (g[i, first[i], H[i, j, first[i]]] / (G[i, first[i]] + 0.000000001)) *
        P[i, j, first[i]] +
        (1 - step(H[i, j, first[i]] - 2)) *
        (1 - P[i, j, first[i]])
      
      Ones[i, j, first[i]] ~ dbern(captureProb[i, first[i], j])
    }
    
    
    # Later primary sessions.
    for (k in (first[i] + 1):K[i]) {
      
      Palive[i, k - 1] <- z[i, k - 1] * phi[group[i], k - 1]
      
      z[i, k] ~ dbern(Palive[i, k - 1] * step(death.occasion[i] - k))
      
      theta[i, k - 1] ~ dunif(-3.141593, 3.141593)
      
      d[i, k - 1] ~ dexp(dlambda[group[i]])
      
      S[i, 1, k] <- S[i, 1, k - 1] + d[i, k - 1] * cos(theta[i, k - 1])
      S[i, 2, k] <- S[i, 2, k - 1] + d[i, k - 1] * sin(theta[i, k - 1])
      
      col_S[i, k] <- max(1, min(n_cols, trunc(S[i, 1, k])))
      row_S[i, k] <- max(1, min(n_rows, trunc(S[i, 2, k])))
      
      habitat_here[i, k] <- habitat_mat[row_S[i, k], col_S[i, k]]
      SG_here[i, k]      <- SG_mat[row_S[i, k], col_S[i, k]]
      
      p_land[i, k] <- habitat_here[i, k] +
        (1 - habitat_here[i, k]) * exp(-lake_penalty)
      
      land_ok[i, k] ~ dbern(p_land[i, k])
      
      same_SG[i, k] <- equals(SG_here[i, k], SG_here[i, k - 1])
      
      p_boundary[i, k] <- same_SG[i, k] +
        (1 - same_SG[i, k]) * exp(-boundary_penalty)
      
      boundary_ok[i, k] ~ dbern(p_boundary[i, k])
      
      g[i, k, 1] <- 0
      
      for (r in 1:R) {
        
        D[i, r, k] <- sqrt(
          pow(S[i, 1, k] - X[r, 1], 2) +
            pow(S[i, 2, k] - X[r, 2], 2)
        )
        
        g[i, k, r + 1] <-
          exp(-pow(D[i, r, k] / sigma[group[i]], kappa[group[i]]))
      }
      
      G[i, k] <- sum(g[i, k, 1:(R + 1)])
      
      for (j in 1:J[i, k]) {
        
        P[i, j, k] <- (1 - exp(-lambda[group[i]] * G[i, k])) * z[i, k]
        
        captureProb[i, k, j] <-
          step(H[i, j, k] - 2) *
          (g[i, k, H[i, j, k]] / (G[i, k] + 0.000000001)) *
          P[i, j, k] +
          (1 - step(H[i, j, k] - 2)) *
          (1 - P[i, j, k])
        
        Ones[i, j, k] ~ dbern(captureProb[i, k, j])
      }
    }
  }
})


# ---- 8. Initial values ------------------------------------------------------
set.seed(RANDOM_SEED)

inits <- make_spatial_initial_values(
  bc = bc,
  ids = ids,
  nind = nind,
  n.prim = n.prim,
  first = first,
  K = K
)

if (anyNA(inits$S)) {
  stop("Initial S array contains NA values after robust filling. Check spatial joins.")
}

message("Initial S array successfully created with no NA values.")

# ---- 9. Build model ---------------------------------------------------------
land_ok <- matrix(1L, nrow = nind, ncol = n.prim)

boundary_ok <- matrix(NA_integer_, nrow = nind, ncol = n.prim)
if (n.prim >= 2) {
  boundary_ok[, 2:n.prim] <- 1L
}

consts <- list(
  R = nrow(X),
  N = N,
  K = K,
  J = J,
  first = first,
  X = as.matrix(X),
  n.prim = n.prim,
  H = H,
  death.occasion = death.occasion,
  group = group,
  first_location = first_location,
  
  # These are constants because they define dimensions/limits
  n_rows = n_rows,
  n_cols = n_cols
)

data <- list(
  Ones = array(1L, dim(H)),
  z = z_data,
  
  # Pseudo-observation data
  land_ok = land_ok,
  boundary_ok = boundary_ok,
  
  # IMPORTANT:
  # These are dynamically indexed inside the model,
  # so they must be supplied as data, not constants.
  habitat_mat = habitat_mat,
  SG_mat = SG_mat
)

model <- nimbleModel(
  code,
  constants = consts,
  data = data,
  inits = inits,
  dimensions = list(
    habitat_mat = dim(habitat_mat),
    SG_mat = dim(SG_mat),
    land_ok = dim(land_ok),
    boundary_ok = dim(boundary_ok),
    Ones = dim(H),
    z = dim(z_data)
  )
)

# ---- 10. Configure MCMC -----------------------------------------------------

if (MONITOR_LATENT) {
  
  valid_S_nodes <- character()
  valid_z_nodes <- character()
  
  # Because primary was shifted by +1, primary index 1 is usually unused.
  # Therefore we only monitor years from first[i]:K[i], and drop any year 1
  # just in case it sneaks in.
  for (i in seq_len(nind)) {
    
    if (is.na(first[i]) || is.na(K[i])) next
    if (K[i] < first[i]) next
    
    valid_years <- seq(from = first[i], to = K[i])
    valid_years <- valid_years[valid_years > 1L]
    
    if (length(valid_years) == 0) next
    
    valid_S_nodes <- c(
      valid_S_nodes,
      paste0("S[", i, ", 1, ", valid_years, "]"),
      paste0("S[", i, ", 2, ", valid_years, "]")
    )
    
    valid_z_nodes <- c(
      valid_z_nodes,
      paste0("z[", i, ", ", valid_years, "]")
    )
  }
  
  latent_monitor_nodes <- c(valid_S_nodes, valid_z_nodes)
  
  message("Monitoring ", length(valid_S_nodes), " valid S nodes.")
  message("Monitoring ", length(valid_z_nodes), " valid z nodes.")
  message("Total latent monitor nodes: ", length(latent_monitor_nodes))
  
  time1_nodes <- grep("\\[, [12], 1\\]$|\\[, 1\\]$", latent_monitor_nodes, value = TRUE)
  
  if (length(time1_nodes) > 0) {
    warning(
      "The latent monitor list contains time-1 nodes. First examples:\n",
      paste(head(time1_nodes, 20), collapse = ", ")
    )
  } else {
    message("No time-1 latent S/z nodes included in latent_monitor_nodes.")
  }
  
  config <- configureMCMC(
    model,
    monitors = c("phi", "kappa", "sigma", "PL", "dmean"),
    thin = 1,
    monitors2 = latent_monitor_nodes,
    thin2 = 20
  )
  
} else {
  
  message("MONITOR_LATENT is FALSE: only main parameters will be monitored.")
  
  config <- configureMCMC(
    model,
    monitors = c("phi", "kappa", "sigma", "PL", "dmean"),
    thin = 1
  )
}

rMCMC <- buildMCMC(config)

# Compile the model first
cModel <- compileNimble(model)

# Then compile the MCMC using the model project
cMCMC <- compileNimble(rMCMC, project = model)


# ---- 11. Run MCMC -----------------------------------------------------------

TEST_RUN <- FALSE

if (TEST_RUN) {
  niter <- MCMC_NITER_TEST
  nburnin <- MCMC_NBURN_TEST
  nchains <- MCMC_NCHAINS_TEST
} else {
  niter <- MCMC_NITER_FULL
  nburnin <- MCMC_NBURN_FULL
  nchains <- MCMC_NCHAINS_FULL
}

output_suffix <- if (TEST_RUN) "TEST" else "FULL"

message(
  "Running MCMC with niter = ", niter,
  ", nburnin = ", nburnin,
  ", nchains = ", nchains
)

system.time(
  run <- runMCMC(
    cMCMC,
    niter = niter,
    nburnin = nburnin,
    nchains = nchains,
    progressBar = TRUE,
    summary = FALSE,
    samplesAsCodaMCMC = TRUE
  )
)

output_path <- file.path(
  out_dir,
  paste0("spatial_cmr_infection_group_habitat_and_SG", output_suffix, ".rds")
)

saveRDS(run, output_path)

message("Saved model output to: ", output_path)


# ---- 12. Post-run validation ------------------------------------------------

cat("\n================ MAIN PARAMETER SAMPLE CHECKS ================\n")

samples_mat <- as.matrix(run$samples)

cat("Number of posterior draws: ", nrow(samples_mat), "\n")
cat("Number of monitored parameters: ", ncol(samples_mat), "\n")
cat("NA values:  ", sum(is.na(samples_mat)), "\n")
cat("NaN values: ", sum(is.nan(samples_mat)), "\n")
cat("Inf values: ", sum(is.infinite(samples_mat)), "\n")

bad_main_cols <- colnames(samples_mat)[
  apply(samples_mat, 2, function(x) {
    any(is.na(x) | is.nan(x) | is.infinite(x))
  })
]

if (length(bad_main_cols) > 0) {
  warning(
    "Some main monitored parameters contain NA/NaN/Inf values:\n",
    paste(head(bad_main_cols, 30), collapse = ", ")
  )
} else {
  message("PASS: No NA/NaN/Inf values in main monitored parameters.")
}


cat("\n================ PARAMETER RANGE CHECKS ================\n")

check_range(samples_mat, "^PL\\[",     lower = 0, upper = 1,   label = "PL")
check_range(samples_mat, "^phi\\[",    lower = 0, upper = 1,   label = "phi")
check_range(samples_mat, "^sigma\\[",  lower = 0, upper = Inf, label = "sigma")
check_range(samples_mat, "^kappa\\[",  lower = 0, upper = Inf, label = "kappa")
check_range(samples_mat, "^dmean\\[",  lower = 0, upper = Inf, label = "dmean")


cat("\n================ PRIOR BOUNDARY CHECKS ================\n")

boundary_check(samples_mat, "PL",     lower = 0.01, upper = 0.99, tolerance = 0.01)
boundary_check(samples_mat, "sigma",  lower = 0.25, upper = 30,   tolerance = 0.5)
boundary_check(samples_mat, "kappa",  lower = 0.25, upper = 10,   tolerance = 0.25)
boundary_check(samples_mat, "dmean",  lower = 0.25, upper = 100,  tolerance = 1)


cat("\n================ LATENT SAMPLE CHECKS ================\n")

if (!is.null(run$samples2)) {
  
  samples2_mat <- as.matrix(run$samples2)
  
  cat("Number of latent posterior draws: ", nrow(samples2_mat), "\n")
  cat("Number of latent monitored nodes: ", ncol(samples2_mat), "\n")
  cat("NA values:  ", sum(is.na(samples2_mat)), "\n")
  cat("NaN values: ", sum(is.nan(samples2_mat)), "\n")
  cat("Inf values: ", sum(is.infinite(samples2_mat)), "\n")
  
  bad_latent_cols <- colnames(samples2_mat)[
    apply(samples2_mat, 2, function(x) {
      any(is.na(x) | is.nan(x) | is.infinite(x))
    })
  ]
  
  if (length(bad_latent_cols) > 0) {
    warning(
      "Some latent monitored nodes contain NA/NaN/Inf. First affected nodes:\n",
      paste(head(bad_latent_cols, 30), collapse = ", ")
    )
    
    cat("\nBreakdown of bad latent nodes:\n")
    cat("Bad S nodes: ", sum(grepl("^S\\[", bad_latent_cols)), "\n")
    cat("Bad z nodes: ", sum(grepl("^z\\[", bad_latent_cols)), "\n")
    
  } else {
    message("PASS: No NA/NaN/Inf values in monitored latent S/z nodes.")
  }
  
} else {
  message("No samples2 object found. Latent S/z nodes were not monitored.")
}


cat("\n================ MCMC DIAGNOSTICS ================\n")

print(summary(run$samples))

if (inherits(run$samples, "mcmc.list") && length(run$samples) > 1) {
  
  gelman <- tryCatch(
    coda::gelman.diag(run$samples, multivariate = FALSE),
    error = function(e) {
      warning("Could not calculate Gelman diagnostics: ", conditionMessage(e))
      NULL
    }
  )
  
  if (!is.null(gelman)) {
    print(gelman)
    
    rhat_values <- gelman$psrf[, "Point est."]
    
    cat("\nRhat summary:\n")
    print(summary(rhat_values))
    
    if (any(rhat_values > 1.1, na.rm = TRUE)) {
      warning("Some Rhat values are > 1.1. Expected in short test runs, not acceptable for final inference.")
    } else {
      message("PASS: No Rhat values > 1.1.")
    }
  }
  
  ess <- tryCatch(
    coda::effectiveSize(run$samples),
    error = function(e) {
      warning("Could not calculate effective sample sizes: ", conditionMessage(e))
      NULL
    }
  )
  
  if (!is.null(ess)) {
    cat("\nEffective sample size summary:\n")
    print(summary(ess))
    
    low_ess <- ess[ess < 100]
    
    if (length(low_ess) > 0) {
      warning(
        length(low_ess),
        " parameters have ESS < 100. Expected in short test runs, not acceptable for final inference."
      )
      print(head(sort(low_ess), 20))
    } else {
      message("PASS: All ESS values >= 100.")
    }
  }
  
} else {
  message("Only one chain detected or samples are not an mcmc.list; skipping Rhat diagnostics.")
}


# ---- 13. Movement contrast --------------------------------------------------

cat("\n================ MOVEMENT CONTRAST ================\n")

samples_df <- as.data.frame(samples_mat)

if (all(c("dmean[1]", "dmean[2]") %in% names(samples_df))) {
  
  dM_difference <- samples_df %>%
    transmute(
      dmean_never_positive = `dmean[1]`,
      dmean_positive_as_cub = `dmean[2]`,
      difference = dmean_positive_as_cub - dmean_never_positive
    )
  
  summary_difference <- dM_difference %>%
    summarise(
      median_difference = median(difference, na.rm = TRUE),
      mean_difference = mean(difference, na.rm = TRUE),
      lower_95 = quantile(difference, 0.025, na.rm = TRUE),
      upper_95 = quantile(difference, 0.975, na.rm = TRUE),
      prob_positive_as_cub_greater = mean(difference > 0, na.rm = TRUE)
    )
  
  print(summary_difference)
  
  write.csv(
    summary_difference,
    file = file.path(plot_dir, paste0("dmean_difference_summary_", output_suffix, ".csv")),
    row.names = FALSE
  )
  
} else {
  warning("Could not find both dmean[1] and dmean[2] in posterior samples.")
}


# ---- 14. Basic plots --------------------------------------------------------

try(
  mcmcplots::mcmcplot(
    run$samples,
    parms = c("kappa", "sigma", "PL", "dmean")
  ),
  silent = TRUE
)

if (all(c("dmean[1]", "dmean[2]") %in% names(samples_df))) {
  
  dM_samples <- samples_df %>%
    dplyr::select(starts_with("dmean")) %>%
    pivot_longer(everything(), names_to = "Infection", values_to = "Estimate") %>%
    mutate(
      Infection = factor(
        Infection,
        levels = c("dmean[1]", "dmean[2]"),
        labels = c("Never positive", "Test-positive as cub")
      )
    )
  
  p_dmean <- ggplot(dM_samples, aes(x = Estimate, colour = Infection, fill = Infection)) +
    geom_density(alpha = 0.3) +
    labs(
      x = "Mean annual activity-centre movement distance",
      y = "Density",
      colour = "Badger status",
      fill = "Badger status"
    ) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      text = element_text(size = 17)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  
  print(p_dmean)
  
  ggsave(
    filename = file.path(plot_dir, paste0("dmean_density_infection_group_", output_suffix, ".png")),
    plot = p_dmean,
    width = 9,
    height = 6,
    dpi = 300
  )
}


message("\nMaster candidate model script complete.")
