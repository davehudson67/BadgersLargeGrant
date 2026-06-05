# ============================================================================
# Animate posterior mean badger activity centres with Pr(alive)
# ============================================================================
library(tidyverse)
library(stringr)
library(ggplot2)
library(gganimate)
library(gifski)
library(sf)

# ---- 1. Find project root ---------------------------------------------------

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
        "Open the NewBadgers_organised folder as your RStudio project.",
        call. = FALSE
      )
    }
    
    current <- parent
  }
}

project_root <- find_project_root()

# ---- 2. User options --------------------------------------------------------

MODEL_SUFFIX <- "FULL"   # use "TEST" if animating a test run

ALIVE_ALPHA_MIN <- 0.05
ALIVE_ALPHA_MAX <- 1

# Set to 0 if you want all monitored states shown with fading.
# Set to 0.5 if you want individuals to disappear once more likely dead than alive.
MIN_PR_ALIVE_TO_SHOW <- 0.05

FPS <- 4
WIDTH <- 850
HEIGHT <- 850

# ---- 3. Paths ---------------------------------------------------------------

model_path <- file.path(
  project_root,
  "04_model_outputs/01_saved_model_runs",
  paste0("spatial_cmr_infection_group_habitat_and_SGFULL", ".rds")
)

spatial_object_path <- file.path(
  project_root,
  "02_data/04_saved_spatial_objects/spatial_grid_and_sett_objects.RData"
)

gif_dir <- file.path(
  project_root,
  "05_figures_and_animations/movement_gifs"
)

dir.create(gif_dir, recursive = TRUE, showWarnings = FALSE)

gif_path <- file.path(
  gif_dir,
  paste0("badger_activity_centres_pr_alive_", MODEL_SUFFIX, ".gif")
)

# ---- 4. Load model output and spatial objects -------------------------------

run <- readRDS(model_path)

if (is.null(run$samples2)) {
  stop(
    "run$samples2 is NULL. Re-run the model with MONITOR_LATENT <- TRUE.",
    call. = FALSE
  )
}

samples2_mat <- as.matrix(run$samples2)

if (!any(grepl("^S\\[", colnames(samples2_mat)))) {
  stop("No S nodes found in run$samples2.", call. = FALSE)
}

if (!any(grepl("^z\\[", colnames(samples2_mat)))) {
  stop("No z nodes found in run$samples2.", call. = FALSE)
}

load(spatial_object_path)

# ---- 5. Recreate individual lookup table -----------------------------------
#
# Your saved run object contains posterior samples but not necessarily the
# individual metadata. For animation labels/colours we recreate the same
# filtering and ordering used in the model script.
#
# This should match the model script. If you later save `ids`, `group`, `first`
# and `K` directly as a metadata RDS, use that instead.

capture_data_path <- file.path(
  project_root,
  "02_data/01_processed_badger_data/badger_capture_diagnostic_cleaned_2024.rds"
)

bc <- readRDS(capture_data_path) %>%
  droplevels()

replace_na_if_exists <- function(dat, cols, value = 0) {
  for (nm in cols) {
    if (!nm %in% names(dat)) {
      dat[[nm]] <- value
    } else {
      dat[[nm]][is.na(dat[[nm]])] <- value
    }
  }
  dat
}

bc <- replace_na_if_exists(bc, c("brock", "GAMMA", "statpak"), value = 0)

cult_cols <- grep("^Cult_", names(bc), value = TRUE)

if (length(cult_cols) == 0) {
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
    ever_positive = as.integer(any(test_positive == 1, na.rm = TRUE)),
    never_positive = as.integer(all(test_positive == 0 | is.na(test_positive))),
    .groups = "drop"
  )

bc <- bc %>%
  select(-any_of(c("positive_as_cub", "ever_positive", "never_positive", "cub_positive"))) %>%
  left_join(individual_summary, by = "tattoo") %>%
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

bc <- bc %>%
  mutate(
    SG = iconv(socg, from = "latin1", to = "UTF-8", sub = ""),
    SG = gsub(" ", "", SG)
  )

settGrid_clean <- settGrid %>%
  st_drop_geometry() %>%
  rename(SG = Sett)

SGs_all <- settGrid_clean$SG

bc <- bc %>%
  filter(SG %in% SGs_all) %>%
  left_join(settGrid_clean, by = "SG") %>%
  mutate(
    date = lubridate::ymd(date),
    primary_year = lubridate::year(date)
  ) %>%
  filter(primary_year > 1981) %>%
  arrange(date) %>%
  mutate(
    primary = factor(primary_year, levels = as.character(1982:2020))
  ) %>%
  filter(!is.na(primary)) %>%
  group_by(tattoo) %>%
  filter(!all(pm == "Yes")) %>%
  ungroup() %>%
  arrange(tattoo, date) %>%
  mutate(
    primary = as.numeric(primary) + 1L,
    trap_season = as.integer(trap_season)
  )

n.prim <- max(bc$primary)
n.sec <- max(bc$trap_season)
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
  tidyr::fill(firstSecondary, .direction = "downup") %>%
  mutate(
    lastSecondary = ifelse(primary == maxPrimary, lastSecondary, NA_integer_)
  ) %>%
  tidyr::fill(lastSecondary, .direction = "downup") %>%
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
  ungroup() %>%
  group_by(tattoo) %>%
  mutate(maxPrimary = ifelse(maxPrimary == death.occasion, maxPrimary - 1L, maxPrimary)) %>%
  ungroup() %>%
  group_by(tattoo) %>%
  filter(min(maxPrimary, na.rm = TRUE) >= min(minPrimary, na.rm = TRUE)) %>%
  ungroup()

ids <- unique(bc$tattoo)

first <- bc %>%
  distinct(tattoo, .keep_all = TRUE) %>%
  pull(minPrimary)

death.occasion <- bc %>%
  distinct(tattoo, .keep_all = TRUE) %>%
  pull(death.occasion)

K <- rep(n.prim, length(ids))
K[death.occasion < n.prim] <- death.occasion[death.occasion < n.prim]

death.primary <- death.occasion
death.primary[death.primary > n.prim] <- NA

for (i in seq_along(death.primary)) {
  if (!is.na(death.primary[i])) {
    K[i] <- death.primary[i] - 1L
  }
}

group <- bc %>%
  distinct(tattoo, .keep_all = TRUE) %>%
  pull(infection_group)

last.prim <- bc %>%
  distinct(tattoo, .keep_all = TRUE) %>%
  pull(maxPrimary)

ord <- order(K - first)

ids <- ids[ord]
first <- first[ord]
K <- K[ord]
group <- group[ord]
last.prim <- last.prim[ord]

individual_lookup <- tibble(
  individual = seq_along(ids),
  tattoo = ids,
  first = first,
  K = K,
  last_prim = last.prim,
  infection_group = group,
  infection_group_label = factor(
    group,
    levels = c(1, 2),
    labels = c("Never positive", "Test-positive as cub")
  )
)

# ---- 6. Extract posterior S -------------------------------------------------

S_cols <- grep("^S\\[", colnames(samples2_mat), value = TRUE)

S_long <- samples2_mat[, S_cols, drop = FALSE] %>%
  as_tibble(.name_repair = "minimal") %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(
    cols = - .draw,
    names_to = "node",
    values_to = "value"
  ) %>%
  filter(!is.na(value), !is.nan(value), is.finite(value)) %>%
  mutate(
    parsed = str_match(node, "^S\\[(\\d+),\\s*(\\d+),\\s*(\\d+)\\]$"),
    individual = as.integer(parsed[, 2]),
    coordinate = as.integer(parsed[, 3]),
    time = as.integer(parsed[, 4])
  ) %>%
  select(.draw, individual, coordinate, time, value)

S_summary <- S_long %>%
  group_by(individual, coordinate, time) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = coordinate,
    values_from = mean_value,
    names_prefix = "S_"
  ) %>%
  rename(
    S_x = S_1,
    S_y = S_2
  )

# ---- 7. Extract posterior z / Pr(alive) -------------------------------------

z_cols <- grep("^z\\[", colnames(samples2_mat), value = TRUE)

z_long <- samples2_mat[, z_cols, drop = FALSE] %>%
  as_tibble(.name_repair = "minimal") %>%
  mutate(.draw = row_number()) %>%
  pivot_longer(
    cols = - .draw,
    names_to = "node",
    values_to = "z"
  ) %>%
  filter(!is.na(z), !is.nan(z), is.finite(z)) %>%
  mutate(
    parsed = str_match(node, "^z\\[(\\d+),\\s*(\\d+)\\]$"),
    individual = as.integer(parsed[, 2]),
    time = as.integer(parsed[, 3])
  ) %>%
  select(.draw, individual, time, z)

z_summary <- z_long %>%
  group_by(individual, time) %>%
  summarise(
    pr_alive = mean(z, na.rm = TRUE),
    .groups = "drop"
  )

# ---- 8. Build animation dataframe ------------------------------------------

animation_data <- S_summary %>%
  left_join(z_summary, by = c("individual", "time")) %>%
  left_join(individual_lookup, by = "individual") %>%
  mutate(
    pr_alive = replace_na(pr_alive, 0),
    year = 1981 + time - 1L
  ) %>%
  filter(
    time >= first,
    pr_alive >= MIN_PR_ALIVE_TO_SHOW
  )

if (nrow(animation_data) == 0) {
  stop("animation_data has zero rows after filtering. Try MIN_PR_ALIVE_TO_SHOW <- 0.", call. = FALSE)
}

# ---- Smooth/interpolate annual activity-centre positions --------------------
# The model estimates annual activity centres. For visualisation only, we
# interpolate between annual posterior mean locations so points move smoothly.
#
# This DOES NOT add biological information. It is just a display choice.

N_FRAMES_BETWEEN_YEARS <- 24

# Use a lower threshold before smoothing so animals can fade in/out rather than
# abruptly appearing/disappearing.
animation_data_for_smoothing <- animation_data %>%
  arrange(individual, time) %>%
  filter(!is.na(S_x), !is.na(S_y), !is.na(year), !is.na(pr_alive))

animation_data_smooth <- animation_data_for_smoothing %>%
  group_by(individual) %>%
  group_modify(~ {
    dat <- .x %>% arrange(year)
    
    # If an individual has only one annual point, keep it as a single static row.
    if (nrow(dat) == 1) {
      dat$frame_year <- dat$year
      return(dat)
    }
    
    segments <- vector("list", nrow(dat) - 1)
    
    for (ii in seq_len(nrow(dat) - 1)) {
      start_row <- dat[ii, ]
      end_row   <- dat[ii + 1, ]
      
      # Number of interpolated steps scales with the gap between years.
      year_gap <- max(1, end_row$year - start_row$year)
      n_steps <- N_FRAMES_BETWEEN_YEARS * year_gap
      
      frame_years <- seq(
        from = start_row$year,
        to = end_row$year,
        length.out = n_steps + 1
      )
      
      # Avoid duplicate frames at joins between segments.
      if (ii < nrow(dat) - 1) {
        frame_years <- frame_years[-length(frame_years)]
      }
      
      prop <- seq(0, 1, length.out = length(frame_years))
      
      segment <- start_row[rep(1, length(frame_years)), ]
      
      segment$frame_year <- frame_years
      segment$S_x <- start_row$S_x + prop * (end_row$S_x - start_row$S_x)
      segment$S_y <- start_row$S_y + prop * (end_row$S_y - start_row$S_y)
      segment$pr_alive <- start_row$pr_alive + prop * (end_row$pr_alive - start_row$pr_alive)
      
      segments[[ii]] <- segment
    }
    
    bind_rows(segments)
  }) %>%
  ungroup() %>%
  filter(!is.na(frame_year), !is.na(S_x), !is.na(S_y)) %>%
  mutate(
    frame_year = as.numeric(frame_year),
    infection_group_label = factor(
      infection_group_label,
      levels = c("Never positive", "Test-positive as cub")
    )
  )

# Quick checks before animating
message("Rows in annual animation_data: ", nrow(animation_data))
message("Rows in smoothed animation_data_smooth: ", nrow(animation_data_smooth))
message("Number of individuals in smoothed data: ", n_distinct(animation_data_smooth$individual))
message("Frame-year range: ", paste(range(animation_data_smooth$frame_year, na.rm = TRUE), collapse = " to "))

print(
  animation_data_smooth %>%
    count(infection_group_label, name = "n_rows")
)

if (nrow(animation_data_smooth) == 0) {
  stop("animation_data_smooth has zero rows. Check filtering and Pr(alive) threshold.")
}

# ---- 9. Build grid background ----------------------------------------------

grid_df <- settGrid_clean %>%
  transmute(
    SG,
    col = col_index,
    row = row_index
  )

# ---- 10. Animation 1: smooth moving points only -----------------------------
gif_path_points <- file.path(
  gif_dir,
  paste0("badger_activity_centres_smooth_points_pr_alive_", MODEL_SUFFIX, ".gif")
)

p_points <- ggplot() +
  geom_point(
    data = grid_df,
    aes(x = col, y = row),
    colour = "grey85",
    size = 0.8,
    alpha = 0.5
  ) +
  geom_point(
    data = animation_data_smooth,
    aes(
      x = S_x,
      y = S_y,
      colour = infection_group_label,
      alpha = pr_alive,
      group = individual
    ),
    size = 2
  ) +
  scale_alpha_continuous(
    limits = c(0, 1),
    range = c(0.02, 0.9),
    name = "Pr(alive)"
  ) +
  coord_equal() +
  theme_classic() +
  labs(
    title = "Posterior mean activity centres: {round(frame_time, 1)}",
    subtitle = "Positions interpolated between annual estimates; transparency shows posterior Pr(alive)",
    x = "Grid column",
    y = "Grid row",
    colour = "Infection group"
  ) +
  gganimate::transition_time(frame_year) +
  gganimate::ease_aes("linear")

animation_points <- animate(
  p_points,
  fps = 6,
  width = WIDTH,
  height = HEIGHT,
  renderer = gifski_renderer(gif_path_points)
)

animation_points

message("Saved smooth points animation to: ", gif_path_points)


# ---- 11. Animation 2: smooth points with short trails ------------------------
#
# Optional version. The trails make movement easier to see, but can get busy.
# shadow_wake() gives a short fading trail rather than drawing all historical
# paths for all individuals.

gif_path_trails <- file.path(
  gif_dir,
  paste0("badger_activity_centres_smooth_trails_pr_alive_", MODEL_SUFFIX, ".gif")
)

p_trails <- ggplot() +
  geom_point(
    data = grid_df,
    aes(x = col, y = row),
    colour = "grey85",
    size = 0.8,
    alpha = 0.5
  ) +
  geom_point(
    data = animation_data_smooth,
    aes(
      x = S_x,
      y = S_y,
      colour = infection_group_label,
      alpha = pr_alive,
      group = individual
    ),
    size = 2
  ) +
  scale_alpha_continuous(
    limits = c(0, 1),
    range = c(0.02, 0.9),
    name = "Pr(alive)"
  ) +
  coord_equal() +
  theme_classic() +
  labs(
    title = "Posterior mean activity centres: {round(frame_time, 1)}",
    subtitle = "Short trails show recent interpolated movement; transparency shows posterior Pr(alive)",
    x = "Grid column",
    y = "Grid row",
    colour = "Infection group"
  ) +
  gganimate::transition_time(frame_year) +
  gganimate::ease_aes("linear") +
  gganimate::shadow_wake(
    wake_length = 0.08,
    alpha = FALSE,
    size = FALSE
  )

animation_trails <- animate(
  p_trails,
  fps = 6,
  width = WIDTH,
  height = HEIGHT,
  renderer = gifski_renderer(gif_path_trails)
)

animation_trails

message("Saved smooth trails animation to: ", gif_path_trails)