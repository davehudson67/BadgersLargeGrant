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

TEST_RUN <- TRUE          # TRUE = short MCMC for debugging
SAMPLE_N <- 100           # use NA_integer_ for full dataset
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
library(mcmcplots)


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

bc <- bc %>%
  mutate(
    SG = iconv(socg, from = "latin1", to = "UTF-8", sub = ""),
    SG = gsub(" ", "", SG)
  )

settGrid <- settGrid %>%
  st_drop_geometry() %>%
  rename(SG = Sett)

stop_if_missing_cols(settGrid, c("SG", "row_index", "col_index"), "settGrid")

if (exists("studyArea")) rm(studyArea)

SGs_all <- settGrid$SG

bc <- bc %>%
  filter(SG %in% SGs_all) %>%
  left_join(settGrid, by = "SG") %>%
  mutate(
    date = ymd(date),
    primary_year = year(date)
  ) %>%
  filter(primary_year > 1981) %>%
  arrange(date) %>%
  mutate(
    primary = factor(primary_year, levels = as.character(1982:2020))
  ) %>%
  dplyr::select(
    tattoo, date, sett, pm, socg, SG, sex, age,
    positive_as_cub, ever_positive, never_positive,
    infection_group, infection_group_label,
    primary, trap_season, row_index, col_index
  ) %>%
  ungroup() %>%
  filter(!is.na(primary))

# Remove individuals represented only by post-mortem rows.
bc <- bc %>%
  group_by(tattoo) %>%
  filter(!all(pm == "Yes")) %>%
  ungroup() %>%
  arrange(tattoo, date)

if (nrow(bc) == 0) {
  stop("No observations remain after spatial/date/pm filtering.")
}

# Keep the original +1 primary indexing for now.
bc <- bc %>%
  mutate(
    primary = as.numeric(primary) + 1L,
    trap_season = as.integer(trap_season)
  )

if (anyNA(bc$trap_season)) {
  stop("trap_season contains NA after conversion to integer.")
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

SGs <- levels(as.factor(bc$SG))

settGrid <- settGrid %>%
  filter(SG %in% SGs)

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

sett_map <- setNames(seq_along(SGs), SGs)
bc$SGnum <- as.integer(sett_map[bc$SG])

if (anyNA(bc$SGnum)) {
  stop("Some SG values could not be mapped to sett indices.")
}

H <- array(1L, dim = c(nind, n.sec, n.prim), dimnames = list(ids, NULL, NULL))

for (row_i in seq_len(nrow(bc))) {
  p <- bc$primary[row_i]
  s <- bc$trap_season[row_i]
  ind <- bc$tattoo[row_i]
  
  H[ind, s, p] <- bc$SGnum[row_i] + 1L
}

R <- length(unique(bc$SGnum))

settGrid$SG <- factor(settGrid$SG, levels = SGs)
settGrid <- settGrid[order(settGrid$SG), ]
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
  pull(SGnum)

# Sort individuals so those with K == first come first, matching the original
# Ergon & Gardner-style code split.
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

if (anyNA(as.matrix(X))) stop("X contains NA values.")
if (anyNA(first_location)) stop("first_location contains NA values.")
if (anyNA(group)) stop("group contains NA values.")
if (any(K < first)) stop("Some individuals have K < first.")


# ---- 7. NIMBLE model --------------------------------------------------------

code <- nimbleCode({
  
  ## PRIORS AND CONSTRAINTS
  
  for (grp in 1:2) {
    
    # Detection kernel shape and scale.
    # Lower bounds are away from zero to avoid 0/0 or 0^0 behaviour in G.
    kappa[grp] ~ dunif(0.25, 10)
    sigma[grp] ~ dunif(0.25, 30)
    
    # Baseline capture probability.
    PL[grp] ~ dunif(0.01, 0.99)
    log.lambda0[grp] <- log(-log(1 - PL[grp]))
    lambda[grp] <- exp(log.lambda0[grp])
    
    # Annual survival.
    for (k in 1:(n.prim - 1)) {
      phi[grp, k] ~ dunif(0.001, 0.999)
    }
    
    # Annual activity-centre movement distance.
    dmean[grp] ~ dunif(0.25, 100)
    dlambda[grp] <- 1 / dmean[grp]
  }
  
  
  ## MODEL
  
  # Individuals seen only in their final/censored primary session.
  for (i in 1:N[1]) {
    
    z[i, first[i]] ~ dbern(1)
    
    S[i, 1, first[i]] <- X[first_location[i], 1]
    S[i, 2, first[i]] <- X[first_location[i], 2]
    
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
  first_location = first_location
)

data <- list(
  Ones = array(1L, dim(H)),
  z = z_data
)

model <- nimbleModel(
  code,
  constants = consts,
  data = data,
  inits = inits
)

# Useful if something fails at model build:
# model$initializeInfo()

cModel <- compileNimble(model)


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
cMCMC <- compileNimble(rMCMC)


# ---- 11. Run MCMC -----------------------------------------------------------

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
  paste0("spatial_cmr_infection_group_adjusted_first_location_", output_suffix, ".rds")
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
