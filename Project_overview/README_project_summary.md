# NewBadgers organised project

This is a cleaned and reorganised copy of the original `NewBadgers` zip.

The untouched original zip is preserved here:

```text
00_original_archive/NewBadgers_original_untouched.zip
```

The active working copy has been reorganised into clearer folders, and many scripts have been renamed to describe what they appear to do. I have not rewritten or refactored the R code itself, so some hard-coded paths, `setwd()` calls, `readRDS()`, `load()`, shapefile paths, and output paths may need updating before scripts run from this new layout.

A full map of old to new locations is here:

```text
01_project_overview/file_rename_and_location_map.csv
```

## Main interpretation

The project appears to have developed from simple multistate badger CMR models into a NIMBLE implementation of an Ergon & Gardner-style spatial robust design model. The most advanced candidate script is the infection-history grouping model with adjusted first centre-of-activity location:

```text
03_scripts/05_ergon_gardner_spatial_cmr/infection_group/02_badger_spatial_cmr_infection_group_adjusted_first_location.R
```

## Suggested workflow order

1. `03_scripts/01_data_preparation/01_merge_clean_capture_diagnostic_data.R`
2. `03_scripts/02_spatial_setup/01_create_spatial_grid_objects.R`
3. `03_scripts/02_spatial_setup/02_plot_spatial_grid_and_setts.R`
4. Optional prototype models in `03_scripts/03_early_multistate_models/`, `03_scripts/04_temporary_emigration_models/`, and `03_scripts/06_spatial_cjs_prototypes/`
5. Main spatial CMR scripts in `03_scripts/05_ergon_gardner_spatial_cmr/`
6. Plotting and animation scripts in `03_scripts/07_plotting_and_animation/`

## Important caution before continuing

Before treating the infection-group model as final, check the infection classification and grouping logic. In the old script, the grouping index is still named `sex` in places even when being used for infection status. Also check whether the intended comparison is truly `test-positive as cub` versus `never positive`, rather than `test-positive as cub` versus all other badgers.

## Folder guide

```text
00_original_archive/          Untouched original zip
01_project_overview/          Summary, run order, and file map
02_data/                      Cleaned data, spatial inputs, saved spatial objects
03_scripts/                   Organised R scripts by project stage
04_model_outputs/             Saved model outputs found in original zip
05_figures_and_animations/    Existing GIF outputs
06_rstudio_and_git_metadata/  RStudio/Git metadata copied where useful
```
