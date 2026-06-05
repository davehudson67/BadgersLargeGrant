# BadgersLargeGrant script guide

This README describes the R scripts currently contained in the `BadgersLargeGrant` project folder. It is intended as a practical guide to what each script does, where each model idea came from, and which scripts are likely to be active versus historical/prototype material.

The project has developed from simple capture-mark-recapture examples, through multistate/social-group movement models, into an Ergon & Gardner-style spatial capture-recapture model for badger movement and survival. The most advanced current model appears to be the infection-history spatial CMR model in:

```text
Scripts/05_ergon_gardner_spatial_cmr/infection_group/
```

## Key terminology used in the scripts

| Term | Meaning in this project |
|---|---|
| `tattoo` | Individual badger ID. |
| `sett` | The exact sett/location recorded for a capture. There are many of these. In the current spatial model, these are mostly metadata unless a full point layer exists for every sett. |
| `socg` / `capture_loc_id` | The coarser social-group/capture-location identifier used by the existing spatial lookup object. This is the detection-location ID used to build the observation model. |
| `row_index`, `col_index` | Matrix/grid indices used as the coordinate system inside the current NIMBLE spatial model. |
| `X` | Matrix of fixed detection/capture locations, with columns `col_index` and `row_index`. |
| `S` | Latent centre of activity for each individual and primary session. In the current model this is in grid-index coordinates, not British National Grid metre coordinates. |
| `H` | Capture-history array: individual × secondary session × primary session. Values index the capture location, with 1 usually meaning not captured. |
| `z` | Alive/dead latent state. |
| `SG_id` | GIS-derived social-group territory ID in the landscape grid. This is separate from `socg`; it is intended for boundary-crossing penalties. |
| `habitat` | GIS-derived habitat code: 1 = land, 0 = lake. |

## Recommended high-level workflow

A clean active workflow would be:

```text
1. Data preparation
   Scripts/01_merge_clean_capture_diagnostic_data.R

2. Spatial grid setup
   Scripts/01_create_spatial_grid_objects.R
   Scripts/02_plot_spatial_grid_and_setts.R
   plus the newer GIS-derived habitat/SG grid workflow

3. Main spatial CMR modelling
   Scripts/05_ergon_gardner_spatial_cmr/infection_group/

4. Diagnostics and outputs
   model output checks within the main model script
   Scripts/07_plotting_and_animation/
```

The early multistate scripts and spatial CJS prototypes are useful for model development history, but they are probably not part of the active final analysis unless you intentionally revisit those model structures.

---

# Script-by-script guide

## Data and spatial setup scripts

### `Scripts/01_merge_clean_capture_diagnostic_data.R`

**Model origin:** Not a model script. Data preparation.

**What it does in plain terms:**

Combines capture records, badger metadata, post-mortem information and diagnostic/infection information into a cleaner badger capture dataset. It standardises dates and infection-related fields, handles missing sex values, and creates a processed data object for later modelling.

**Role in the project:** Active or near-active. This should be treated as the starting point for the modelling workflow.

**Outputs likely used later:**

```text
Data/badger_capture_diagnostic_cleaned_2024.rds
```

**Notes:** This is the place to check infection-test definitions before final modelling.

---

### `Scripts/01_create_spatial_grid_objects.R`

**Model origin:** Not a statistical model. Spatial preprocessing for grid-index-based spatial CMR.

**What it does in plain terms:**

Reads the study-area grid and sett-location shapefiles, converts map coordinates into matrix row/column indices, and creates the spatial objects that the early spatial CMR scripts use. It creates grid-style spatial inputs such as `settGrid`, `studyArea`, and related matrix/index objects.

**Role in the project:** Historically important, but now partly superseded by the newer habitat/SG grid. The concept is still important: all model coordinates must line up with the same grid.

**Caution:** If the GIS grid extent has changed, old `row_index` and `col_index` values from this script should not be reused without rebuilding them.

---

### `Scripts/02_plot_spatial_grid_and_setts.R`

**Model origin:** Not a statistical model. Spatial checking/visualisation.

**What it does in plain terms:**

Plots the study-area grid and sett/capture-location positions so you can check whether the grid, study area and spatial points line up properly.

**Role in the project:** Useful diagnostic script. It should be retained because spatial indexing mistakes are one of the easiest ways to break the model.

---

### `Scripts/03_social_group_name_changes_checks.R`

**Model origin:** Not a formal model. Exploratory data check.

**What it does in plain terms:**

Looks at how often individuals appear to change social group/capture location across their capture histories. It calculates individual-level movement/change summaries and summaries by social group.

**Role in the project:** Useful exploratory script for understanding whether social-group switching is common enough to model.

---

## Early multistate and temporary-emigration models

These scripts are development/prototype material. They helped explore how to represent badger movement between social groups or capture states before the project moved towards spatial robust-design CMR.

### `Scripts/01_paper_code_multistate_reference.R`

**Model origin:** Spatial Cormack-Jolly-Seber-style reference/simulation code.

**Likely reference family:** Cormack-Jolly-Seber models: Cormack (1964), Jolly (1965), Seber (1965), extended here to a spatial CJS setting with locations and dispersal.

**What it does in plain terms:**

Simulates animals moving through an irregularly shaped study area, with survival, recapture and spatial displacement. It then fits or prepares a spatial CJS-style model using known first capture information.

**Role in the project:** Reference/prototype only. Useful for understanding the origins of the spatial survival model ideas, but not the main badger analysis.

---

### `Scripts/01_temporary_emigration_multistate_model.R`

**Model origin:** Temporary-emigration model from the Bayesian Population Analysis framework; comments refer to “Temporary Emigration Model P282 in BPA book”.

**Likely reference family:** Robust design and temporary emigration models, especially Pollock’s robust design and the Kéry & Schaub Bayesian Population Analysis examples.

**What it does in plain terms:**

Treats individuals as able to be available or unavailable for capture between occasions. This is useful when animals may still be alive but temporarily outside the sampled area or otherwise not catchable.

**Role in the project:** Prototype/development script. It explores a related problem, but the current spatial CMR model deals with movement more explicitly through latent activity centres.

---

### `Scripts/02_simple_3state_multistate_model.R`

**Model origin:** Simple multistate capture-recapture model.

**Likely reference family:** Multistate CMR models, including Arnason/Schwarz/Brownie-style models for movement among states.

**What it does in plain terms:**

Simplifies each badger’s status into three states: alive at the same site, alive at a new site, or dead. It estimates survival, recapture probability, and movement to a new site.

**Role in the project:** Early prototype. Useful conceptually, but too simple for the final spatial questions.

---

### `Scripts/03_multisite_model_shorthand.R`

**Model origin:** Multisite/multistate capture-recapture model.

**Likely reference family:** Multistate CMR movement models, where animals can move among multiple discrete sites/states.

**What it does in plain terms:**

Extends the simple three-state idea to many sites. Instead of just “same site” versus “new site”, individuals can move among a set of named locations/social groups. The model estimates survival, recapture and movement probabilities among sites.

**Role in the project:** Prototype. Good for understanding social-group transition ideas, but it is a discrete-state model rather than the current continuous-space activity-centre model.

---

### `Scripts/04_multisite_model_longhand.R`

**Model origin:** Multisite/multistate capture-recapture model.

**Likely reference family:** Multistate CMR models.

**What it does in plain terms:**

A longer, more explicit version of the multisite model. It lays out transition probabilities among many sites and a dead state. It is easier to inspect conceptually but less scalable than a compact/shorthand implementation.

**Role in the project:** Prototype/reference. Useful if you need to revisit discrete social-group transition models.

---

### `Scripts/05_simulated_3state_cjs_model.R`

**Model origin:** Simulated three-state CJS/multistate example.

**Likely reference family:** Cormack-Jolly-Seber and multistate CMR models.

**What it does in plain terms:**

Simulates capture-recapture data for a small three-state system, then fits a NIMBLE model. It is mainly a teaching/debugging example for understanding state-transition and observation matrices.

**Role in the project:** Example/prototype only.

---

## Ergon & Gardner-style spatial CMR scripts

These are the main spatial model development scripts. They are based on a spatial robust-design capture-recapture framework in which each animal has a latent centre of activity that can move through time. Detection depends on distance between the latent centre and fixed capture locations. This lets the model separate apparent disappearance into mortality versus movement/dispersal more explicitly than standard CJS models.

**Model origin:** Ergon & Gardner-style spatial robust-design CMR. The reference scripts and filenames refer to Ergon & Gardner 2013/2014-style code. The relevant published model family is usually cited as Ergon & Gardner’s robust-design spatial capture-recapture approach for separating mortality and emigration/dispersal.

### `Scripts/05_ergon_gardner_spatial_cmr/baseline_and_reference/01_ergon_gardner_2013_working_nimble_reference.R`

**Model origin:** NIMBLE translation/adaptation of an Ergon & Gardner reference model.

**What it does in plain terms:**

Provides a working reference version of the spatial robust-design model in NIMBLE. It includes latent activity centres, survival, detection, movement and dispersal components.

**Role in the project:** Reference script. Useful for checking how the badger-specific models relate back to the original model structure.

---

### `Scripts/05_ergon_gardner_spatial_cmr/baseline_and_reference/02_ergon_gardner_2013_reference_model.R`

**Model origin:** JAGS reference model for Ergon & Gardner-style spatial robust-design CMR.

**What it does in plain terms:**

Holds the original/reference model specification. It is not tailored to the badger data, but it shows the model structure that later badger scripts adapt.

**Role in the project:** Keep as a reference. Do not edit unless documenting changes.

---

### `Scripts/05_ergon_gardner_spatial_cmr/badger_baseline/01_first_badger_adaptation_of_ergon_gardner_model.R`

**Model origin:** First badger-data adaptation of the Ergon & Gardner model.

**What it does in plain terms:**

Takes the reference spatial CMR idea and begins applying it to badger capture data. It builds capture-history arrays, fixed capture-location coordinates, alive/dead states, movement distances and detection probabilities.

**Role in the project:** Development history. Useful for tracing how the final model was built.

---

### `Scripts/05_ergon_gardner_spatial_cmr/badger_baseline/02_badger_spatial_cmr_baseline_no_covariates.R`

**Model origin:** Baseline badger spatial CMR model based on Ergon & Gardner.

**What it does in plain terms:**

Fits a spatial CMR model to the badger data without the infection-history grouping. It estimates survival, detection, space-use/detection parameters and movement for the whole analysed badger set.

**Role in the project:** Important baseline. This is useful for comparing against later models with infection group or sex effects.

---

### `Scripts/05_ergon_gardner_spatial_cmr/sex_effects/01_badger_spatial_cmr_sex_effects.R`

**Model origin:** Sex-effect extension of the baseline Ergon & Gardner-style spatial CMR model.

**What it does in plain terms:**

Extends the badger spatial model so parameters can vary by sex. This is useful for asking whether male and female badgers differ in movement, survival, detection or space use.

**Role in the project:** Model variant. It may be useful as a comparison model but is not the current infection-history focus.

---

### `Scripts/05_ergon_gardner_spatial_cmr/infection_group/01_badger_spatial_cmr_infection_group_initial.R`

**Model origin:** Infection-history extension of the Ergon & Gardner-style spatial CMR model.

**What it does in plain terms:**

Begins separating badgers into infection-history groups, especially never-positive versus test-positive as cub, and allows spatial movement/survival parameters to differ between those groups.

**Role in the project:** Early infection-group version.

---

### `Scripts/05_ergon_gardner_spatial_cmr/infection_group/02_badger_spatial_cmr_infection_group_adjusted_first_location.R`

**Model origin:** Infection-history spatial CMR model with adjusted first activity-centre location.

**What it does in plain terms:**

Fits the Ergon & Gardner-style spatial model while comparing movement and survival between infection-history groups. The first activity centre is anchored to the first observed capture-location/social-group location rather than being freely estimated from a broad prior. This helps stabilise the model.

**Role in the project:** Important candidate model.

---

### `Scripts/05_ergon_gardner_spatial_cmr/infection_group/badger_spatial_cmr_infection_group_adjusted_first_location_v2_no_S_NA.R`

**Model origin:** Current/master candidate infection-history spatial CMR model.

**What it does in plain terms:**

This is the most developed candidate script in the folder. It prepares infection-history groups, builds the spatial capture-history arrays, fixes the first latent activity centre to the first observed capture location, models movement between years, estimates survival and detection, and monitors latent activity centres for later animation.

**Role in the project:** Treat this as the main active model script unless replaced by a newer landscape-aware version. Recent development has focused on adding GIS-derived lake/habitat and social-group boundary matrices to this model.

**Important current caution:** The existing script works in row/column grid-index coordinates. If the state-space grid changes, the capture-location indices, `X`, `S`, `habitat_mat`, and `SG_mat` must all be rebuilt from the same grid.

---

### `Scripts/05_ergon_gardner_spatial_cmr/development_and_debugging/01_test_log_probability_eg_badger_model.R`

**Model origin:** Debugging support for Ergon & Gardner-style badger model.

**What it does in plain terms:**

Checks whether the model has valid log probabilities and helps identify where bad likelihoods, impossible states, `NA`, `NaN` or `-Inf` values arise.

**Role in the project:** Keep as debugging support.

---

### `Scripts/05_ergon_gardner_spatial_cmr/development_and_debugging/02_badger_spatial_cmr_covariate_development.R`

**Model origin:** Covariate-development version of the spatial CMR model.

**What it does in plain terms:**

Explores adding covariates to the baseline spatial CMR structure. This appears to be a development script rather than a final analysis script.

**Role in the project:** Archive/development.

---

### `Scripts/05_ergon_gardner_spatial_cmr/development_and_debugging/03_initial_values_for_baseline_badger_model.R`

**Model origin:** Initial-value support for the baseline spatial CMR model.

**What it does in plain terms:**

Works out sensible starting values for latent activity centres and alive/dead states so the NIMBLE model can initialise without missing values or impossible starting locations.

**Role in the project:** Useful debugging/helper script.

---

### `Scripts/05_ergon_gardner_spatial_cmr/development_and_debugging/04_initial_values_for_sex_group_badger_model.R`

**Model origin:** Initial-value support for the sex-effect spatial CMR model.

**What it does in plain terms:**

Same purpose as the baseline initial-value script, but adapted for the sex-effect model variant.

**Role in the project:** Useful if revisiting sex-effect models.

---

## Spatial CJS prototypes

These scripts explore spatial CJS-style alternatives. Unlike the Ergon & Gardner robust-design SCR model, these are more directly connected to CJS survival modelling and spatial/discrete locations.

### `Scripts/06_spatial_cjs_prototypes/01_spatial_cjs_badgers_initial.R`

**Model origin:** Spatial CJS model with irregular study area.

**Likely reference family:** Cormack-Jolly-Seber survival models extended with spatial displacement.

**What it does in plain terms:**

Builds a spatial capture-history object for badgers and tests a spatial CJS-style model where individuals survive, are detected, and have locations through time.

**Role in the project:** Prototype.

---

### `Scripts/06_spatial_cjs_prototypes/02_spatial_cjs_badgers_social_group.R`

**Model origin:** Discrete social-group/multistate CJS prototype.

**What it does in plain terms:**

Represents movement among social groups as transitions between discrete states. It was likely used to test whether social-group changes could be modelled directly.

**Role in the project:** Prototype.

---

### `Scripts/06_spatial_cjs_prototypes/03_spatial_cjs_badgers_first_location.R`

**Model origin:** Spatial CJS prototype with first-location treatment.

**What it does in plain terms:**

Similar to the initial spatial CJS prototype, but focuses on how to use or condition on the first observed location of each individual.

**Role in the project:** Prototype that likely informed the later “adjusted first location” Ergon & Gardner model.

---

### `Scripts/06_spatial_cjs_prototypes/04_spatial_cjs_badgers_full_social_group.R`

**Model origin:** Full social-group multistate CJS prototype.

**What it does in plain terms:**

Builds a larger multistate model where social groups are the possible alive states and death is the absorbing state. It estimates transition probabilities among social groups and death.

**Role in the project:** Prototype/reference for discrete social-group movement modelling.

---

## Plotting and animation scripts

### `Scripts/07_plotting_and_animation/01_estimated_prevalence_plot.R`

**Model origin:** Not a model. Data-summary/plotting script.

**What it does in plain terms:**

Summarises infection-test prevalence through time and plots estimated or observed prevalence by occasion/year.

**Role in the project:** Useful descriptive output.

---

### `Scripts/07_plotting_and_animation/02_movement_animation_baseline.R`

**Model origin:** Output visualisation for the baseline spatial CMR model.

**What it does in plain terms:**

Reads posterior samples for latent activity centres `S[i, coordinate, year]`, summarises them, and animates estimated movement through time.

**Role in the project:** Useful after running the baseline model with latent `S` nodes monitored.

---

### `Scripts/07_plotting_and_animation/03_movement_animation_infection_group.R`

**Model origin:** Output visualisation for the infection-group spatial CMR model.

**What it does in plain terms:**

Creates movement animations from the posterior activity-centre samples, colouring or grouping individuals by infection-history status.

**Role in the project:** Useful for communicating model outputs, especially movement differences between never-positive and test-positive-as-cub badgers.

---

# Suggested simplification of the folder structure

The current folder is much clearer than the original, but I would still split it into active, reference and archive material.

Suggested structure:

```text
BadgersLargeGrant/
├── README.md
├── BadgersLargeGrant.Rproj
├── Data/
│   ├── processed/
│   ├── spatial_original/
│   ├── spatial_model_grid/
│   └── derived_objects/
├── Scripts/
│   ├── 01_data_preparation/
│   ├── 02_spatial_setup/
│   ├── 03_main_models/
│   ├── 04_diagnostics/
│   ├── 05_plotting_and_animation/
│   └── archive_model_development/
├── Figures_and_animations/
├── Model_outputs/
└── Project_overview/
```

A particularly useful active-script set would be:

```text
Scripts/01_data_preparation/01_merge_clean_capture_diagnostic_data.R
Scripts/02_spatial_setup/01_create_spatial_grid_objects.R
Scripts/02_spatial_setup/02_plot_spatial_grid_and_setts.R
Scripts/03_main_models/01_spatial_cmr_infection_group_landscape_model.R
Scripts/04_diagnostics/01_check_model_output.R
Scripts/05_plotting_and_animation/03_movement_animation_infection_group.R
```

Everything else can go into:

```text
Scripts/archive_model_development/
```

That would make it much easier to know what to run and what is historical background.

---

# Reference/model-origin notes

These are the main model families represented in the folder:

- **Cormack-Jolly-Seber models:** survival and recapture models for marked individuals. Classic origins: Cormack (1964), Jolly (1965), Seber (1965).
- **Multistate/multisite capture-recapture models:** models where living animals can occupy different states/sites and transition among them. Associated model family includes Arnason-style and Brownie/Schwarz-style multistate extensions.
- **Temporary-emigration / robust-design models:** models where animals may be alive but temporarily unavailable for capture. Associated with Pollock’s robust design and later robust-design CMR developments.
- **Ergon & Gardner-style spatial robust-design CMR:** spatial capture-recapture model with latent activity centres, movement/dispersal, detection as a function of distance, and survival. This is the model family behind the active badger spatial CMR scripts.
- **NIMBLE implementation:** most active model scripts are written in NIMBLE, which allows Bayesian hierarchical models with custom latent-state structures.

