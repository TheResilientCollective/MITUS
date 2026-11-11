# Data files

This folder holds the data the model uses, in one folder per location:

- `{loc}/`: one state or DC, such as `CA/`
- `ST/`: data shared by all states
- `US/`: national data
- `extdata/`: location lists and other tables

Call `model_load` with a location code to load that location's data. To load data kept outside the package, call `model_load` with the `data_dir` argument, naming a directory with the same layout.

## File names
Files are named `{loc}_{kind}`, often followed by a date, such as `CA_ModelInputs_11-12-21.rds`. The date records when the file was made, not the years its data cover. Where a folder holds several files of one kind, the latest is used.

## Files used when loading, calibrating or projecting a location

| Files | Contents | Used as |
|---|---|---|
| `{loc}/{loc}_ModelInputs_*`, `US/US_Inputs_*` | The location's model input object: background mortality, initial population, births, treatment default and quality, net migration, and immigration with the TB burden among immigrants | Model inputs |
| `ST/ST_SingleYearLifeTable` | Mortality rates by single year of age and year, one table per state | Model inputs: scales the national mortality rates to a state |
| `ST/ST_PopCountsByAge` | Population by single year of age at census years, for every state | Model inputs: aging rates |
| `US/US_MortalityCountsByAge`, `US/US_PopCountsByAge` | National deaths and population by single year of age and year | Model inputs: mortality weighting, and the US aging rates |
| `ST/ST_CalibDat_*` | Calibration targets for every state | Calibration targets |
| `US/US_CalibDat_*` | National calibration targets | Calibration targets |
| `ST/STdeathbyAge` | All-cause deaths by age group and year, for every state | Calibration targets |
| `ST/ST_ParamInit_*`, `US/US_ParamInit_*` | Every model parameter: its background value, its prior, and whether it is calibrated | Parameters |
| `ST/ST_StartVal_*`, `US/US_StartVal_*` | Starting points for calibration, on the calibration's normal scale | Calibration |
| `{loc}/{loc}_Optim_all_*`, `US/US_Optim_all_*` | Calibration results: one parameter set per starting point on the calibration's normal scale, with its negative log posterior in the last column | Calibration results |
| `{loc}/{loc}_Param_*`, `US/US_Param_*` | Calibrated parameter sets in model units | Parameters for projection |
| `{loc}/{loc}_results_1`, `US/US_results_1` | Base-case model outputs | Targeted testing and treatment scenarios |
| `extdata/county_ID.csv` | Locations other than states, with their codes, FIPS codes and states' FIPS codes | Locations |

## Data for a new location
A location without a folder here needs a folder named by its location code, `{loc}/`, holding:

| File | Contents | Needed to |
|---|---|---|
| `{loc}_ModelInputs` | The model input object, laid out as a state's | Load the location |
| `{loc}_SingleYearLifeTable` | Mortality rates by single year of age and year, for years the national counts in `US/` cover | Load the location |
| `{loc}_PopCountsByAge` | Population by single year of age at census years | Load the location |
| `{loc}_CalibDat` | Calibration targets, laid out as one state's entry in `ST/ST_CalibDat_*` | Calibrate |
| `{loc}_deathbyAge` | All-cause deaths by age group and year | Calibrate |
| `{loc}_Param` | A calibrated parameter set in model units | Project without calibrating |

Yearly values in the model input object start in 1950, labelled by year or given one per year. They may end before 2100; the years after their last observation are forecast.

A location that is not a state is also listed in `extdata/county_ID.csv`. Parameter definitions and calibration starting points come from `ST/` unless the folder holds its own `ParamInit` and `StartVal` files.

## Other files
Other files in these folders are not used when loading, calibrating or projecting a location. They include:

- `{loc}/calibration_targets/`, `{loc}/calibration_outputs/` and `ST/calibration_targets/`: files for the Tabby2 web application's calibration plots
- `ST/ST_immigration`, `ST/ST_netmigration`, `ST/ST_decennialpop` and `ST/ST_TBburdenimmig`: sources for building the state model input files
