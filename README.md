# MITUS
Modeling Interventions of Tuberculosis in the United States 

### Project Goals
The goal of the MITUS package is to provide a deterministic model to investigate the impacts of various TB control strategies in the United States. The model has been calibrated to national data and to data for each state and the District of Columbia. Built-in control strategies are in three broad categories: 
- pre-defined scenarios
- programmatic changes to the LTBI and active TB diagnosis and care cascades
- targeted testing and treatment of customizable risk groups

This fork of [PPML/MITUS](https://github.com/PPML/MITUS) calibrates the model to states and counties from recent data.

## Installation 
```
# install.packages("devtools")
devtools::install()  # from a clone of this repository
```

## Model data
The model uses two kinds of data for each location.

- **Model inputs** describe the population and its history: mortality, births, immigration, treatment and the initial population. `model_load` turns them into complete series from 1950 to 2100. Values for years the data cover are observations, or estimates derived from them; the years after an input's last observation are forecast by the rules at the top of `R/input_series.R`.
- **Calibration targets** are the outcomes the model is fitted to: TB cases, TB deaths, population and LTBI prevalence, and published studies of TB natural history. Calibration compares a simulation with them; the simulation does not use them.

`inst/README.md` describes the data files and where they are kept.

## Calibrating a state or county
`model_load` loads a location's model inputs. A state's inputs are in the package; another location's are read from `data_dir`. `last_weight_year` sets the last year of data the calibration fits.

```
model_load("CA")
model_load("SanDiego", data_dir = "path/to/county/inputs", last_weight_year = 2025)
```

`optim_b_st` fits the model's parameters to the location's calibration targets, starting from a row of the default starting values. It saves each round's `optim` result to the working directory.

```
optim_b_st(StartVal_st, samp_i = 1, loc = "CA")
```

`calib_fit_for` simulates a fitted parameter set and pairs each calibration target with its simulated value, and `plot_calib_fit` draws that table to a PDF.

```
load("Opt_CA_r9_1_2026-11-09.rda")  # o9, the last round
targets <- calib_targets("CA")
fit <- calib_fit_for(o9$par, "CA", targets)
plot_calib_fit(fit, "CA_fit.pdf", title = "California")
```
