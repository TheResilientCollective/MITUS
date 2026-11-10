# Shared fixtures for the calibration likelihood tests.

# Importance weights of 1 for every year from 1950 to 2019.
flat_wts <- function() stats::setNames(rep(1, 70), 1950:2019)

# The 10-29-21 all-state calibration data, with its 2015-2019 nativity bucket
# and its 1999-2003 recent-arrival bucket labelled with the years they cover.
state_calib_data <- function() {
  d <- readRDS(system.file("ST", "ST_CalibDat_10-29-21.rds", package = "MITUS"))
  names(d$cases_nat_st_5yr)[names(d$cases_nat_st_5yr) == "X2015.1019"] <- "X2015.2019"
  d$rt_fb_cases_sm$yr_txt <- sub("99_13", "99_03", d$rt_fb_cases_sm$yr_txt)
  d
}

# All-state all-cause deaths by age.
state_deaths_by_age <- function() readRDS(system.file("ST", "STdeathbyAge.rds", package = "MITUS"))

# California's calibration targets, from the corrected 10-29-21 data.
ca_targets <- function() {
  calib_targets_state(state_calib_data(), state_deaths_by_age(), 5, "California")
}

# Restore the global environment when the calling test ends: globals created
# during the test are removed and globals it changed get their old values back.
local_restore_globals <- function(env = parent.frame()) {
  g <- globalenv()
  saved <- as.list(g, all.names = TRUE)
  withr::defer(
    {
      rm(list = setdiff(ls(g, all.names = TRUE), names(saved)), envir = g)
      list2env(saved, envir = g)
    },
    envir = env
  )
}

# Set globals for the calling test only.
local_globals <- function(values, env = parent.frame()) {
  local_restore_globals(env)
  list2env(values, envir = globalenv())
}

# Run model_load() for the calling test only.
local_model_load <- function(loc, env = parent.frame()) {
  local_restore_globals(env)
  suppressWarnings(model_load(loc))
}

# Progression parameters with a different value in row 3, the row the
# published-study likelihoods use.
fake_progression <- function() {
  mpfast <- matrix(0.9, 11, 4)
  mpfast[3, ] <- 0.1
  mrslow <- matrix(0.01, 11, 4)
  mrslow[3, ] <- c(0.0001, 0.0002, 0.0003, 0.0004)
  list(Mpfast = mpfast, Mrslow = mrslow, rfast = 0.1, rRecov = 0.002)
}

# Annual simulation outputs for 1950-2019 with every named output, filled with
# positive values in which each part is smaller than its whole.
fake_outputs <- function() {
  set.seed(3)
  outputs <- func_ResNam()
  out <- matrix(stats::runif(70 * length(outputs), 1, 2) * 1e-3, 70, length(outputs),
                dimnames = list(NULL, outputs))
  parts <- grep("^(NOTIF_US_|N_US_LTBI_|N_FB_LTBI_|INCID_REC_|TLTBI_INITS_)", outputs, value = TRUE)
  out[, parts] <- out[, parts] * 0.2
  out
}
