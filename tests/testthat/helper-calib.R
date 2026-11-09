# Shared fixtures for the calibration likelihood tests.

# Importance weights of 1 for every year from 1950 to 2019.
flat_wts <- function() stats::setNames(rep(1, 70), 1950:2019)

# California's calibration targets. find_loc_file() warns that it picked the
# latest of several calibration data files.
ca_targets <- function() suppressWarnings(calib_targets("CA"))

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
