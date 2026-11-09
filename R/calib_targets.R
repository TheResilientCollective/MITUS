# Calibration targets for a single location. The calibration likelihood
# compares one location's simulation with one list of targets, built here from
# the calibration data files.

# Targets each location has its own values for, named as in the targets list,
# with the field each comes from in the calibration data files.
location_target_fields <- c(
  cases_yr = "cases_yr_st",
  cases_yr_ag_nat = "cases_yr_ag_nat_st",
  cases_yr_ag_nat_5yr = "cases_yr_ag_nat_st_5yr",
  cases_nat_5yr = "cases_nat_st_5yr",
  hr_cases_sm = "hr_cases_sm",
  rt_fb_cases_sm = "rt_fb_cases_sm",
  rct_cases_sm = "rct_cases_sm",
  TLTBI_volume = "TLTBI_volume_state",
  tbdeaths = "tbdeaths",
  homeless_pop = "homeless_pop",
  pop_50_10 = "pop_50_10",
  pop_00_19 = "pop_00_19"
)

# Tables in the all-state calibration data that hold rows for every state,
# with the column naming each row's state. Every other location field there is
# a list with one element per state.
state_table_name_cols <- c(
  cases_nat_st_5yr = "State",
  hr_cases_sm = "st",
  rt_fb_cases_sm = "st",
  rct_cases_sm = "st"
)

# National targets. A county's own calibration data may replace any of them.
national_target_fields <- c(
  "cases_prop_change_53_94",
  "tx_outcomes",
  "TLTBI_dist",
  "LTBI_prev_US_11_IGRA",
  "LTBI_prev_FB_11_IGRA",
  "tbdeaths_age_yr"
)

# Published study data every location is calibrated against.
shared_target_fields <- c("borgdorff_data", "ferebee_data", "sutherland_data")

#' Calibration targets for a state or county
#'
#' Reads the location's calibration data files and gathers its calibration
#' targets into one list: the location's own targets, the national targets and
#' the published study data. A state's own targets are picked out of the
#' all-state calibration data. A county's come from its own calibration data,
#' which may also replace national targets.
#'
#' @param loc state USPS code, or county code from \code{county_ID.csv}
#' @param data_dir directory searched for data files before the package's
#'   \code{inst/}; defaults to the one given to \code{model_load}
#' @return named list of calibration targets
#' @noRd
calib_targets <- function(loc, data_dir = get_data_dir()) {
  loc_info <- resolve_location(loc)
  if (is.null(loc_info) || loc_info$loc_type == "national") {
    stop("calib_targets() takes a state or county code, not '", loc, "'", call. = FALSE)
  }
  if (loc_info$loc_type == "county") {
    return(calib_targets_county(
      county_dat = readRDS(find_loc_file(loc, "CalibDat", data_dir = data_dir)),
      deaths_by_age = readRDS(find_loc_file(loc, "deathbyAge", data_dir = data_dir)),
      national_dat = readRDS(find_loc_file("ST", "CalibDat", data_dir = data_dir)),
      loc = loc
    ))
  }
  calib_targets_state(
    calib_dat = readRDS(
      find_loc_file(loc, "CalibDat", fallback_loc = "ST", data_dir = data_dir)
    ),
    deaths_by_age = readRDS(
      find_loc_file(loc, "deathbyAge", fallback_loc = "ST", data_dir = data_dir)
    ),
    st = loc_info$st,
    state_name = loc_info$name
  )
}

#' Calibration targets for one state, from the all-state calibration data
#'
#' @param calib_dat all-state calibration data, as read from
#'   \code{ST_CalibDat_*.rds}
#' @param deaths_by_age all-state all-cause deaths by age, a list with one
#'   matrix per state
#' @param st the state's position in \code{stateID}
#' @param state_name the state's name, as in \code{stateID}
#' @return named list of calibration targets
#' @noRd
calib_targets_state <- function(calib_dat, deaths_by_age, st, state_name) {
  own <- lapply(location_target_fields, function(field) {
    value <- calib_dat[[field]]
    if (field %in% names(state_table_name_cols)) {
      value[value[[state_table_name_cols[[field]]]] == state_name, , drop = FALSE]
    } else {
      value[[st]]
    }
  })
  targets <- c(
    own,
    list(deaths_by_age = deaths_by_age[[st]]),
    calib_dat[national_target_fields],
    calib_dat[shared_target_fields]
  )
  check_calib_targets(targets, state_name)
}

#' Calibration targets for one county
#'
#' The county's calibration data holds its own targets, laid out as one state's
#' part of the all-state calibration data. National targets the county does not
#' provide come from the all-state calibration data, and a message names them.
#'
#' @param county_dat the county's calibration data
#' @param deaths_by_age the county's all-cause deaths by age
#' @param national_dat all-state calibration data, for national targets and
#'   published study data
#' @param loc county code, used in messages
#' @return named list of calibration targets
#' @noRd
calib_targets_county <- function(county_dat, deaths_by_age, national_dat, loc) {
  own <- county_dat[location_target_fields]
  names(own) <- names(location_target_fields)
  for (field in names(state_table_name_cols)) {
    table <- county_dat[[field]]
    if (!is.null(table) && length(unique(table[[state_table_name_cols[[field]]]])) != 1) {
      stop("Calibration data for ", loc, " has rows for more than one location in ", field,
           call. = FALSE)
    }
  }
  provided <- names(Filter(Negate(is.null), county_dat))
  national <- national_dat[national_target_fields]
  from_county <- intersect(national_target_fields, provided)
  national[from_county] <- county_dat[from_county]
  fallback <- setdiff(national_target_fields, provided)
  if (length(fallback) > 0) {
    message("Calibrating ", loc, " against national targets for ", toString(fallback))
  }
  targets <- c(
    own,
    list(deaths_by_age = deaths_by_age),
    national,
    national_dat[shared_target_fields]
  )
  check_calib_targets(targets, loc)
}

#' Check that a targets list is complete
#'
#' @param targets named list of calibration targets
#' @param label location name or code, used in errors
#' @return \code{targets}, unchanged
#' @noRd
check_calib_targets <- function(targets, label) {
  expected <- c(
    names(location_target_fields), "deaths_by_age", national_target_fields,
    shared_target_fields
  )
  missing <- expected[vapply(expected, function(field) is.null(targets[[field]]), logical(1))]
  if (length(missing) > 0) {
    stop("Calibration targets for ", label, " are missing ", toString(missing), call. = FALSE)
  }
  tables <- names(location_target_fields)[location_target_fields %in% names(state_table_name_cols)]
  empty <- tables[vapply(tables, function(field) nrow(targets[[field]]) == 0, logical(1))]
  if (length(empty) > 0) {
    stop("Calibration targets for ", label, " have no rows in ", toString(empty), call. = FALSE)
  }
  if (nrow(targets$rct_cases_sm) != 1) {
    stop("Calibration targets for ", label, " need one row in rct_cases_sm", call. = FALSE)
  }
  targets
}
