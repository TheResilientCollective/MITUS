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

# Location targets stored as tables with rows for every state in the all-state
# calibration data, with the column that names each row's location. Every other
# location target is stored there as a list with one element per state.
location_table_name_cols <- c(
  cases_nat_5yr = "State",
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

# Column names, or element names for vectors, that the likelihood reads from
# each target.
target_required_names <- list(
  cases_yr = c("year", "cases"),
  cases_yr_ag_nat = c("year", "total"),
  cases_yr_ag_nat_5yr = "usb",
  cases_nat_5yr = "usb",
  hr_cases_sm = c("ss", "fract_sm"),
  rt_fb_cases_sm = "p_recent",
  rct_cases_sm = "fract_sm",
  tbdeaths = c("Year", "Deaths"),
  pop_50_10 = c("year", "usb", "pop"),
  pop_00_19 = c("age_group", "usb"),
  tx_outcomes = c("year", "pct_disc", "pct_dead", "sample_size"),
  TLTBI_dist = c("FB", "HR"),
  LTBI_prev_US_11_IGRA = c("SSp", "SSn"),
  LTBI_prev_FB_11_IGRA = c("SSp", "SSn"),
  tbdeaths_age_yr = "year",
  borgdorff_data = c("Years", "mean"),
  ferebee_data = "At_Risk",
  sutherland_data = "At_Risk"
)

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
  own <- lapply(names(location_target_fields), function(target) {
    value <- calib_dat[[location_target_fields[[target]]]]
    if (target %in% names(location_table_name_cols)) {
      value[value[, location_table_name_cols[[target]]] == state_name, , drop = FALSE]
    } else {
      value[[st]]
    }
  })
  names(own) <- names(location_target_fields)
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
  for (target in names(location_table_name_cols)) {
    table <- own[[target]]
    name_col <- location_table_name_cols[[target]]
    if (is.null(table)) next
    if (!name_col %in% colnames(table)) {
      stop("Calibration data for ", loc, " has no ", name_col, " column in ",
           location_target_fields[[target]], call. = FALSE)
    }
    if (length(unique(table[, name_col])) > 1) {
      stop("Calibration data for ", loc, " has rows for more than one location in ",
           location_target_fields[[target]], call. = FALSE)
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

#' Names a target's values are read by: column names for a table or array,
#' element names for a vector
#' @noRd
target_names <- function(x) if (is.null(dim(x))) names(x) else colnames(x)

#' Check that a targets list is complete and laid out as the likelihood reads it
#'
#' Stops with an error naming the problem when a target is missing, a table has
#' no rows, a required column is missing, the 5-year case tables do not have
#' one bucket for each simulated bucket, or the annual cases lack 1993.
#'
#' @param targets named list of calibration targets
#' @param label location name or code, used in errors
#' @return \code{targets}, unchanged
#' @noRd
check_calib_targets <- function(targets, label) {
  fail <- function(...) stop("Calibration targets for ", label, " ", ..., call. = FALSE)
  expected <- c(
    names(location_target_fields), "deaths_by_age", national_target_fields,
    shared_target_fields
  )
  missing <- expected[vapply(expected, function(field) is.null(targets[[field]]), logical(1))]
  if (length(missing) > 0) {
    fail("are missing ", toString(missing))
  }
  tables <- names(location_table_name_cols)
  empty <- tables[vapply(tables, function(field) nrow(targets[[field]]) == 0, logical(1))]
  if (length(empty) > 0) {
    fail("have no rows in ", toString(empty))
  }
  for (target in names(target_required_names)) {
    absent <- setdiff(target_required_names[[target]], target_names(targets[[target]]))
    if (length(absent) > 0) {
      fail("have no ", toString(absent), " in ", target)
    }
  }
  if (nrow(targets$rct_cases_sm) != 1) {
    fail("need one row in rct_cases_sm")
  }
  n_buckets <- length(case_bucket_years)
  if (length(case_bucket_cols(targets$cases_nat_5yr)) != n_buckets) {
    fail("need ", n_buckets, " 5-year buckets in cases_nat_5yr")
  }
  # The age table has one row per bucket for each of the two nativity groups
  bucket_rows <- c(
    cases_yr_ag_nat_5yr = 2 * n_buckets, hr_cases_sm = n_buckets, rt_fb_cases_sm = n_buckets
  )
  for (target in names(bucket_rows)) {
    if (nrow(targets[[target]]) != bucket_rows[[target]]) {
      fail("need ", n_buckets, " 5-year buckets in ", target)
    }
  }
  if (!1993 %in% as_year(targets$cases_yr[, "year"])) {
    fail("have no 1993 row in cases_yr")
  }
  targets
}
