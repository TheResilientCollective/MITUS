# Calibration targets for a single location. The calibration likelihood
# compares one location's simulation with one set of targets, built here from
# the calibration data files.
#
# The targets are a list of three tables:
#   obs        observations, one row per target, period, nativity and group,
#              with columns target, start, end, nativity, group, value, total
#   estimates  single estimates with intervals or fractions, with columns
#              target, start, end, group, estimate, lower, upper
#   studies    the published study data every location is calibrated against
# Each row's start and end give the calendar years its value covers.

# Location targets in the calibration data files, named as in the source list,
# with the field each comes from in the files.
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

# Column names, or element names for vectors, that the conversion reads from
# each source table.
target_required_names <- list(
  cases_yr = c("year", "cases"),
  cases_yr_ag_nat = c("year", "total"),
  cases_yr_ag_nat_5yr = c("usb", "yr_txt"),
  cases_nat_5yr = "usb",
  hr_cases_sm = c("ss", "fract_sm", "yr_txt"),
  rt_fb_cases_sm = c("p_recent", "yr_txt"),
  rct_cases_sm = "fract_sm",
  tbdeaths = c("Year", "Deaths"),
  pop_50_10 = c("year", "usb", "pop"),
  pop_00_19 = c("age_group", "usb"),
  tx_outcomes = c("year", "pct_disc", "pct_dead", "sample_size"),
  TLTBI_dist = c("FB", "HR"),
  LTBI_prev_US_11_IGRA = c("age_group", "SSp", "SSn"),
  LTBI_prev_FB_11_IGRA = c("age_group", "SSp", "SSn"),
  tbdeaths_age_yr = "year",
  borgdorff_data = c("Years", "mean"),
  ferebee_data = "At_Risk",
  sutherland_data = "At_Risk"
)

# Years of targets whose files do not record them: the national case trend
# runs 1953-1993, LTBI treatment starts are estimated for 2002, the homeless
# population for 2010, LTBI prevalence comes from the 2011 survey, and recent
# transmission is compared over 2015-2018.
undated_target_years <- list(
  cases_decline = 1953:1993,
  ltbi_tx = 2002,
  homeless_pop = 2010,
  ltbi_prev = 2011,
  recent_transmission = 2015:2018
)

# Targets whose observations are counts over multi-year buckets.
bucketed_targets <- c("cases_nativity", "cases_age", "cases_homeless", "cases_recent_entry")

#' Calibration targets for a state or county
#'
#' Reads the location's calibration data files and converts its calibration
#' targets to long tables. A state's own targets are picked out of the
#' all-state calibration data. A county's come from its own calibration data,
#' which may also replace national targets.
#'
#' @param loc state USPS code, or county code from \code{county_ID.csv}
#' @param data_dir directory searched for data files before the package's
#'   \code{inst/}; defaults to the one given to \code{model_load}
#' @param calib_data_file all-state calibration data file to use in place of
#'   the one \code{find_loc_file} picks
#' @return list of three tables, each row giving the calendar years its value
#'   covers in \code{start} and \code{end}: \code{obs}, the observations, with
#'   columns target, start, end, nativity, group, value and total;
#'   \code{estimates}, single estimates with their intervals, with columns
#'   target, start, end, group, estimate, lower and upper; and \code{studies},
#'   the published study data every location is calibrated against
#' @export
calib_targets <- function(loc, data_dir = get_data_dir(), calib_data_file = NULL) {
  loc_info <- resolve_location(loc)
  if (is.null(loc_info) || loc_info$loc_type == "national") {
    stop("calib_targets() takes a state or county code, not '", loc, "'", call. = FALSE)
  }
  all_states <- function(fallback_loc) {
    if (!is.null(calib_data_file)) {
      return(readRDS(calib_data_file))
    }
    readRDS(find_loc_file(fallback_loc, "CalibDat", fallback_loc = "ST", data_dir = data_dir))
  }
  if (loc_info$loc_type == "county") {
    return(calib_targets_county(
      county_dat = readRDS(find_loc_file(loc, "CalibDat", data_dir = data_dir)),
      deaths_by_age = readRDS(find_loc_file(loc, "deathbyAge", data_dir = data_dir)),
      national_dat = all_states("ST"),
      loc = loc
    ))
  }
  calib_targets_state(
    calib_dat = all_states(loc),
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
#' @return list of the tables \code{obs}, \code{estimates} and \code{studies}
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
  source <- c(
    own,
    list(deaths_by_age = deaths_by_age[[st]]),
    calib_dat[national_target_fields],
    calib_dat[shared_target_fields]
  )
  long_targets(check_target_source(source, state_name), state_name)
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
#' @return list of the tables \code{obs}, \code{estimates} and \code{studies}
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
  source <- c(
    own,
    list(deaths_by_age = deaths_by_age),
    national,
    national_dat[shared_target_fields]
  )
  long_targets(check_target_source(source, loc), loc)
}

#' Names a source table's values are read by: column names for a table or
#' array, element names for a vector
#' @noRd
target_names <- function(x) if (is.null(dim(x))) names(x) else colnames(x)

#' Check that one location's source targets are complete and have the columns
#' the conversion reads
#'
#' @param source list of one location's targets as laid out in the calibration
#'   data files
#' @param label location name or code, used in errors
#' @return \code{source}, unchanged
#' @noRd
check_target_source <- function(source, label) {
  fail <- function(...) stop("Calibration targets for ", label, " ", ..., call. = FALSE)
  expected <- c(
    names(location_target_fields), "deaths_by_age", national_target_fields,
    shared_target_fields
  )
  missing <- expected[vapply(expected, function(field) is.null(source[[field]]), logical(1))]
  if (length(missing) > 0) {
    fail("are missing ", toString(missing))
  }
  tables <- names(location_table_name_cols)
  empty <- tables[vapply(tables, function(field) nrow(source[[field]]) == 0, logical(1))]
  if (length(empty) > 0) {
    fail("have no rows in ", toString(empty))
  }
  for (target in names(target_required_names)) {
    absent <- setdiff(target_required_names[[target]], target_names(source[[target]]))
    if (length(absent) > 0) {
      fail("have no ", toString(absent), " in ", target)
    }
  }
  if (nrow(source$rct_cases_sm) != 1) {
    fail("need one row in rct_cases_sm")
  }
  source
}

#' Calendar years named by period labels
#'
#' Reads labels such as \code{"X1995.1999"}, \code{"95-99"} or \code{"96_00"},
#' taking two-digit years from 50 onward as 19xx and earlier ones as 20xx.
#'
#' @param labels period labels
#' @return integer matrix with columns \code{start} and \code{end}
#' @noRd
periods_from_labels <- function(labels) {
  years <- regmatches(labels, gregexpr("[0-9]+", labels))
  bad <- lengths(years) != 2
  if (any(bad)) {
    stop("Period labels must name two years: ", toString(labels[bad]), call. = FALSE)
  }
  to_year <- function(y) {
    y <- as.integer(y)
    ifelse(y < 100, ifelse(y >= 50, 1900L, 2000L) + y, y)
  }
  periods <- t(vapply(years, to_year, integer(2)))
  colnames(periods) <- c("start", "end")
  backwards <- periods[, "end"] < periods[, "start"]
  if (any(backwards)) {
    stop("Period labels end before they start: ", toString(labels[backwards]), call. = FALSE)
  }
  periods
}

#' Rows of the observations table
#' @noRd
obs_rows <- function(target, start, end = start, value, nativity = "all", group = "all",
                     total = NA_real_) {
  data.frame(
    target = target, start = as_year(start), end = as_year(end), nativity = nativity,
    group = group, value = as.numeric(value), total = as.numeric(total),
    stringsAsFactors = FALSE
  )
}

#' Rows of the estimates table
#' @noRd
estimate_rows <- function(target, start, end = start, estimate, group = "all",
                          lower = NA_real_, upper = NA_real_) {
  data.frame(
    target = target, start = as_year(start), end = as_year(end), group = group,
    estimate = as.numeric(estimate), lower = as.numeric(lower), upper = as.numeric(upper),
    stringsAsFactors = FALSE
  )
}

# Nativity groups and their codes in the calibration data files.
nativity_codes <- c(usb = 1, nusb = 0)

#' Convert one location's source targets to the long target tables
#'
#' @param source list of one location's targets, checked by
#'   \code{check_target_source}
#' @param label location name or code, used in errors
#' @return list of the tables \code{obs}, \code{estimates} and \code{studies}
#' @noRd
long_targets <- function(source, label) {
  obs <- rbind(
    obs_rows("cases_total", source$cases_yr[, "year"], value = source$cases_yr[, "cases"]),
    case_decline_rows(source, label),
    nativity_case_rows(source$cases_nat_5yr),
    age_case_rows(source$cases_yr_ag_nat_5yr, label),
    homeless_case_rows(source$hr_cases_sm),
    recent_entry_rows(source, label),
    tx_outcome_rows(source$tx_outcomes),
    ltbi_prev_rows(source),
    obs_rows("tb_deaths_total", source$tbdeaths[, "Year"], value = source$tbdeaths[, "Deaths"]),
    tb_death_age_rows(source$tbdeaths_age_yr, label),
    death_rows(source$deaths_by_age, label),
    census_rows(source$pop_50_10),
    population_rows(source$pop_00_19)
  )
  rownames(obs) <- NULL
  volume <- source$TLTBI_volume
  homeless <- source$homeless_pop
  window <- undated_target_years$recent_transmission
  estimates <- rbind(
    estimate_rows("ltbi_tx_total", undated_target_years$ltbi_tx, estimate = volume[1],
                  lower = volume[2], upper = volume[3]),
    estimate_rows("ltbi_tx_share", undated_target_years$ltbi_tx, group = c("nusb", "homeless"),
                  estimate = source$TLTBI_dist[c("FB", "HR")]),
    estimate_rows("homeless_pop", undated_target_years$homeless_pop, estimate = homeless[1],
                  lower = homeless[2], upper = homeless[3]),
    estimate_rows("recent_transmission", min(window), max(window),
                  estimate = source$rct_cases_sm[, "fract_sm"])
  )
  rownames(estimates) <- NULL
  check_long_targets(
    list(obs = obs, estimates = estimates, studies = source[shared_target_fields]), label
  )
}

#' Expected cases over 1953-1993: the national year-to-year changes, scaled to
#' end at the location's 1993 cases
#' @noRd
case_decline_rows <- function(source, label) {
  years <- undated_target_years$cases_decline
  d <- source$cases_yr
  cases_end <- as.numeric(d[as_year(d[, "year"]) == max(years), "cases"])
  if (length(cases_end) != 1) {
    stop("Calibration targets for ", label, " have no ", max(years), " row in cases_yr",
         call. = FALSE)
  }
  change <- source$cases_prop_change_53_94
  obs_rows("cases_decline", years, value = cumprod(change) / prod(change) * cases_end)
}

#' Cases by nativity in 5-year buckets
#' @noRd
nativity_case_rows <- function(d) {
  cols <- case_bucket_cols(d)
  periods <- periods_from_labels(cols)
  rows <- lapply(names(nativity_codes), function(nativity) {
    counts <- as.numeric(unlist(d[d[, "usb"] == nativity_codes[[nativity]], cols]))
    obs_rows("cases_nativity", periods[, "start"], periods[, "end"], counts, nativity = nativity)
  })
  do.call(rbind, rows)
}

#' Cases by nativity and age in 5-year buckets, NA where suppressed
#' @noRd
age_case_rows <- function(d, label) {
  ages <- case_age_cols(d)
  if (length(ages) != length(data_age_groups)) {
    stop("Calibration targets for ", label, " need ", length(data_age_groups),
         " age groups in cases_yr_ag_nat_5yr", call. = FALSE)
  }
  rows <- lapply(names(nativity_codes), function(nativity) {
    group <- d[d[, "usb"] == nativity_codes[[nativity]], , drop = FALSE]
    periods <- periods_from_labels(as.character(group[, "yr_txt"]))
    n_ages <- length(ages)
    obs_rows(
      "cases_age",
      start = rep(periods[, "start"], each = n_ages),
      end = rep(periods[, "end"], each = n_ages),
      value = as.vector(t(as.matrix(group[, ages]))),
      nativity = nativity,
      group = rep(data_age_groups, times = nrow(group))
    )
  })
  do.call(rbind, rows)
}

#' Cases among people experiencing homelessness in 5-year buckets: the
#' smoothed fraction times all cases
#' @noRd
homeless_case_rows <- function(d) {
  periods <- periods_from_labels(as.character(d[, "yr_txt"]))
  obs_rows("cases_homeless", periods[, "start"], periods[, "end"],
           value = d[, "fract_sm"] * d[, "ss"], total = d[, "ss"])
}

#' Non-US-born cases among recent arrivals in 5-year buckets: the smoothed
#' fraction of recent arrivals times the bucket's non-US-born cases
#' @noRd
recent_entry_rows <- function(source, label) {
  d <- source$rt_fb_cases_sm
  periods <- periods_from_labels(as.character(d[, "yr_txt"]))
  annual <- source$cases_yr_ag_nat
  years <- as_year(annual[, "year", "nusb"])
  nusb <- vapply(seq_len(nrow(periods)), function(i) {
    bucket <- periods[i, "start"]:periods[i, "end"]
    if (!all(bucket %in% years)) {
      stop("Calibration targets for ", label, " lack annual non-US-born cases for ",
           toString(setdiff(bucket, years)), call. = FALSE)
    }
    sum(annual[years %in% bucket, "total", "nusb"])
  }, numeric(1))
  obs_rows("cases_recent_entry", periods[, "start"], periods[, "end"],
           value = d[, "p_recent"] * nusb, nativity = "nusb", total = nusb)
}

# Treatment outcomes in the order the model reports them.
tx_outcome_groups <- c("completed", "discontinued", "died")

#' Treatment outcomes each year: completed, discontinued and died
#' @noRd
tx_outcome_rows <- function(d) {
  share <- cbind(1 - rowSums(d[, c("pct_disc", "pct_dead")]), d[, "pct_disc"], d[, "pct_dead"])
  counts <- share * d[, "sample_size"]
  n_outcomes <- length(tx_outcome_groups)
  obs_rows(
    "tx_outcomes",
    start = rep(d[, "year"], each = n_outcomes),
    value = as.vector(t(counts)),
    group = rep(tx_outcome_groups, times = nrow(d)),
    total = rep(d[, "sample_size"], each = n_outcomes)
  )
}

#' IGRA-positive survey counts by nativity and survey age group
#' @noRd
ltbi_prev_rows <- function(source) {
  surveys <- list(usb = source$LTBI_prev_US_11_IGRA, nusb = source$LTBI_prev_FB_11_IGRA)
  rows <- lapply(names(surveys), function(nativity) {
    d <- surveys[[nativity]]
    obs_rows("ltbi_prev", undated_target_years$ltbi_prev, value = d[, "SSp"],
             nativity = nativity, group = as.character(d[, "age_group"]),
             total = d[, "SSp"] + d[, "SSn"])
  })
  do.call(rbind, rows)
}

#' Deaths with TB by age group each year
#' @noRd
tb_death_age_rows <- function(d, label) {
  ages <- setdiff(colnames(d), "year")
  if (length(ages) != length(data_age_groups)) {
    stop("Calibration targets for ", label, " need ", length(data_age_groups),
         " age groups in tbdeaths_age_yr", call. = FALSE)
  }
  n_ages <- length(ages)
  obs_rows(
    "tb_deaths_age",
    start = rep(d[, "year"], each = n_ages),
    value = as.vector(t(as.matrix(d[, ages]))),
    group = rep(data_age_groups, times = nrow(d))
  )
}

#' All-cause deaths each year, in total and by age group
#'
#' The deaths table has one row per year: the year, deaths in the data's ten
#' age groups, then the total.
#' @noRd
death_rows <- function(d, label) {
  n_ages <- length(data_age_groups)
  if (ncol(d) != n_ages + 2) {
    stop("Calibration targets for ", label, " need the year, ", n_ages,
         " age groups and the total in the deaths by age table", call. = FALSE)
  }
  ages <- d[, 1 + seq_len(n_ages), drop = FALSE]
  rbind(
    obs_rows("deaths_total", d[, 1], value = d[, ncol(d)]),
    obs_rows(
      "deaths_age",
      start = rep(d[, 1], each = n_ages),
      value = as.vector(t(ages)),
      group = rep(data_age_groups, times = nrow(d))
    )
  )
}

#' Census populations by nativity
#' @noRd
census_rows <- function(d) {
  nativity <- names(nativity_codes)[match(d[, "usb"], nativity_codes)]
  obs_rows("population_census", d[, "year"], value = d[, "pop"], nativity = nativity)
}

#' Annual populations by nativity and age group
#'
#' The annual table follows the model's 11 age groups but puts everyone 85 and
#' over in group 10, leaving group 11 empty; group 11 is dropped.
#' @noRd
population_rows <- function(d) {
  years <- suppressWarnings(as.integer(colnames(d)))
  year_cols <- which(!is.na(years))
  d <- d[d[, "age_group"] %in% seq_along(data_age_groups), , drop = FALSE]
  rows <- lapply(year_cols, function(col) {
    obs_rows(
      "population", years[col], value = d[, col],
      nativity = names(nativity_codes)[match(d[, "usb"], nativity_codes)],
      group = data_age_groups[d[, "age_group"]]
    )
  })
  do.call(rbind, rows)
}

#' Check the long target tables
#'
#' Stops with an error when a target has no rows, or when a bucketed target's
#' buckets overlap, leave gaps or differ in length, or when the age buckets do
#' not match the nativity buckets.
#'
#' @param targets list of the tables \code{obs}, \code{estimates} and
#'   \code{studies}
#' @param label location name or code, used in errors
#' @return \code{targets}, unchanged
#' @noRd
check_long_targets <- function(targets, label) {
  fail <- function(...) stop("Calibration targets for ", label, " ", ..., call. = FALSE)
  needed_obs <- c(
    "cases_total", "cases_decline", bucketed_targets, "tx_outcomes", "ltbi_prev",
    "tb_deaths_total", "tb_deaths_age", "deaths_total", "deaths_age", "population_census",
    "population"
  )
  absent <- setdiff(needed_obs, targets$obs$target)
  if (length(absent) > 0) {
    fail("have no rows for ", toString(absent))
  }
  for (target in bucketed_targets) {
    rows <- targets$obs[targets$obs$target == target, , drop = FALSE]
    for (nativity in unique(rows$nativity)) {
      periods <- target_periods(rows[rows$nativity == nativity, , drop = FALSE])
      if (!periods_are_buckets(periods)) {
        fail("have ", target, " periods that overlap, leave gaps or differ in length: ",
             toString(paste(periods$start, periods$end, sep = "-")))
      }
    }
  }
  for (nativity in names(nativity_codes)) {
    of <- function(target) {
      target_periods(targets$obs[targets$obs$target == target &
                                   targets$obs$nativity == nativity, , drop = FALSE])
    }
    if (!identical(of("cases_age"), of("cases_nativity"))) {
      fail("have cases_age periods that differ from cases_nativity periods for ", nativity)
    }
  }
  targets
}

#' Distinct periods of target rows, earliest first
#' @noRd
target_periods <- function(rows) {
  periods <- unique(rows[, c("start", "end")])
  periods <- periods[order(periods$start, periods$end), , drop = FALSE]
  rownames(periods) <- NULL
  periods
}

#' Whether periods are contiguous, non-overlapping and of equal length
#' @noRd
periods_are_buckets <- function(periods) {
  n <- nrow(periods)
  lengths_match <- length(unique(periods$end - periods$start)) == 1
  n == 1 || (lengths_match && all(periods$start[-1] == periods$end[-n] + 1))
}

#' Stop when a multi-year target period runs past the last weighted year
#'
#' The calibration simulates through the last weighted year, so a bucket that
#' ends later would be compared with a partial simulation.
#'
#' @param targets list of the tables \code{obs}, \code{estimates} and
#'   \code{studies}
#' @param wts importance weights by year
#' @noRd
check_target_horizon <- function(targets, wts) {
  last_year <- max(as.integer(names(wts)))
  rows <- rbind(targets$obs[, c("target", "start", "end")],
                targets$estimates[, c("target", "start", "end")])
  late <- rows[rows$end > rows$start & rows$end > last_year, , drop = FALSE]
  if (nrow(late) > 0) {
    stop("Calibration targets run past the last weighted year, ", last_year, ", in ",
         toString(unique(late$target)), call. = FALSE)
  }
}

#' Observation rows of one target, optionally for one nativity group and some
#' groups, earliest first
#' @noRd
target_obs <- function(targets, target, nativity = NULL, group = NULL) {
  d <- targets$obs[targets$obs$target == target, , drop = FALSE]
  if (!is.null(nativity)) d <- d[d$nativity == nativity, , drop = FALSE]
  if (!is.null(group)) d <- d[d$group %in% group, , drop = FALSE]
  if (nrow(d) == 0) {
    stop("No calibration target rows for ", target, call. = FALSE)
  }
  d[order(d$start), , drop = FALSE]
}

#' Estimate rows of one target, optionally for some groups
#' @noRd
target_estimate <- function(targets, target, group = NULL) {
  d <- targets$estimates[targets$estimates$target == target, , drop = FALSE]
  if (!is.null(group)) d <- d[match(group, d$group), , drop = FALSE]
  if (nrow(d) == 0 || anyNA(d$target)) {
    stop("No calibration estimate for ", target, call. = FALSE)
  }
  d
}
