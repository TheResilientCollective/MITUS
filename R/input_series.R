# Input series for one location: the model inputs as complete series from 1950
# to 2100, built once per location by model_load. None of them depend on the
# model parameters; param_init applies the parameters to them.
#
# Each input is read by its years. A series is observed through its last year,
# and the years after it are forecast here. An input that runs through 2100
# holds its own forecast after its data, so it is taken as observed through the
# year in legacy_observed_through.

# First and last year of the input series.
input_series_start_year <- 1950
input_series_end_year <- 2100

# Last observed year of an input that runs through 2100. These are the years the
# previous code treated as the end of the data: it forecast mortality from 2017,
# and applied the 2020 pandemic reduction, a correction to forecast immigration,
# from March 2020.
legacy_observed_through <- c(mortality = 2016, immigration = 2019)

# Yearly fractional declines that forecast background mortality, one for each
# age group: 0-4, 5-14, 15-24, ..., 85-94 and 95 and over.
mortality_declines <- c(0.0159, 0.0159, 0.0090, 0.0090, 0.0090, 0.0107, 0.0107,
                        0.0083, 0.0083, 0.0069, 0.0069)

# Yearly fractional decline that forecasts the TB prevalence trend among immigrants.
immigrant_tb_decline <- 0.015

#' Input series for a location
#'
#' Builds the location's model inputs as complete series from 1950 to 2100.
#' Years after an input's last observed year are forecast with fixed rules:
#' background mortality declines by \code{mortality_declines}, the immigrant TB
#' prevalence trend declines by \code{immigrant_tb_decline}, treatment default
#' and quality hold their last values, and the aging denominators hold the value
#' at the end of their smoothed curve.
#'
#' @param loc location code, as given to \code{model_load}
#' @param inputs the location's model inputs, as \code{model_load} reads them
#' @return list of
#'   \describe{
#'     \item{loc}{the location code}
#'     \item{observed_through}{the last observed year of the mortality,
#'       immigration and immigrant TB prevalence inputs}
#'     \item{init_pop}{the initial population}
#'     \item{births}{monthly births, smoothed from the yearly input}
#'     \item{mortality}{yearly background mortality rates, one row per year
#'       and one column per age group}
#'     \item{mortality_monthly}{monthly background mortality rates, smoothed
#'       from \code{mortality}}
#'     \item{immigration}{yearly immigration volume, in millions, named by year}
#'     \item{immigration_age_dist}{age distribution of immigrants, one column
#'       per year}
#'     \item{immigrant_tb_trend}{yearly trend in TB prevalence among
#'       immigrants, named by year}
#'     \item{immigrant_tb_rr_age}{relative risk of active TB among immigrants
#'       by age group}
#'     \item{tx_default, tx_quality}{yearly treatment default and treatment
#'       quality from 1993, named by year}
#'     \item{aging_denom}{monthly aging denominators}
#'     \item{net_migration}{net migration rates by age group and nativity}
#'   }
#' @noRd
build_input_series <- function(loc, inputs) {
  immig <- inputs[["ImmigInputs"]]
  mortality <- location_mortality(loc, inputs[["BgMort"]])
  immigration <- full_length(annual_values(immig[["TotByYear"]]), "immigration")
  burden <- annual_values(immig[["TBBurdenImmig"]])
  list(
    loc = loc,
    observed_through = c(
      mortality = attr(mortality, "observed_through"),
      immigration = observed_end(immigration, legacy_observed_through[["immigration"]]),
      immigrant_tb_trend = last_year(burden)
    ),
    init_pop = inputs[["InitPop"]],
    births = SmoCurve(full_length(annual_values(inputs[["Births"]]), "births")),
    mortality = mortality,
    mortality_monthly = sapply(seq_len(ncol(mortality)), function(i) SmoCurve(mortality[, i])),
    immigration = immigration,
    immigration_age_dist = immigration_age_dist(immig[["AgeDist"]]),
    immigrant_tb_trend = immigrant_tb_trend(burden, loc),
    immigrant_tb_rr_age = immig[["RR_Active_TB_Age"]],
    tx_default = held_to_end(inputs[["TxInputs"]][[1]]),
    tx_quality = held_to_end(inputs[["TxInputs"]][[2]]),
    aging_denom = age_denom(loc, year_month_to_idx(input_series_end_year, 1)),
    net_migration = inputs[["NetMigrState"]]
  )
}

#' Yearly values named by year
#'
#' An input without year names is in the legacy layout, one value per year from
#' 1950.
#' @noRd
annual_values <- function(x) {
  years <- if (is.null(names(x))) {
    seq(input_series_start_year, length.out = length(x))
  } else {
    as.integer(names(x))
  }
  if (years[1] != input_series_start_year || any(diff(years) != 1)) {
    stop("Input series must run year by year from ", input_series_start_year, call. = FALSE)
  }
  stats::setNames(as.numeric(x), years)
}

#' Last year of yearly values named by year
#' @noRd
last_year <- function(x) as.integer(names(x)[length(x)])

#' Last observed year of yearly values: their last year, or \code{legacy} when
#' they run through \code{input_series_end_year}
#' @noRd
observed_end <- function(x, legacy) {
  if (last_year(x) >= input_series_end_year) legacy else last_year(x)
}

#' Yearly values through \code{input_series_end_year}, for an input MITUS does
#' not yet forecast
#' @noRd
full_length <- function(x, series) {
  if (last_year(x) < input_series_end_year) {
    stop("The ", series, " input ends in ", last_year(x), "; MITUS has no forecast rule for ",
         series, " yet", call. = FALSE)
  }
  x[as.character(input_series_start_year:input_series_end_year)]
}

#' Age distribution of immigrants with one column per year, 1950-2100
#'
#' Columns without year names are in the legacy layout, one per year from 1950.
#' @noRd
immigration_age_dist <- function(age_dist) {
  years <- if (is.null(colnames(age_dist))) {
    seq(input_series_start_year, length.out = ncol(age_dist))
  } else {
    as.integer(colnames(age_dist))
  }
  if (max(years) < input_series_end_year) {
    stop("The immigrant age distribution input ends in ", max(years),
         "; MITUS has no forecast rule for it yet", call. = FALSE)
  }
  age_dist[, match(input_series_start_year:input_series_end_year, years)]
}

#' Yearly background mortality rates for a location, 1950-2100
#'
#' The input's rates, with the years the location's life tables cover replaced
#' by mortality weighted to the location, and the years after the last observed
#' year forecast by \code{mortality_declines}.
#'
#' @param mortality the inputs' mortality table: a year column, then one column
#'   per age group
#' @return rates with one row per year and one column per age group, and the
#'   last observed year as attribute \code{observed_through}
#' @noRd
location_mortality <- function(loc, mortality) {
  mortality <- as.matrix(mortality)
  years <- as.integer(mortality[, 1])
  rates <- matrix(NA_real_, length(input_series_start_year:input_series_end_year),
                  ncol(mortality) - 1)
  rates[match(years, input_series_start_year:input_series_end_year), ] <- mortality[, -1]
  weighted <- weight_mort(loc)
  rates[year_to_idx(as.integer(rownames(weighted))), ] <- weighted
  observed <- if (max(years) >= input_series_end_year) {
    legacy_observed_through[["mortality"]]
  } else {
    max(years)
  }
  rates <- forecast_mortality(rates, observed + 1)
  attr(rates, "observed_through") <- observed
  rates
}

#' Forecast yearly mortality by the fixed \code{mortality_declines}
#'
#' @param mortality mortality rates, one row per year from 1950 and one column
#'   per age group
#' @param from first forecast year
#' @return \code{mortality} with the years from \code{from} replaced
#' @noRd
forecast_mortality <- function(mortality, from) {
  for (j in year_to_idx(from):nrow(mortality)) {
    mortality[j, ] <- mortality[j - 1, ] * (1 - mortality_declines)
  }
  mortality
}

#' Trend in TB prevalence among immigrants, 1950-2100
#'
#' The yearly burden as a crude rate, declined by \code{immigrant_tb_decline}
#' each year after its data; for ND, doubled over 2010-2020.
#'
#' @param burden yearly TB burden among immigrants, named by year
#' @return yearly trend, named by year
#' @noRd
immigrant_tb_trend <- function(burden, loc) {
  past <- burden * (90 / 1e5)
  later <- seq(last_year(burden) + 1, input_series_end_year)
  future <- Reduce(function(prev, year) prev * (1 - immigrant_tb_decline), later,
                   accumulate = TRUE, init = past[[length(past)]])[-1]
  trend <- stats::setNames(c(unname(past), future), c(names(past), later))
  if (loc == "ND") {
    trend[as.character(2010:2020)] <- seq(from = trend[["2010"]], to = trend[["2010"]] * 2,
                                          length.out = 11)
  }
  trend
}

#' Yearly values of a table with a year column and a value column, held at
#' their last value through \code{input_series_end_year}, named by year
#' @noRd
held_to_end <- function(table) {
  years <- table[, 1]
  values <- table[, 2]
  later <- seq(max(years) + 1, input_series_end_year)
  stats::setNames(c(values, rep(values[length(values)], length(later))), c(years, later))
}

#' Whether the location's immigration input is forecast from data ending before
#' 2020, the inputs the 2020 pandemic immigration reduction corrects
#' @noRd
immigration_forecast_in_2020 <- function(series) {
  series$observed_through[["immigration"]] < 2020
}

#' Input series of the loaded location
#'
#' @return the list from \code{build_input_series}, as built by
#'   \code{model_load}
#' @noRd
input_series <- function() {
  if (!exists("input_series", envir = .mitus_state, inherits = FALSE)) {
    stop("Run model_load() before building model parameters", call. = FALSE)
  }
  get("input_series", envir = .mitus_state, inherits = FALSE)
}
