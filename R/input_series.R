# This file builds a location's input series: yearly or monthly values, from
# 1950 to 2100, for each model input that does not depend on the model
# parameters.
#
# build_input_series() takes a location's input object and returns the series.
# The forecast rules, set out below, are the only part of this file that
# produces values not derived from an input.
#
# Rules for inputs:
# - Yearly values are labelled by year, through names or a year column.
#   Unlabelled yearly values are taken as one per year from 1950.
# - An input that reaches 2100 is used as given. An input that stops earlier is
#   used through its last year, and every later year is forecast by the rule
#   for that input.
# - Background mortality takes its rates from the location's life table for
#   the years the life table covers, in place of the input's rates for those
#   years.
# - The package's own state and national model input files hold PPML's
#   forecasts after their data. For these files only, mortality and
#   immigration are used through the years in legacy_input_through, and later
#   years are forecast by the rules.

# First and last year of every input series.
input_series_start_year <- 1950
input_series_end_year <- 2100

# Last year used from the package's own state and national model input files.
legacy_input_through <- list(
  state = c(mortality = 2016, immigration = 2019),
  US = c(mortality = 2016, immigration = 2017)
)

# ---- Forecast rules ---------------------------------------------------------
# An input that stops before 2100 is extended from its last year through 2100
# by one of these rules. They are modelling assumptions, not data.
#
#   background mortality            declines each year by mortality_declines
#   immigration volume              grows each year by immigration_growth
#   immigrant age distribution      holds the last year's distribution
#   births                          hold the last value
#   immigrant TB prevalence trend   declines each year by immigrant_tb_decline
#   treatment default and quality   hold their last values
#   aging denominators              hold the value at the end of their smoothed curve

# Yearly fractional declines in background mortality, one for each age group:
# 0-4, 5-14, 15-24, ..., 85-94 and 95 and over.
mortality_declines <- c(0.0159, 0.0159, 0.0090, 0.0090, 0.0090, 0.0107, 0.0107,
                        0.0083, 0.0083, 0.0069, 0.0069)

# Yearly fractional growth in immigration volume.
immigration_growth <- 0.005

# Yearly fractional decline in the TB prevalence trend among immigrants.
immigrant_tb_decline <- 0.015
# ------------------------------------------------------------------------------

#' Build a location's input series
#'
#' Reads each input by its years and extends it through 2100 by the rules at
#' the top of this file. Mortality and the aging denominators also use the
#' location's life table and census population, through \code{weight_mort} and
#' \code{age_denom}.
#'
#' @param loc location code; selects the life table, the census population and
#'   the adjustments made only for ND
#' @param inputs the location's input object: a list with \code{BgMort},
#'   \code{InitPop}, \code{Births}, \code{TxInputs}, \code{NetMigrState} and
#'   \code{ImmigInputs}
#' @param legacy whether \code{inputs} is one of the package's own state and
#'   national model input files, used through \code{legacy_input_through}
#' @return list of
#'   \describe{
#'     \item{loc}{the location code}
#'     \item{input_through}{the last year taken from the mortality, immigration
#'       and immigrant TB prevalence inputs; later years are forecast}
#'     \item{init_pop}{initial population, as given}
#'     \item{births}{monthly births, smoothed from yearly births through 2100}
#'     \item{mortality}{yearly background mortality rates, one row per year and
#'       one column per age group}
#'     \item{mortality_monthly}{monthly background mortality rates, smoothed
#'       from \code{mortality}}
#'     \item{immigration}{yearly immigration volume in millions, named by year}
#'     \item{immigration_age_dist}{age distribution of immigrants, one column
#'       per year}
#'     \item{immigrant_tb_trend}{yearly trend in TB prevalence among
#'       immigrants, named by year}
#'     \item{immigrant_tb_rr_age}{relative risk of active TB among immigrants by
#'       age group, as given}
#'     \item{tx_default, tx_quality}{yearly treatment default and treatment
#'       quality from 1993, named by year}
#'     \item{aging_denom}{monthly aging denominators}
#'     \item{net_migration}{net migration rates by age group and nativity, as
#'       given}
#'   }
#' @noRd
build_input_series <- function(loc, inputs, legacy = FALSE) {
  through <- if (legacy) legacy_input_through[[if (loc == "US") "US" else "state"]]
  immig <- inputs[["ImmigInputs"]]
  mortality <- location_mortality(loc, inputs[["BgMort"]], through[["mortality"]])
  immigration <- through_year(annual_values(immig[["TotByYear"]]), through[["immigration"]])
  births <- annual_values(inputs[["Births"]])
  burden <- annual_values(immig[["TBBurdenImmig"]])
  list(
    loc = loc,
    input_through = c(
      mortality = attr(mortality, "input_through"),
      immigration = last_year(immigration),
      immigrant_tb_trend = last_year(burden)
    ),
    init_pop = inputs[["InitPop"]],
    births = SmoCurve(forecast_births(births)),
    mortality = mortality,
    mortality_monthly = sapply(seq_len(ncol(mortality)), function(i) SmoCurve(mortality[, i])),
    immigration = forecast_immigration(immigration),
    immigration_age_dist = immigration_age_dist(immig[["AgeDist"]]),
    immigrant_tb_trend = immigrant_tb_trend(burden, loc),
    immigrant_tb_rr_age = immig[["RR_Active_TB_Age"]],
    tx_default = held_to_end(inputs[["TxInputs"]][[1]]),
    tx_quality = held_to_end(inputs[["TxInputs"]][[2]]),
    aging_denom = age_denom(loc, year_month_to_idx(input_series_end_year, 1)),
    net_migration = inputs[["NetMigrState"]]
  )
}

#' Label yearly values by year
#'
#' Takes a numeric vector, named by year or unnamed, and returns it as a
#' numeric vector named by year. Unnamed values are labelled one per year from
#' 1950. Stops unless the years start in 1950 and run without gaps.
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

#' Return the last year of yearly values named by year
#' @noRd
last_year <- function(x) as.integer(names(x)[length(x)])

#' Return yearly values named by year through \code{year}, or all of them when
#' \code{year} is NULL
#' @noRd
through_year <- function(x, year) {
  if (is.null(year)) x else x[as.integer(names(x)) <= year]
}

#' Extend yearly values through 2100
#'
#' Takes yearly values named by year and a rule, a function of the last value
#' and the number of years to add that returns the added values. Returns the
#' values for 1950 to 2100, named by year.
#' @noRd
extended <- function(x, rule) {
  later <- seq_len(input_series_end_year - last_year(x))
  if (length(later) > 0) {
    x <- c(x, stats::setNames(rule(x[[length(x)]], length(later)), last_year(x) + later))
  }
  x[as.character(input_series_start_year:input_series_end_year)]
}

#' Extend yearly immigration volume through 2100, growing it by
#' \code{immigration_growth} each year after its last year
#' @noRd
forecast_immigration <- function(immigration) {
  extended(immigration, function(last, n) {
    Reduce(function(prev, year) prev * (1 + immigration_growth), seq_len(n),
           accumulate = TRUE, init = last)[-1]
  })
}

#' Extend yearly births through 2100, holding them at their last value
#' @noRd
forecast_births <- function(births) {
  extended(births, function(last, n) rep(last, n))
}

#' Return the immigrant age distribution with one column per year from 1950 to
#' 2100
#'
#' Takes a matrix with one row per age group and one column per year, labelled
#' by year or unlabelled from 1950. Years after its last column take the last
#' column's distribution.
#' @noRd
immigration_age_dist <- function(age_dist) {
  years <- if (is.null(colnames(age_dist))) {
    seq(input_series_start_year, length.out = ncol(age_dist))
  } else {
    as.integer(colnames(age_dist))
  }
  columns <- match(input_series_start_year:input_series_end_year, years)
  columns[is.na(columns)] <- which.max(years)
  age_dist[, columns]
}

#' Build a location's yearly background mortality rates, 1950 to 2100
#'
#' Takes the input's mortality table and replaces the years covered by
#' \code{weight_mort(loc)} with its weighted rates. Uses the result through
#' \code{through}, or through its last year when \code{through} is NULL, and
#' forecasts any later years by \code{forecast_mortality}.
#'
#' @param loc location code
#' @param mortality the input's mortality table: a year column, then one column
#'   per age group
#' @param through last year to use, or NULL
#' @return matrix with one row per year from 1950 and one column per age group,
#'   with the last year used as attribute \code{input_through}
#' @noRd
location_mortality <- function(loc, mortality, through = NULL) {
  mortality <- as.matrix(mortality)
  all_years <- input_series_start_year:input_series_end_year
  rates <- matrix(NA_real_, length(all_years), ncol(mortality) - 1)
  rates[match(as.integer(mortality[, 1]), all_years), ] <- mortality[, -1]
  weighted <- weight_mort(loc)
  rates[year_to_idx(as.integer(rownames(weighted))), ] <- weighted
  last <- max(all_years[stats::complete.cases(rates)])
  if (!is.null(through)) last <- min(last, through)
  if (last < input_series_end_year) {
    rates <- forecast_mortality(rates, last + 1)
  }
  attr(rates, "input_through") <- last
  rates
}

#' Forecast yearly mortality by \code{mortality_declines}
#'
#' @param mortality matrix with one row per year from 1950 and one column per
#'   age group
#' @param from first year to forecast
#' @return \code{mortality} with each year from \code{from} set to the previous
#'   year's rates less \code{mortality_declines}
#' @noRd
forecast_mortality <- function(mortality, from) {
  for (j in year_to_idx(from):nrow(mortality)) {
    mortality[j, ] <- mortality[j - 1, ] * (1 - mortality_declines)
  }
  mortality
}

#' Build the trend in TB prevalence among immigrants, 1950 to 2100
#'
#' Converts the yearly TB burden among immigrants to a crude rate, then extends
#' it by declining \code{immigrant_tb_decline} each year after its last year.
#' For ND, sets 2010 to 2020 to a line rising from the 2010 value to double it.
#'
#' @param burden yearly TB burden among immigrants, named by year
#' @param loc location code
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

#' Extend a table of yearly values through 2100, holding its last value
#'
#' Takes a table with a year column and a value column, and returns its values
#' from its first year through 2100, named by year.
#' @noRd
held_to_end <- function(table) {
  years <- table[, 1]
  values <- table[, 2]
  later <- seq(max(years) + 1, input_series_end_year)
  stats::setNames(c(values, rep(values[length(values)], length(later))), c(years, later))
}

#' Return whether an input series' immigration is forecast in every year from
#' 2020, that is, taken from its input only through a year before 2020
#' @noRd
immigration_forecast_in_2020 <- function(series) {
  series$input_through[["immigration"]] < 2020
}

#' Return the input series stored in the package state
#'
#' Stops if no input series is stored.
#' @noRd
input_series <- function() {
  if (!exists("input_series", envir = .mitus_state, inherits = FALSE)) {
    stop("Run model_load() before building model parameters", call. = FALSE)
  }
  get("input_series", envir = .mitus_state, inherits = FALSE)
}
