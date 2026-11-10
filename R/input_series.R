# Input series for one location: the model inputs as complete series from 1950
# to 2100, built once per location by model_load. None of them depend on the
# model parameters; param_init applies the parameters to them.

# Last year of the input series.
input_series_end_year <- 2100

# Yearly fractional declines that forecast background mortality, applied from
# 2017 to the year column and the first ten age groups of the mortality table.
mortality_forecast_start <- 2017
mortality_declines <- c(0.0159, 0.0159, 0.0090, 0.0090, 0.0090, 0.0107, 0.0107,
                        0.0083, 0.0083, 0.0069, 0.0069)

#' Input series for a location
#'
#' Builds the location's model inputs as complete series from 1950 to 2100.
#' Years past an input's data are forecast with fixed rules: background
#' mortality declines at fixed yearly rates from 2017, the immigrant TB
#' prevalence trend declines 1.5\% a year after its data, treatment default and
#' quality hold their last values, and the aging denominators hold their
#' January 2020 values.
#'
#' @param loc location code, as given to \code{model_load}
#' @param inputs the location's model inputs, as \code{model_load} reads them
#' @return list of
#'   \describe{
#'     \item{loc}{the location code}
#'     \item{init_pop}{the initial population}
#'     \item{births}{monthly births, smoothed from the yearly input}
#'     \item{mortality}{yearly background mortality rates, one row per year
#'       and one column per age group}
#'     \item{mortality_monthly}{monthly background mortality rates, smoothed
#'       from \code{mortality}}
#'     \item{immigration}{yearly immigration volume, in millions}
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
  years <- year_to_idx(1950):year_to_idx(input_series_end_year)
  mortality <- forecast_mortality(weighted_mortality(loc, inputs))
  immig <- inputs[["ImmigInputs"]]
  list(
    loc = loc,
    init_pop = inputs[["InitPop"]],
    births = SmoCurve(inputs[["Births"]]),
    mortality = mortality,
    mortality_monthly = sapply(seq_len(ncol(mortality)), function(i) SmoCurve(mortality[, i])),
    immigration = immig[[1]][years],
    immigration_age_dist = immig[["AgeDist"]],
    immigrant_tb_trend = immigrant_tb_trend(inputs, loc),
    immigrant_tb_rr_age = immig[["RR_Active_TB_Age"]],
    tx_default = held_to_end(inputs[["TxInputs"]][[1]]),
    tx_quality = held_to_end(inputs[["TxInputs"]][[2]]),
    aging_denom = age_denom(loc, year_month_to_idx(input_series_end_year, 1)),
    net_migration = inputs[["NetMigrState"]]
  )
}

#' Yearly background mortality rates from the inputs, with the years the life
#' tables cover replaced by mortality weighted to the location
#'
#' @return the inputs' mortality table, a year column followed by one column
#'   per age group
#' @noRd
weighted_mortality <- function(loc, inputs) {
  mortality <- as.matrix(inputs[["BgMort"]])
  if (loc == "US") {
    mortality[1:68, 2:12] <- weight_mort(loc)
  } else {
    mortality[10:67, 2:12] <- weight_mort(loc)
  }
  mortality
}

#' Forecast yearly mortality from \code{mortality_forecast_start} by the fixed
#' \code{mortality_declines}
#'
#' @param mortality mortality table with a year column, one row per year from
#'   1950
#' @return the age-group columns of the table, with the forecast years replaced
#' @noRd
forecast_mortality <- function(mortality) {
  cols <- seq_along(mortality_declines)
  for (j in year_to_idx(mortality_forecast_start):nrow(mortality)) {
    mortality[j, cols] <- mortality[j - 1, cols] * (1 - mortality_declines)
  }
  mortality[, -1]
}

#' Trend in TB prevalence among immigrants: the inputs' yearly burden, declined
#' 1.5\% a year after its data; for ND, doubled over 2010-2020
#' @noRd
immigrant_tb_trend <- function(inputs, loc) {
  trend <- crude_rate(inputs, loc)
  if (loc == "ND") {
    trend[year_to_idx(2010):year_to_idx(2020)] <- seq(
      from = trend[year_to_idx(2010)], to = trend[year_to_idx(2010)] * 2, length.out = 11
    )
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
