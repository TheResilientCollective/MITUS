# Log-likelihood terms comparing one location's simulation with its calibration
# targets, one function per target. Each term is a log-likelihood less its
# value at a reference fit: the targets themselves, or for Dirichlet-multinomial
# terms the targets plus a small count. A term is near zero for a good fit. The
# beta and Dirichlet-multinomial densities peak slightly away from the
# reference, so those terms can be slightly positive.
#
# Simulated values come in as annual series with one row per simulated year,
# the first row being the model's base year. Simulated counts are in millions,
# as the model reports them. Targets come from calib_targets(): each term reads
# its rows, and the years each row covers, from them. wts are the importance
# weights by year set by model_load().

# The model's 11 age groups.
model_age_groups <- c(
  "0_4", "5_14", "15_24", "25_34", "35_44", "45_54", "55_64", "65_74", "75_84", "85_94", "95p"
)

# The calibration data's ten age groups: the model's, with 85-94 and 95+ pooled.
data_age_groups <- c(model_age_groups[1:9], "85p")

# IGRA sensitivity and specificity for US-born people, non-US-born children
# under 5, and non-US-born people 5 and over.
igra_sens <- c(usb = 0.780, nusb_young = 0.712, nusb = 0.789)
igra_spec <- c(usb = 0.979, nusb_young = 0.989, nusb = 0.985)

#' Weighted normal log-likelihood, less its value at \code{sim = obs}
#'
#' @param obs observed values
#' @param sim simulated values on the scale of \code{obs}
#' @param sd standard deviations
#' @param w weights
#' @noRd
llik_normal <- function(obs, sim, sd, w) {
  sum(stats::dnorm(obs, sim, sd, log = TRUE) * w) - sum(stats::dnorm(obs, obs, sd, log = TRUE) * w)
}

#' Standard deviation of a normal whose 95% interval is \code{value} +/-
#' \code{relative} times \code{value}
#' @noRd
sd_from_interval <- function(value, relative) value * relative / 1.96

#' Weighted Dirichlet-multinomial log-likelihood, less its value at
#' \code{sim = obs + pad}
#'
#' @param obs observed counts, one row per year or bucket
#' @param sim simulated values, one row per year or bucket; only their
#'   proportions within each row matter
#' @param rho correlation parameter
#' @param w weights, one per row
#' @param pad small count added to \code{obs} for the reference fit, so that
#'   zero counts give finite densities
#' @noRd
llik_dirmult <- function(obs, sim, rho, w, pad) {
  sum(dDirMult(M = sim, n = obs, Rho = rho) * w) -
    sum(dDirMult(M = obs + pad, n = obs, Rho = rho) * w)
}

#' Calendar years as integers, whether stored as numbers, text or factors
#' @noRd
as_year <- function(years) as.integer(as.character(years))

#' Weights for the given years
#' @noRd
year_wts <- function(wts, years) wts[as.character(years)]

#' Which years fall within the importance weights
#' @noRd
weighted_years <- function(years, wts) as_year(years) %in% as.integer(names(wts))

#' Weights of multi-year periods: the mean of the yearly weights over each
#' period's years
#'
#' @param start,end first and last year of each period
#' @noRd
period_wts <- function(wts, start, end) {
  vapply(seq_along(start), function(i) mean(year_wts(wts, start[i]:end[i])), numeric(1))
}

#' Sum simulated values over periods of years
#'
#' @param sim annual simulated values, a vector or a matrix with one row per year
#' @param start,end first and last year of each period
#' @return one value per period, or a matrix with one row per period
#' @noRd
sum_periods <- function(sim, start, end) {
  rows <- lapply(seq_along(start), function(i) year_to_idx(start[i]:end[i]))
  if (is.null(dim(sim))) {
    return(vapply(rows, function(r) sum(sim[r]), numeric(1)))
  }
  t(vapply(rows, function(r) colSums(sim[r, , drop = FALSE]), numeric(ncol(sim))))
}

#' Combine the model's 85-94 and 95+ age groups into one 85+ column
#'
#' @param by_age matrix with one column per model age group
#' @return matrix with ten columns, the last for 85+
#' @noRd
merge_85_plus <- function(by_age) {
  merged <- by_age[, -11, drop = FALSE]
  merged[, 10] <- merged[, 10] + by_age[, 11]
  merged
}

#' Columns of a 5-year case table holding bucket counts, oldest first
#' @noRd
case_bucket_cols <- function(d) {
  grep("^X[0-9]{4}\\.[0-9]{4}$", colnames(d), value = TRUE)
}

#' Columns of the 5-year age table holding age groups, youngest first
#' @noRd
case_age_cols <- function(d) {
  grep("^X[0-9]", colnames(d), value = TRUE)
}

#' Observed values of target rows as a matrix with one row per year or period
#' and one column per group
#'
#' @param d target rows
#' @param groups the groups, in column order
#' @return list of the periods and the matrix of values
#' @noRd
obs_by_group <- function(d, groups) {
  periods <- target_periods(d)
  values <- matrix(NA_real_, nrow(periods), length(groups))
  values[cbind(match(d$start, periods$start), match(d$group, groups))] <- d$value
  list(periods = periods, values = values)
}

#' Total cases each year
#'
#' Normal, centred on the reported cases with a 95% interval of +/- 10%. Years
#' with no reported cases are left out. The weight of the latest weighted year
#' is doubled; when that year has no reported cases, no year is doubled.
#'
#' @param sim_cases simulated cases, in millions
#' @noRd
llik_cases_total <- function(sim_cases, targets, wts) {
  d <- target_obs(targets, "cases_total")
  d <- d[weighted_years(d$start, wts), , drop = FALSE]
  w <- year_wts(wts, d$start)
  w[length(w)] <- 2 * w[length(w)]
  sim <- sim_cases[year_to_idx(d$start)] * 1e6
  reported <- !is.na(d$value) & d$value != 0
  obs <- d$value[reported]
  llik_normal(obs, sim[reported], sd_from_interval(obs, 0.1), w[reported])
}

#' Decline in total cases, 1953-1993
#'
#' Normal around the expected cases, with a 95% interval of +/- 20%.
#'
#' @param sim_cases simulated cases, in millions
#' @noRd
llik_cases_decline <- function(sim_cases, targets, wts) {
  d <- target_obs(targets, "cases_decline")
  sim <- sim_cases[year_to_idx(d$start)] * 1e6
  llik_normal(d$value, sim, sd_from_interval(d$value, 0.2), year_wts(wts, d$start))
}

#' Cases in one nativity group, in 5-year buckets
#'
#' Normal, with a 95% interval of +/- 5%.
#'
#' @param sim_cases simulated cases in the nativity group, in millions
#' @param nativity \code{"usb"} or \code{"nusb"}
#' @noRd
llik_cases_nativity_5yr <- function(sim_cases, nativity, targets, wts) {
  d <- target_obs(targets, "cases_nativity", nativity)
  sim <- sum_periods(sim_cases, d$start, d$end) * 1e6
  w <- period_wts(wts, d$start, d$end)
  llik_normal(d$value, sim, sd_from_interval(d$value, 0.05), w)
}

#' Trend in cases by nativity across the two most recent 5-year buckets
#'
#' Normal on the change in log cases between the two buckets, with standard
#' deviation 0.005.
#'
#' @param sim_usb,sim_nusb simulated US-born and non-US-born cases, on any scale
#' @param usb_trend target for the US-born trend in place of the data's, or NA
#' @noRd
llik_cases_nativity_trend <- function(sim_usb, sim_nusb, targets, usb_trend = NA) {
  trend <- function(nativity, sim) {
    d <- utils::tail(target_obs(targets, "cases_nativity", nativity), 2)
    c(obs = diff(log(d$value)), sim = diff(log(sum_periods(sim, d$start, d$end))))
  }
  nusb <- trend("nusb", sim_nusb)
  usb <- trend("usb", sim_usb)
  obs <- c(nusb = nusb[["obs"]], usb = usb[["obs"]])
  if (!is.na(usb_trend)) {
    obs[["usb"]] <- usb_trend
  }
  llik_normal(obs, c(nusb[["sim"]], usb[["sim"]]), 0.005, 1)
}

#' Age distribution of cases in one nativity group, in 5-year buckets
#'
#' Dirichlet-multinomial. Age groups with suppressed counts are pooled into one
#' group holding the bucket's remaining cases, taken from the bucket's cases by
#' nativity.
#'
#' @param sim_by_age simulated cases in the nativity group, in millions, one
#'   column per model age group
#' @param nativity \code{"usb"} or \code{"nusb"}
#' @param rho correlation parameter
#' @noRd
llik_cases_age_5yr <- function(sim_by_age, nativity, targets, wts, rho = 0.1) {
  by_age <- obs_by_group(target_obs(targets, "cases_age", nativity), data_age_groups)
  periods <- by_age$periods
  obs <- by_age$values
  totals_d <- target_obs(targets, "cases_nativity", nativity)
  totals <- totals_d$value[match(periods$start, totals_d$start)]
  sim <- sum_periods(merge_85_plus(sim_by_age), periods$start, periods$end)
  w <- period_wts(wts, periods$start, periods$end)
  terms <- vapply(seq_len(nrow(obs)), function(i) {
    o <- obs[i, ]
    s <- sim[i, ]
    suppressed <- is.na(o)
    if (any(suppressed)) {
      o <- c(o[!suppressed], totals[i] - sum(o[!suppressed]))
      s <- c(s[!suppressed], sum(s) - sum(s[!suppressed]))
    }
    llik_dirmult(o, s * 1e6, rho, w[i], pad = 0.01)
  }, numeric(1))
  sum(terms)
}

#' Cases among people experiencing homelessness, in 5-year buckets
#'
#' Normal around the smoothed counts, with a 95% interval of +/- 10%.
#'
#' @param sim_cases simulated homeless cases, in millions
#' @noRd
llik_cases_homeless <- function(sim_cases, targets, wts) {
  d <- target_obs(targets, "cases_homeless")
  sim <- sum_periods(sim_cases, d$start, d$end) * 1e6
  w <- period_wts(wts, d$start, d$end)
  llik_normal(d$value, sim, sd_from_interval(d$value, 0.1), w)
}

#' Fraction of cases from recent transmission
#'
#' Beta with effective sample size 1000.
#'
#' @param sim_fraction simulated fraction of cases from recent transmission over
#'   the target's years
#' @noRd
llik_recent_transmission <- function(sim_fraction, targets) {
  obs <- target_estimate(targets, "recent_transmission")$estimate
  stats::dbeta(obs, sim_fraction * 1000, (1 - sim_fraction) * 1000, log = TRUE) -
    stats::dbeta(obs, obs * 1000, (1 - obs) * 1000, log = TRUE)
}

#' Non-US-born cases among recent arrivals, in 5-year buckets
#'
#' Normal, with a 95% interval of +/- 10%.
#'
#' @param sim_cases simulated cases among recent arrivals, in millions
#' @noRd
llik_cases_recent_entry <- function(sim_cases, targets, wts) {
  d <- target_obs(targets, "cases_recent_entry")
  sim <- sum_periods(sim_cases, d$start, d$end) * 1e6
  w <- period_wts(wts, d$start, d$end)
  llik_normal(d$value, sim, sd_from_interval(d$value, 0.1), w)
}

#' Treatment outcomes each year: completed, discontinued and died
#'
#' Dirichlet-multinomial.
#'
#' @param sim_outcomes simulated outcomes, three columns in that order, on any
#'   scale
#' @param rho correlation parameter
#' @noRd
llik_tx_outcomes <- function(sim_outcomes, targets, wts, rho = 0.01) {
  d <- target_obs(targets, "tx_outcomes")
  outcomes <- obs_by_group(d[weighted_years(d$start, wts), , drop = FALSE], tx_outcome_groups)
  years <- outcomes$periods$start
  sim <- sim_outcomes[year_to_idx(years), , drop = FALSE]
  llik_dirmult(outcomes$values, sim, rho, year_wts(wts, years), pad = 0.01)
}

#' LTBI treatment starts
#'
#' Normal, centred on the estimate, with standard deviation (upper - lower) /
#' 1.96 from its interval.
#'
#' @param sim_starts simulated treatment starts in the estimate's year, in
#'   millions
#' @noRd
llik_ltbi_tx_total <- function(sim_starts, targets) {
  est <- target_estimate(targets, "ltbi_tx_total")
  llik_normal(est$estimate, sim_starts * 1e6, (est$upper - est$lower) / 1.96, 1)
}

#' Fractions of LTBI treatment starts that are non-US-born and homeless
#'
#' Beta with effective sample size 100.
#'
#' @param sim_fractions simulated non-US-born and homeless fractions
#' @noRd
llik_ltbi_tx_dist <- function(sim_fractions, targets) {
  obs <- target_estimate(targets, "ltbi_tx_share", c("nusb", "homeless"))$estimate
  sum(stats::dbeta(obs, sim_fractions * 100, (1 - sim_fractions) * 100, log = TRUE)) -
    sum(stats::dbeta(obs, obs * 100, (1 - obs) * 100, log = TRUE))
}

#' Expected IGRA results for people with and without LTBI
#'
#' @param ltbi,no_ltbi people with and without LTBI, one value per age group
#' @param sens,spec IGRA sensitivity and specificity
#' @return matrix of expected positive and negative results, one row per age
#'   group
#' @noRd
igra_results <- function(ltbi, no_ltbi, sens, spec) {
  outer(ltbi, c(sens, 1 - sens)) + outer(no_ltbi, c(1 - spec, spec))
}

#' IGRA-positive fraction by age in the survey year
#'
#' Beta on the survey's positive and negative counts in each age group from
#' 5-14 to 75+, the survey's eight age groups. Counts four times.
#'
#' @param tested expected positive and negative results, one row per model age
#'   group
#' @param survey survey rows, one per survey age group, with the positive count
#'   as \code{value} and the number tested as \code{total}
#' @noRd
llik_ltbi_prev <- function(tested, survey) {
  # Pool the model's 75-84, 85-94 and 95+ into the survey's 75+
  tested[9, ] <- colSums(tested[9:11, ])
  survey_ages <- 2:9
  sim <- tested[survey_ages, 1] / rowSums(tested[survey_ages, ])
  positive <- survey$value
  negative <- survey$total - survey$value
  obs <- positive / survey$total
  (sum(stats::dbeta(sim, positive, negative, log = TRUE)) -
     sum(stats::dbeta(obs, positive, negative, log = TRUE))) * 4
}

#' IGRA-positive fraction by age among US-born people
#'
#' @param sim_ltbi,sim_pop simulated US-born people with LTBI and in total in the
#'   survey year, by model age group, on any common scale
#' @noRd
llik_ltbi_prev_usb <- function(sim_ltbi, sim_pop, targets) {
  tested <- igra_results(sim_ltbi, sim_pop - sim_ltbi, igra_sens[["usb"]], igra_spec[["usb"]])
  llik_ltbi_prev(tested, target_obs(targets, "ltbi_prev", "usb"))
}

#' IGRA-positive fraction by age among non-US-born people
#'
#' Children under 5 take their own IGRA accuracy. The survey starts at age 5, so
#' their results do not enter the likelihood.
#'
#' @param sim_ltbi,sim_pop simulated non-US-born people with LTBI and in total in
#'   the survey year, by model age group, on any common scale
#' @noRd
llik_ltbi_prev_nusb <- function(sim_ltbi, sim_pop, targets) {
  no_ltbi <- sim_pop - sim_ltbi
  tested <- rbind(
    igra_results(sim_ltbi[1], no_ltbi[1], igra_sens[["nusb_young"]], igra_spec[["nusb_young"]]),
    igra_results(sim_ltbi[-1], no_ltbi[-1], igra_sens[["nusb"]], igra_spec[["nusb"]])
  )
  llik_ltbi_prev(tested, target_obs(targets, "ltbi_prev", "nusb"))
}

#' Deaths with TB each year
#'
#' Normal, with a 95% interval of +/- 10%. Years without data are left out.
#'
#' @param sim_deaths simulated deaths with TB, in millions, one column per model
#'   age group
#' @noRd
llik_tb_deaths_total <- function(sim_deaths, targets, wts) {
  d <- target_obs(targets, "tb_deaths_total")
  d <- d[weighted_years(d$start, wts) & !is.na(d$value), , drop = FALSE]
  sim <- rowSums(sim_deaths[year_to_idx(d$start), , drop = FALSE]) * 1e6
  llik_normal(d$value, sim, sd_from_interval(d$value, 0.1), year_wts(wts, d$start))
}

#' Age distribution of deaths with TB each year
#'
#' Dirichlet-multinomial.
#'
#' @param sim_deaths simulated deaths with TB, one column per model age group, on
#'   any scale
#' @param rho correlation parameter
#' @noRd
llik_tb_deaths_age <- function(sim_deaths, targets, wts, rho = 0.005) {
  d <- target_obs(targets, "tb_deaths_age")
  by_age <- obs_by_group(d[weighted_years(d$start, wts), , drop = FALSE], data_age_groups)
  years <- by_age$periods$start
  sim <- merge_85_plus(sim_deaths[year_to_idx(years), , drop = FALSE])
  llik_dirmult(by_age$values, sim, rho, year_wts(wts, years), pad = 0.005)
}

#' Probability of not having developed TB at each time since infection
#'
#' Uses the progression parameters of the 15-24 age group, row 3 of
#' \code{Mpfast} and \code{Mrslow}, averaged over the four progression risk
#' groups in the proportions of the package constant \code{dist_gen}.
#'
#' @param progression list of \code{Mpfast}, \code{Mrslow}, \code{rfast} and
#'   \code{rRecov} from the model parameters, with rates per month
#' @param times years since infection
#' @noRd
tb_free_probability <- function(progression, times) {
  pfast <- progression$Mpfast[3, ]
  rslow <- progression$Mrslow[3, ] * 12
  rfast <- progression$rfast * 12
  rrecov <- progression$rRecov * 12
  p <- matrix(NA_real_, length(pfast), length(times))
  for (i in seq_along(pfast)) {
    p[i, ] <- pfast[i] * (1 - (1 - exp(-(rfast + rrecov) * times)) * (rfast / (rfast + rrecov))) +
      (1 - pfast[i]) *
        (1 - (1 - exp(-(rslow[i] + rrecov) * times)) * (rslow[i] / (rslow[i] + rrecov)))
  }
  as.numeric(t(p) %*% colSums(dist_gen))
}

#' Likelihood of a published study, or -10^4 when it errors or is NaN or -Inf
#'
#' @param value an expression giving the log-likelihood
#' @noRd
study_llik <- function(value) {
  llik <- tryCatch(value, error = function(e) -Inf)
  if (is.nan(llik) || llik == -Inf) -1e4 else llik
}

#' Borgdorff's cumulative TB incidence after infection
#'
#' Binomial on each interval's share of cases, with effective sample size
#' 854.5463 / 4.
#'
#' @param progression see \code{tb_free_probability}
#' @noRd
llik_borgdorff <- function(progression, targets) {
  ss <- 854.5463 / 4
  d <- targets$studies$borgdorff_data
  obs <- diff(-d[, "mean"])
  adj <- sum(obs * ss * log(obs) + (1 - obs) * ss * log(1 - obs))
  study_llik({
    p_free <- tb_free_probability(progression, d[, "Years"])
    sim <- diff(-(1 - (1 - p_free) / (1 - p_free)[nrow(d)]))
    sum(obs * ss * log(sim) + (1 - obs) * ss * log(1 - sim)) - adj
  })
}

#' Annual TB incidence after infection in a published cohort
#'
#' Binomial on each year's cases among those still TB-free, divided by 4.
#'
#' @param progression see \code{tb_free_probability}
#' @param d cohort data, one row per year: the year, people at risk
#'   (\code{At_Risk}), and cases in the third column
#' @noRd
llik_cohort <- function(progression, d) {
  at_risk <- d[, "At_Risk"]
  cases <- d[, 3]
  # Years with no cases have no reference rate, so count 0.01 cases there
  fit_cases <- replace(cases, cases == 0, 0.01)
  adj <- sum(
    fit_cases * log(fit_cases / at_risk) + (at_risk - fit_cases) * log(1 - fit_cases / at_risk)
  )
  n_years <- nrow(d)
  study_llik({
    p_free <- tb_free_probability(progression, 0:n_years)
    rate <- -log(1 - diff(-p_free) / p_free[-(n_years + 1)])
    sum(cases * log(rate) + (at_risk - cases) * log(1 - rate)) / 4 - adj / 4
  })
}

#' Ferebee's TB incidence after infection
#' @noRd
llik_ferebee <- function(progression, targets) {
  llik_cohort(progression, targets$studies$ferebee_data)
}

#' Sutherland's TB incidence after infection
#' @noRd
llik_sutherland <- function(progression, targets) {
  llik_cohort(progression, targets$studies$sutherland_data)
}

#' Tiemersma's duration of untreated TB and its case fatality
#'
#' Normal on a duration of 3 years, with a 95% interval of +/- 0.5 years, and
#' beta on a case fatality of 0.9, with effective sample size 50.
#'
#' @param self_cure,mortality monthly rates of self-cure and death from
#'   untreated TB
#' @noRd
llik_tiemersma <- function(self_cure, mortality) {
  self_cure <- self_cure * 12
  mortality <- mortality * 12
  duration <- 1 / (self_cure + mortality)
  fatality <- mortality / (self_cure + mortality)
  adj <- stats::dnorm(3, 3, 0.5 / 1.96, log = TRUE) +
    stats::dbeta(0.9, 0.9 * 50, (1 - 0.9) * 50, log = TRUE)
  stats::dnorm(3, duration, 0.5 / 1.96, log = TRUE) +
    stats::dbeta(0.9, fatality * 50, (1 - fatality) * 50, log = TRUE) - adj
}

#' Weighted years of the annual population rows, earliest first
#' @noRd
population_years <- function(targets, wts) {
  years <- unique(target_obs(targets, "population")$start)
  sort(years[weighted_years(years, wts)])
}

#' Total population of one nativity group: each census after the first, and
#' the latest weighted year of the annual population
#'
#' @param nativity \code{"usb"} or \code{"nusb"}
#' @return list of years, population, and the population of the last census,
#'   which sets the interval width
#' @noRd
pop_total_obs <- function(targets, nativity, wts) {
  census <- target_obs(targets, "population_census", nativity)[-1, , drop = FALSE]
  annual <- target_obs(targets, "population", nativity)
  latest <- max(population_years(targets, wts))
  list(
    years = c(census$start, latest),
    pop = c(census$value, sum(annual$value[annual$start == latest])),
    last_census_pop = census$value[nrow(census)]
  )
}

#' Total population of one nativity group
#'
#' Normal, with a 95% interval of +/- 5% of the last census.
#'
#' @param sim_pop simulated population of the group, in millions
#' @param nativity \code{"usb"} or \code{"nusb"}
#' @noRd
llik_pop_total <- function(sim_pop, nativity, targets, wts) {
  tg <- pop_total_obs(targets, nativity, wts)
  sim <- sim_pop[year_to_idx(tg$years)] * 1e6
  llik_normal(
    tg$pop, sim, sd_from_interval(tg$last_census_pop, 0.05), year_wts(wts, tg$years)
  )
}

#' Age distribution of the population by nativity, over the latest years
#'
#' Multinomial with effective sample size \code{ess}, on the population summed
#' over the latest \code{n_years} weighted years.
#'
#' @param sim_usb,sim_nusb simulated US-born and non-US-born population, one
#'   column per model age group, on any common scale
#' @param ess effective sample size
#' @param n_years number of latest years pooled
#' @noRd
llik_pop_age <- function(sim_usb, sim_nusb, targets, wts, ess = 500, n_years = 3) {
  years <- utils::tail(population_years(targets, wts), n_years)
  pooled <- function(nativity) {
    d <- target_obs(targets, "population", nativity)
    d <- d[d$start %in% years, , drop = FALSE]
    vapply(data_age_groups, function(age) sum(d$value[d$group == age]), numeric(1))
  }
  pop_usb <- pooled("usb")
  pop_nusb <- pooled("nusb")
  obs <- cbind(pop_usb / sum(pop_usb), pop_nusb / sum(pop_nusb))
  adj <- (sum(log(obs[, 1]) * obs[, 1]) + sum(log(obs[, 2]) * obs[, 2])) * ess
  by_age <- cbind(
    colSums(sim_usb[year_to_idx(years), , drop = FALSE]),
    colSums(sim_nusb[year_to_idx(years), , drop = FALSE])
  )
  by_age <- t(merge_85_plus(t(by_age)))
  sim <- cbind(by_age[, 1] / sum(by_age[, 1]), by_age[, 2] / sum(by_age[, 2]))
  (sum(log(sim[, 1]) * obs[, 1]) + sum(log(sim[, 2]) * obs[, 2])) * ess - adj
}

#' Population experiencing homelessness
#'
#' Normal, with standard deviation (upper - lower) / 2 / 1.96 from the
#' estimate's interval.
#'
#' @param sim_homeless simulated homeless population in the estimate's year, in
#'   millions, as is the target
#' @noRd
llik_homeless_pop <- function(sim_homeless, targets) {
  est <- target_estimate(targets, "homeless_pop")
  llik_normal(est$estimate, sim_homeless, (est$upper - est$lower) / 2 / 1.96, 1)
}

#' All-cause deaths in the latest weighted year
#'
#' Normal, with a 95% interval of +/- 10%.
#'
#' @param sim_deaths simulated deaths, in millions, one column per model age group
#' @noRd
llik_deaths_total <- function(sim_deaths, targets, wts) {
  d <- target_obs(targets, "deaths_total")
  year <- max(d$start[weighted_years(d$start, wts)])
  obs <- d$value[d$start == year]
  sim <- rowSums(sim_deaths[year_to_idx(year), , drop = FALSE]) * 1e6
  llik_normal(obs, sim, sd_from_interval(obs, 0.1), year_wts(wts, year))
}

#' Age distribution of all-cause deaths in the latest two weighted years
#'
#' Dirichlet-multinomial.
#'
#' @param sim_deaths simulated deaths, one column per model age group, on any
#'   scale
#' @param rho correlation parameter
#' @noRd
llik_deaths_age <- function(sim_deaths, targets, wts, rho = 0.01) {
  d <- target_obs(targets, "deaths_age")
  years <- utils::tail(sort(unique(d$start[weighted_years(d$start, wts)])), 2)
  by_age <- obs_by_group(d[d$start %in% years, , drop = FALSE], data_age_groups)
  sim <- merge_85_plus(sim_deaths[year_to_idx(years), , drop = FALSE])
  llik_dirmult(by_age$values, sim, rho, year_wts(wts, years), pad = 0.1)
}

#' Distribution of the population across mortality risk groups, 2015-2016
#'
#' Dirichlet-multinomial in each age group, against the distribution in the
#' package constant \code{dist_gen}, which the model assumes; there is no
#' location data for this target.
#'
#' @param sim_by_group simulated population in millions, four mortality risk
#'   group columns for each model age group in turn
#' @param rho correlation parameter
#' @noRd
llik_mort_dist <- function(sim_by_group, wts, rho = 0.1) {
  years <- 2015:2016
  obs <- matrix(rowSums(dist_gen), length(years), 4, byrow = TRUE)
  rows <- year_to_idx(years)
  terms <- vapply(seq_along(model_age_groups), function(age) {
    sim <- sim_by_group[rows, (1:4) + 4 * (age - 1)] * 1e6
    llik_dirmult(obs, sim, rho, year_wts(wts, years), pad = 0.01)
  }, numeric(1))
  sum(terms)
}
