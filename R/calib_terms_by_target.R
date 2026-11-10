# Log-likelihood terms comparing one location's simulation with its calibration
# targets, one function per target. Each term takes the target's rows of the fit
# table from calib_fit(): the observed values and the simulated values they are
# compared with. Each term is a log-likelihood less its value at a reference
# fit: the targets themselves, or for Dirichlet-multinomial terms the targets
# plus a small count. A term is near zero for a good fit. The beta and
# Dirichlet-multinomial densities peak slightly away from the reference, so
# those terms can be slightly positive.
#
# wts are the importance weights by year set by model_load().

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

#' One column of target or fit rows as a matrix with one row per year or period
#' and one column per group
#'
#' @param d target or fit rows
#' @param groups the groups, in column order
#' @param column the column of \code{d} holding the values
#' @return list of the periods and the matrix of values
#' @noRd
obs_by_group <- function(d, groups, column = "value") {
  periods <- target_periods(d)
  values <- matrix(NA_real_, nrow(periods), length(groups))
  values[cbind(match(d$start, periods$start), match(d$group, groups))] <- d[[column]]
  list(periods = periods, values = values)
}

#' Rows of fit rows within the importance weights
#' @noRd
weighted_rows <- function(f, wts) f[weighted_years(f$start, wts), , drop = FALSE]

#' Total cases each year
#'
#' Normal, centred on the reported cases with a 95% interval of +/- 10%. Years
#' with no reported cases are left out. The weight of the latest weighted year
#' is doubled; when that year has no reported cases, no year is doubled.
#'
#' @param f fit rows of the target
#' @noRd
llik_cases_total <- function(f, wts) {
  f <- weighted_rows(f, wts)
  w <- year_wts(wts, f$start)
  w[length(w)] <- 2 * w[length(w)]
  reported <- !is.na(f$observed) & f$observed != 0
  obs <- f$observed[reported]
  llik_normal(obs, f$simulated[reported], sd_from_interval(obs, 0.1), w[reported])
}

#' Decline in total cases, 1953-1993
#'
#' Normal around the expected cases, with a 95% interval of +/- 20%.
#'
#' @param f fit rows of the target
#' @noRd
llik_cases_decline <- function(f, wts) {
  llik_normal(f$observed, f$simulated, sd_from_interval(f$observed, 0.2), year_wts(wts, f$start))
}

#' Cases in one nativity group, in 5-year buckets
#'
#' Normal, with a 95% interval of +/- 5%.
#'
#' @param f fit rows of the target for one nativity group
#' @noRd
llik_cases_nativity_5yr <- function(f, wts) {
  llik_normal(f$observed, f$simulated, sd_from_interval(f$observed, 0.05),
              period_wts(wts, f$start, f$end))
}

#' Trend in cases by nativity across the two most recent 5-year buckets
#'
#' Normal on the change in log cases between the two buckets, with standard
#' deviation 0.005.
#'
#' @param f_usb,f_nusb fit rows of cases by nativity for each nativity group
#' @param usb_trend target for the US-born trend in place of the data's, or NA
#' @noRd
llik_cases_nativity_trend <- function(f_usb, f_nusb, usb_trend = NA) {
  trend <- function(f) {
    f <- utils::tail(f[order(f$start), , drop = FALSE], 2)
    c(obs = diff(log(f$observed)), sim = diff(log(f$simulated)))
  }
  nusb <- trend(f_nusb)
  usb <- trend(f_usb)
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
#' @param f fit rows of cases by age for one nativity group
#' @param f_totals fit rows of cases by nativity for the same group
#' @param rho correlation parameter
#' @noRd
llik_cases_age_5yr <- function(f, f_totals, wts, rho = 0.1) {
  observed <- obs_by_group(f, data_age_groups, "observed")
  periods <- observed$periods
  obs <- observed$values
  sim <- obs_by_group(f, data_age_groups, "simulated")$values
  totals <- f_totals$observed[match(periods$start, f_totals$start)]
  w <- period_wts(wts, periods$start, periods$end)
  terms <- vapply(seq_len(nrow(obs)), function(i) {
    o <- obs[i, ]
    s <- sim[i, ]
    suppressed <- is.na(o)
    if (any(suppressed)) {
      o <- c(o[!suppressed], totals[i] - sum(o[!suppressed]))
      s <- c(s[!suppressed], sum(s) - sum(s[!suppressed]))
    }
    llik_dirmult(o, s, rho, w[i], pad = 0.01)
  }, numeric(1))
  sum(terms)
}

#' Cases among people experiencing homelessness, in 5-year buckets
#'
#' Normal around the smoothed counts, with a 95% interval of +/- 10%.
#'
#' @param f fit rows of the target
#' @noRd
llik_cases_homeless <- function(f, wts) {
  llik_normal(f$observed, f$simulated, sd_from_interval(f$observed, 0.1),
              period_wts(wts, f$start, f$end))
}

#' Fraction of cases from recent transmission
#'
#' Beta with effective sample size 1000.
#'
#' @param f fit row of the target
#' @noRd
llik_recent_transmission <- function(f) {
  obs <- f$observed
  sim <- f$simulated
  stats::dbeta(obs, sim * 1000, (1 - sim) * 1000, log = TRUE) -
    stats::dbeta(obs, obs * 1000, (1 - obs) * 1000, log = TRUE)
}

#' Non-US-born cases among recent arrivals, in 5-year buckets
#'
#' Normal, with a 95% interval of +/- 10%.
#'
#' @param f fit rows of the target
#' @noRd
llik_cases_recent_entry <- function(f, wts) {
  llik_normal(f$observed, f$simulated, sd_from_interval(f$observed, 0.1),
              period_wts(wts, f$start, f$end))
}

#' Treatment outcomes each year: completed, discontinued and died
#'
#' Dirichlet-multinomial.
#'
#' @param f fit rows of the target
#' @param rho correlation parameter
#' @noRd
llik_tx_outcomes <- function(f, wts, rho = 0.01) {
  f <- weighted_rows(f, wts)
  observed <- obs_by_group(f, tx_outcome_groups, "observed")
  sim <- obs_by_group(f, tx_outcome_groups, "simulated")$values
  years <- observed$periods$start
  llik_dirmult(observed$values, sim, rho, year_wts(wts, years), pad = 0.01)
}

#' LTBI treatment starts
#'
#' Normal, centred on the estimate, with standard deviation (upper - lower) /
#' 1.96 from its interval.
#'
#' @param f fit row of the target
#' @noRd
llik_ltbi_tx_total <- function(f, targets) {
  est <- target_estimate(targets, "ltbi_tx_total")
  llik_normal(f$observed, f$simulated, (est$upper - est$lower) / 1.96, 1)
}

#' Fractions of LTBI treatment starts that are non-US-born and homeless
#'
#' Beta with effective sample size 100.
#'
#' @param f fit rows of the target
#' @noRd
llik_ltbi_tx_dist <- function(f) {
  f <- f[match(c("nusb", "homeless"), f$group), , drop = FALSE]
  obs <- f$observed
  sim <- f$simulated
  sum(stats::dbeta(obs, sim * 100, (1 - sim) * 100, log = TRUE)) -
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

#' Expected IGRA-positive fraction in each of the survey's age groups
#'
#' The survey's eight age groups run from 5-14 to 75+; the model's 75-84, 85-94
#' and 95+ are pooled into 75+.
#'
#' @param tested expected positive and negative results, one row per model age
#'   group
#' @noRd
igra_positive_fraction <- function(tested) {
  tested[9, ] <- colSums(tested[9:11, ])
  survey_ages <- 2:9
  tested[survey_ages, 1] / rowSums(tested[survey_ages, ])
}

#' IGRA-positive fraction by age in the survey year
#'
#' Beta on the survey's positive and negative counts in each age group. Counts
#' four times.
#'
#' @param f fit rows of the target for one nativity group, with the positive
#'   count as \code{observed}, the number tested as \code{total}, and the
#'   expected positive count as \code{simulated}
#' @noRd
llik_ltbi_prev <- function(f) {
  positive <- f$observed
  negative <- f$total - f$observed
  sim <- f$simulated / f$total
  obs <- positive / f$total
  (sum(stats::dbeta(sim, positive, negative, log = TRUE)) -
     sum(stats::dbeta(obs, positive, negative, log = TRUE))) * 4
}

#' Deaths with TB each year
#'
#' Normal, with a 95% interval of +/- 10%. Years without data are left out.
#'
#' @param f fit rows of the target
#' @noRd
llik_tb_deaths_total <- function(f, wts) {
  f <- weighted_rows(f, wts)
  f <- f[!is.na(f$observed), , drop = FALSE]
  llik_normal(f$observed, f$simulated, sd_from_interval(f$observed, 0.1), year_wts(wts, f$start))
}

#' Age distribution of deaths with TB each year
#'
#' Dirichlet-multinomial.
#'
#' @param f fit rows of the target
#' @param rho correlation parameter
#' @noRd
llik_tb_deaths_age <- function(f, wts, rho = 0.005) {
  f <- weighted_rows(f, wts)
  observed <- obs_by_group(f, data_age_groups, "observed")
  sim <- obs_by_group(f, data_age_groups, "simulated")$values
  years <- observed$periods$start
  llik_dirmult(observed$values, sim, rho, year_wts(wts, years), pad = 0.005)
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
#' @param f fit rows of the target
#' @noRd
llik_borgdorff <- function(f) {
  ss <- 854.5463 / 4
  obs <- f$observed
  sim <- f$simulated
  adj <- sum(obs * ss * log(obs) + (1 - obs) * ss * log(1 - obs))
  study_llik(sum(obs * ss * log(sim) + (1 - obs) * ss * log(1 - sim)) - adj)
}

#' Annual TB incidence after infection in a published cohort
#'
#' Binomial on each year's cases among those still TB-free, divided by 4.
#'
#' @param f fit rows of the target, with each year's cases as \code{observed},
#'   the people at risk as \code{total}, and the simulated rate as
#'   \code{simulated}
#' @noRd
llik_cohort <- function(f) {
  cases <- f$observed
  at_risk <- f$total
  rate <- f$simulated
  # Years with no cases have no reference rate, so count 0.01 cases there
  fit_cases <- replace(cases, cases == 0, 0.01)
  adj <- sum(
    fit_cases * log(fit_cases / at_risk) + (at_risk - fit_cases) * log(1 - fit_cases / at_risk)
  )
  study_llik(sum(cases * log(rate) + (at_risk - cases) * log(1 - rate)) / 4 - adj / 4)
}

#' Tiemersma's duration of untreated TB and its case fatality
#'
#' Normal on a duration of 3 years, with a 95% interval of +/- 0.5 years, and
#' beta on a case fatality of 0.9, with effective sample size 50.
#'
#' @param f fit rows of the target: the duration and the case fatality
#' @noRd
llik_tiemersma <- function(f) {
  value <- function(column, quantity) f[[column]][f$group == quantity]
  duration <- value("observed", "duration")
  fatality <- value("observed", "fatality")
  adj <- stats::dnorm(duration, duration, 0.5 / 1.96, log = TRUE) +
    stats::dbeta(fatality, fatality * 50, (1 - fatality) * 50, log = TRUE)
  sim_fatality <- value("simulated", "fatality")
  stats::dnorm(duration, value("simulated", "duration"), 0.5 / 1.96, log = TRUE) +
    stats::dbeta(fatality, sim_fatality * 50, (1 - sim_fatality) * 50, log = TRUE) - adj
}

#' Total population of one nativity group: each census after the first, and
#' the latest weighted year of the annual population
#'
#' Normal, with a 95% interval of +/- 5% of the last census.
#'
#' @param f_census fit rows of the census populations for the group
#' @param f_total fit rows of the group's annual total population
#' @noRd
llik_pop_total <- function(f_census, f_total, wts) {
  census <- f_census[order(f_census$start), , drop = FALSE][-1, , drop = FALSE]
  annual <- weighted_rows(f_total, wts)
  latest <- annual[annual$start == max(annual$start), , drop = FALSE]
  years <- c(census$start, latest$start)
  sd <- sd_from_interval(census$observed[nrow(census)], 0.05)
  llik_normal(c(census$observed, latest$observed), c(census$simulated, latest$simulated), sd,
              year_wts(wts, years))
}

#' Age distribution of the population by nativity, over the latest years
#'
#' Multinomial with effective sample size \code{ess}, on the population summed
#' over the latest \code{n_years} weighted years.
#'
#' @param f_usb,f_nusb fit rows of the annual population by age for each
#'   nativity group
#' @param ess effective sample size
#' @param n_years number of latest years pooled
#' @noRd
llik_pop_age <- function(f_usb, f_nusb, wts, ess = 500, n_years = 3) {
  years <- utils::tail(sort(unique(weighted_rows(f_usb, wts)$start)), n_years)
  pooled <- function(f, column) {
    f <- f[f$start %in% years, , drop = FALSE]
    shares <- vapply(data_age_groups, function(age) sum(f[[column]][f$group == age]), numeric(1))
    shares / sum(shares)
  }
  obs <- cbind(pooled(f_usb, "observed"), pooled(f_nusb, "observed"))
  sim <- cbind(pooled(f_usb, "simulated"), pooled(f_nusb, "simulated"))
  adj <- (sum(log(obs[, 1]) * obs[, 1]) + sum(log(obs[, 2]) * obs[, 2])) * ess
  (sum(log(sim[, 1]) * obs[, 1]) + sum(log(sim[, 2]) * obs[, 2])) * ess - adj
}

#' Population experiencing homelessness
#'
#' Normal, with standard deviation (upper - lower) / 2 / 1.96 from the
#' estimate's interval.
#'
#' @param f fit row of the target
#' @noRd
llik_homeless_pop <- function(f, targets) {
  est <- target_estimate(targets, "homeless_pop")
  llik_normal(f$observed, f$simulated, (est$upper - est$lower) / 2 / 1.96, 1)
}

#' All-cause deaths in the latest weighted year
#'
#' Normal, with a 95% interval of +/- 10%.
#'
#' @param f fit rows of the target
#' @noRd
llik_deaths_total <- function(f, wts) {
  f <- weighted_rows(f, wts)
  f <- f[f$start == max(f$start), , drop = FALSE]
  llik_normal(f$observed, f$simulated, sd_from_interval(f$observed, 0.1), year_wts(wts, f$start))
}

#' Age distribution of all-cause deaths in the latest two weighted years
#'
#' Dirichlet-multinomial.
#'
#' @param f fit rows of the target
#' @param rho correlation parameter
#' @noRd
llik_deaths_age <- function(f, wts, rho = 0.01) {
  f <- weighted_rows(f, wts)
  f <- f[f$start %in% utils::tail(sort(unique(f$start)), 2), , drop = FALSE]
  observed <- obs_by_group(f, data_age_groups, "observed")
  sim <- obs_by_group(f, data_age_groups, "simulated")$values
  years <- observed$periods$start
  llik_dirmult(observed$values, sim, rho, year_wts(wts, years), pad = 0.1)
}

#' Distribution of the population across mortality risk groups
#'
#' Dirichlet-multinomial in each model age group, against the distribution the
#' model assumes; there is no location data for this target.
#'
#' @param f fit rows of the target
#' @param rho correlation parameter
#' @noRd
llik_mort_dist <- function(f, wts, rho = 0.1) {
  terms <- vapply(model_age_groups, function(age) {
    groups <- paste0(age, "_M", 1:4)
    rows <- f[f$group %in% groups, , drop = FALSE]
    observed <- obs_by_group(rows, groups, "observed")
    sim <- obs_by_group(rows, groups, "simulated")$values
    years <- observed$periods$start
    llik_dirmult(observed$values, sim, rho, year_wts(wts, years), pad = 0.01)
  }, numeric(1))
  sum(terms)
}
