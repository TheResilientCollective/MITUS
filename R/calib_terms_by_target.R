# Log-likelihood terms comparing one location's simulation with its calibration
# targets. Each term is zero when the simulation matches its targets exactly,
# and negative otherwise.
#
# Simulated values come in as annual series with one row per simulated year,
# the first row being the model's base year. Targets come from calib_targets(),
# and wts are the importance weights by year set by model_load().

# The model's 11 age groups. Calibration data groups the last two as 85+.
model_age_groups <- c(
  "0_4", "5_14", "15_24", "25_34", "35_44", "45_54", "55_64", "65_74", "75_84", "85_94", "95p"
)

# Simulated years summed into each 5-year bucket of case targets, and the year
# whose importance weight each bucket carries.
case_bucket_years <- list(1995:1999, 2000:2004, 2005:2009, 2010:2014, 2015:2019)
case_bucket_wt_years <- c(1998, 2003, 2008, 2013, 2019)

# Weighting years for the homeless and recent-entry case buckets.
smoothed_bucket_wt_years <- c(1994, 1999, 2004, 2009, 2014)

# Simulated years summed into each bucket of recent-entry cases.
recent_entry_bucket_years <- list(1994:1998, 1999:2003, 2004:2008, 2009:2013, 2014:2018)

# IGRA sensitivity and specificity for US-born people, non-US-born children
# under 5, and non-US-born people 5 and over.
igra_sens <- c(usb = 0.780, nusb_young = 0.712, nusb = 0.789)
igra_spec <- c(usb = 0.979, nusb_young = 0.989, nusb = 0.985)

#' Weighted normal log-likelihood, relative to a perfect fit
#'
#' @param obs observed values
#' @param sim simulated values on the scale of \code{obs}
#' @param sd standard deviations
#' @param w weights
#' @noRd
llik_normal <- function(obs, sim, sd, w) {
  sum(stats::dnorm(obs, sim, sd, log = TRUE) * w) - sum(stats::dnorm(obs, obs, sd, log = TRUE) * w)
}

#' Weighted Dirichlet-multinomial log-likelihood, relative to a perfect fit
#'
#' @param obs observed counts, one row per year or bucket
#' @param sim simulated values, one row per year or bucket; only their
#'   proportions within each row matter
#' @param rho correlation parameter
#' @param w weights, one per row
#' @param pad small count added to \code{obs} for the perfect fit, so that
#'   zero counts give finite densities
#' @noRd
llik_dirmult <- function(obs, sim, rho, w, pad) {
  sum(dDirMult(M = sim, n = obs, Rho = rho) * w) -
    sum(dDirMult(M = obs + pad, n = obs, Rho = rho) * w)
}

#' Weights for the given years
#' @noRd
year_wts <- function(wts, years) wts[as.character(years)]

#' Which years fall within the importance weights
#' @noRd
weighted_years <- function(years, wts) as.integer(years) %in% as.integer(names(wts))

#' Sum simulated values over each bucket of years
#'
#' @param sim annual simulated values, a vector or a matrix with one row per year
#' @param buckets list of year vectors
#' @return one value per bucket, or a matrix with one row per bucket
#' @noRd
sum_buckets <- function(sim, buckets) {
  if (is.null(dim(sim))) {
    return(vapply(buckets, function(years) sum(sim[year_to_idx(years)]), numeric(1)))
  }
  sums <- vapply(
    buckets,
    function(years) colSums(sim[year_to_idx(years), , drop = FALSE]),
    numeric(ncol(sim))
  )
  t(sums)
}

#' Combine the model's 85-94 and 95+ age groups into one 85+ column
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

#' Total cases each year
#'
#' Normal, centred on the reported cases with a 95% interval of +/- 10%. Years
#' with no reported cases are left out, and the most recent year counts double.
#'
#' @param sim_cases simulated cases, in millions
#' @noRd
llik_cases_total <- function(sim_cases, targets, wts) {
  d <- targets$cases_yr
  keep <- weighted_years(d[, "year"], wts)
  years <- as.integer(d[keep, "year"])
  obs <- as.numeric(d[keep, "cases"])
  w <- year_wts(wts, years)
  w[length(w)] <- 2 * w[length(w)]
  sim <- sim_cases[year_to_idx(years)] * 1e6
  reported <- !is.na(obs) & obs != 0
  llik_normal(obs[reported], sim[reported], obs[reported] * 0.1 / 1.96, w[reported])
}

#' Decline in total cases, 1953-1993
#'
#' The national year-to-year changes in cases, scaled to end at the location's
#' 1993 cases. Normal, with a 95% interval of +/- 20%.
#'
#' @param sim_cases simulated cases, in millions
#' @noRd
llik_cases_decline <- function(sim_cases, targets, wts) {
  years <- 1953:1993
  d <- targets$cases_yr
  cases_1993 <- as.numeric(d[d[, "year"] == 1993, "cases"])
  change <- targets$cases_prop_change_53_94
  obs <- cumprod(change) / prod(change) * cases_1993
  sim <- sim_cases[year_to_idx(years)] * 1e6
  llik_normal(obs, sim, obs * 0.2 / 1.96, year_wts(wts, years))
}

#' Cases by nativity, in 5-year buckets
#'
#' Normal, with a 95% interval of +/- 5%.
#'
#' @param sim_cases simulated cases in the nativity group, in millions
#' @param usb 1 for US-born, 0 for non-US-born
#' @noRd
llik_cases_nat_5yr <- function(sim_cases, usb, targets, wts) {
  d <- targets$cases_nat_5yr
  obs <- as.numeric(unlist(d[d$usb == usb, case_bucket_cols(d)]))
  sim <- sum_buckets(sim_cases, case_bucket_years) * 1e6
  llik_normal(obs, sim, obs * 0.05 / 1.96, year_wts(wts, case_bucket_wt_years))
}

#' Trend in cases by nativity across the two most recent 5-year buckets
#'
#' Normal on the change in log cases, with standard deviation 0.005.
#'
#' @param sim_nusb,sim_usb simulated non-US-born and US-born cases
#' @param usb_trend target for the US-born trend in place of the data's, or NA
#' @noRd
llik_cases_nat_trend <- function(sim_nusb, sim_usb, targets, usb_trend = NA) {
  d <- targets$cases_nat_5yr
  last_two <- utils::tail(case_bucket_cols(d), 2)
  obs <- c(
    diff(log(as.numeric(unlist(d[d$usb == 0, last_two])))),
    diff(log(as.numeric(unlist(d[d$usb == 1, last_two]))))
  )
  if (!is.na(usb_trend)) {
    obs[2] <- usb_trend
  }
  buckets <- utils::tail(case_bucket_years, 2)
  sim <- c(
    diff(log(sum_buckets(sim_nusb, buckets))),
    diff(log(sum_buckets(sim_usb, buckets)))
  )
  llik_normal(obs, sim, 0.005, 1)
}

#' Age distribution of cases in one nativity group, in 5-year buckets
#'
#' Dirichlet-multinomial. Age groups with suppressed counts are pooled into one
#' group holding the bucket's remaining cases.
#'
#' @param sim_by_age simulated cases in the nativity group, one column per model
#'   age group
#' @param usb 1 for US-born, 0 for non-US-born
#' @param rho correlation parameter
#' @noRd
llik_cases_age_5yr <- function(sim_by_age, usb, targets, wts, rho = 0.1) {
  d <- targets$cases_yr_ag_nat_5yr
  obs <- as.matrix(d[d[, "usb"] == usb, case_age_cols(d)])
  totals_d <- targets$cases_nat_5yr
  totals <- as.numeric(unlist(totals_d[totals_d$usb == usb, case_bucket_cols(totals_d)]))
  sim <- sum_buckets(merge_85_plus(sim_by_age), case_bucket_years)
  w <- year_wts(wts, case_bucket_wt_years)
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
  d <- targets$hr_cases_sm
  obs <- d$fract_sm * d$ss
  sim <- sum_buckets(sim_cases, case_bucket_years) * 1e6
  llik_normal(obs, sim, obs * 0.1 / 1.96, year_wts(wts, smoothed_bucket_wt_years))
}

#' Fraction of cases from recent transmission
#'
#' Beta with effective sample size 1000.
#'
#' @param sim_fraction simulated fraction of cases from recent transmission
#' @noRd
llik_recent_transmission <- function(sim_fraction, targets) {
  obs <- as.numeric(targets$rct_cases_sm$fract_sm)
  stats::dbeta(obs, sim_fraction * 1000, (1 - sim_fraction) * 1000, log = TRUE) -
    stats::dbeta(obs, obs * 1000, (1 - obs) * 1000, log = TRUE)
}

#' Non-US-born cases among recent entrants, in 5-year buckets
#'
#' The smoothed fraction of recent entrants times the non-US-born cases in each
#' bucket. Normal, with a 95% interval of +/- 10%.
#'
#' @param sim_cases simulated cases among recent entrants, in millions
#' @noRd
llik_cases_recent_entry <- function(sim_cases, targets, wts) {
  annual <- targets$cases_yr_ag_nat
  years <- annual[, "year", "nusb"]
  nusb <- vapply(
    recent_entry_bucket_years,
    function(bucket) sum(annual[years %in% bucket, "total", "nusb"]),
    numeric(1)
  )
  obs <- targets$rt_fb_cases_sm$p_recent * nusb
  sim <- sum_buckets(sim_cases, recent_entry_bucket_years) * 1e6
  llik_normal(obs, sim, obs * 0.1 / 1.96, year_wts(wts, smoothed_bucket_wt_years))
}

#' Treatment outcomes each year: completed, discontinued and died
#'
#' Dirichlet-multinomial.
#'
#' @param sim_outcomes simulated outcomes, three columns in that order
#' @param rho correlation parameter
#' @noRd
llik_tx_outcomes <- function(sim_outcomes, targets, wts, rho = 0.01) {
  d <- targets$tx_outcomes
  d <- d[weighted_years(d[, "year"], wts), , drop = FALSE]
  years <- as.integer(d[, "year"])
  share <- cbind(1 - rowSums(d[, c("pct_disc", "pct_dead")]), d[, "pct_disc"], d[, "pct_dead"])
  obs <- share * d[, "sample_size"]
  sim <- sim_outcomes[year_to_idx(years), , drop = FALSE]
  llik_dirmult(obs, sim, rho, year_wts(wts, years), pad = 0.01)
}

#' LTBI treatment starts in 2002
#'
#' Normal, centred on the estimate, with standard deviation set by its interval.
#'
#' @param sim_starts simulated treatment starts, in millions
#' @noRd
llik_ltbi_tx_total <- function(sim_starts, targets) {
  volume <- targets$TLTBI_volume
  llik_normal(volume[1], sim_starts * 1e6, diff(volume[2:3]) / 1.96, 1)
}

#' Fractions of LTBI treatment starts that are non-US-born and homeless, 2002
#'
#' Beta with effective sample size 100.
#'
#' @param sim_fractions simulated non-US-born and homeless fractions
#' @noRd
llik_ltbi_tx_dist <- function(sim_fractions, targets) {
  obs <- targets$TLTBI_dist[c("FB", "HR")]
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

#' IGRA-positive fraction by age in 2011
#'
#' Beta, from the survey's positive and negative counts in each age group from
#' 5 to 75+. Counts four times.
#'
#' @param tested expected positive and negative results by model age group
#' @param prev survey counts, columns \code{SSp} and \code{SSn}
#' @noRd
llik_ltbi_prev <- function(tested, prev) {
  tested[9, ] <- colSums(tested[9:11, ])
  ages <- 2:9
  sim <- tested[ages, 1] / rowSums(tested[ages, ])
  obs <- prev$SSp / (prev$SSp + prev$SSn)
  (sum(stats::dbeta(sim, prev$SSp, prev$SSn, log = TRUE)) -
     sum(stats::dbeta(obs, prev$SSp, prev$SSn, log = TRUE))) * 4
}

#' IGRA-positive fraction by age among US-born people, 2011
#'
#' @param sim_ltbi,sim_pop simulated US-born people with LTBI and in total, by
#'   model age group
#' @noRd
llik_ltbi_prev_usb <- function(sim_ltbi, sim_pop, targets) {
  tested <- igra_results(sim_ltbi, sim_pop - sim_ltbi, igra_sens[["usb"]], igra_spec[["usb"]])
  llik_ltbi_prev(tested, targets$LTBI_prev_US_11_IGRA)
}

#' IGRA-positive fraction by age among non-US-born people, 2011
#'
#' @param sim_ltbi,sim_pop simulated non-US-born people with LTBI and in total,
#'   by model age group
#' @noRd
llik_ltbi_prev_nusb <- function(sim_ltbi, sim_pop, targets) {
  no_ltbi <- sim_pop - sim_ltbi
  tested <- rbind(
    igra_results(sim_ltbi[1], no_ltbi[1], igra_sens[["nusb_young"]], igra_spec[["nusb_young"]]),
    igra_results(sim_ltbi[-1], no_ltbi[-1], igra_sens[["nusb"]], igra_spec[["nusb"]])
  )
  llik_ltbi_prev(tested, targets$LTBI_prev_FB_11_IGRA)
}

#' Deaths with TB each year
#'
#' Normal, with a 95% interval of +/- 10%. Years without data are left out.
#'
#' @param sim_deaths simulated deaths with TB by model age group, in millions
#' @noRd
llik_tb_deaths_total <- function(sim_deaths, targets, wts) {
  d <- targets$tbdeaths
  keep <- weighted_years(d$Year, wts)
  years <- as.integer(d$Year[keep])
  obs <- as.numeric(d$Deaths[keep])
  sim <- rowSums(sim_deaths[year_to_idx(years), , drop = FALSE]) * 1e6
  reported <- !is.na(obs)
  llik_normal(
    obs[reported], sim[reported], obs[reported] * 0.1 / 1.96, year_wts(wts, years)[reported]
  )
}

#' Age distribution of deaths with TB each year
#'
#' Dirichlet-multinomial.
#'
#' @param sim_deaths simulated deaths with TB by model age group
#' @param rho correlation parameter
#' @noRd
llik_tb_deaths_age <- function(sim_deaths, targets, wts, rho = 0.005) {
  d <- targets$tbdeaths_age_yr
  keep <- weighted_years(d$year, wts)
  years <- as.integer(d$year[keep])
  obs <- as.matrix(d[keep, -1])
  sim <- merge_85_plus(sim_deaths[year_to_idx(years), , drop = FALSE])
  llik_dirmult(obs, sim, rho, year_wts(wts, years), pad = 0.005)
}

#' Probability of not having developed TB after each time since infection
#'
#' Averaged over the model's progression risk groups.
#'
#' @param progression list of \code{Mpfast}, \code{Mrslow}, \code{rfast} and
#'   \code{rRecov} from the model parameters
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

#' Likelihood of a published study, or -10^4 when it cannot be evaluated
#'
#' @param value an expression giving the log-likelihood
#' @noRd
study_llik <- function(value) {
  llik <- tryCatch(value, error = function(e) -Inf)
  if (is.nan(llik) || llik == -Inf) -1e4 else llik
}

#' Borgdorff's cumulative TB incidence after infection
#'
#' Binomial on each year's share of cases, with effective sample size
#' 854.5463 / 4.
#'
#' @param progression see \code{tb_free_probability}
#' @noRd
llik_borgdorff <- function(progression, targets) {
  ss <- 854.5463 / 4
  d <- targets$borgdorff_data
  obs <- diff(-d[, 2])
  adj <- sum(obs * ss * log(obs) + (1 - obs) * ss * log(1 - obs))
  study_llik({
    p_free <- tb_free_probability(progression, d[, 1])
    sim <- diff(-(1 - (1 - p_free) / (1 - p_free)[nrow(d)]))
    sum(obs * ss * log(sim) + (1 - obs) * ss * log(1 - sim)) - adj
  })
}

#' Annual TB incidence after infection in a published cohort
#'
#' Binomial on each year's cases among those still TB-free, divided by 4.
#'
#' @param progression see \code{tb_free_probability}
#' @param d cohort data: year, people at risk, cases
#' @noRd
llik_cohort <- function(progression, d) {
  at_risk <- d[, 2]
  cases <- d[, 3]
  # Years with no cases have no perfect-fit rate, so count 0.01 cases there
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
  llik_cohort(progression, targets$ferebee_data)
}

#' Sutherland's TB incidence after infection
#' @noRd
llik_sutherland <- function(progression, targets) {
  llik_cohort(progression, targets$sutherland_data)
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

#' Years of the annual population table within the importance weights
#' @noRd
pop_years <- function(pop_00_19, wts) {
  years <- suppressWarnings(as.integer(colnames(pop_00_19)))
  sort(years[!is.na(years) & years %in% as.integer(names(wts))])
}

#' Total population of one nativity group: each census after the first, and
#' the latest weighted year of the annual table
#'
#' @param usb 1 for US-born, 0 for non-US-born
#' @return list of years, population, and the population of the last census,
#'   which sets the interval width
#' @noRd
pop_total_obs <- function(targets, usb, wts) {
  census <- targets$pop_50_10
  census <- census[census[, "usb"] == usb, , drop = FALSE]
  census <- census[order(census[, "year"]), , drop = FALSE][-1, , drop = FALSE]
  annual <- targets$pop_00_19
  latest <- max(pop_years(annual, wts))
  annual <- annual[annual[, "usb"] == usb & annual[, "age_group"] != 11, , drop = FALSE]
  list(
    years = c(as.integer(census[, "year"]), latest),
    pop = c(as.numeric(census[, "pop"]), sum(annual[, as.character(latest)])),
    scale = as.numeric(census[nrow(census), "pop"])
  )
}

#' Total population of one nativity group
#'
#' Normal, with a 95% interval of +/- 5% of the last census.
#'
#' @param sim_pop simulated population of the group, in millions
#' @param usb 1 for US-born, 0 for non-US-born
#' @noRd
llik_pop_total <- function(sim_pop, usb, targets, wts) {
  tg <- pop_total_obs(targets, usb, wts)
  sim <- sim_pop[year_to_idx(tg$years)] * 1e6
  llik_normal(tg$pop, sim, tg$scale * 0.05 / 1.96, year_wts(wts, tg$years))
}

#' Age distribution of the population by nativity, over the latest years
#'
#' Multinomial with effective sample size \code{ess}, on the population summed
#' over the latest \code{n_years} weighted years.
#'
#' @param sim_usb,sim_nusb simulated US-born and non-US-born population by model
#'   age group
#' @param ess effective sample size
#' @param n_years number of latest years pooled
#' @noRd
llik_pop_age <- function(sim_usb, sim_nusb, targets, wts, ess = 500, n_years = 3) {
  annual <- targets$pop_00_19
  years <- utils::tail(pop_years(annual, wts), n_years)
  annual <- annual[annual[, "age_group"] != 11, , drop = FALSE]
  pop_usb <- rowSums(annual[annual[, "usb"] == 1, as.character(years), drop = FALSE])
  pop_nusb <- rowSums(annual[annual[, "usb"] == 0, as.character(years), drop = FALSE])
  obs <- cbind(pop_usb / sum(pop_usb), pop_nusb / sum(pop_nusb))
  adj <- (sum(log(obs[, 1]) * obs[, 1]) + sum(log(obs[, 2]) * obs[, 2])) * ess
  sim <- cbind(
    colSums(sim_usb[year_to_idx(years), , drop = FALSE]),
    colSums(sim_nusb[year_to_idx(years), , drop = FALSE])
  )
  sim <- rbind(sim[1:9, ], sim[10, ] + sim[11, ])
  sim <- cbind(sim[, 1] / sum(sim[, 1]), sim[, 2] / sum(sim[, 2]))
  (sum(log(sim[, 1]) * obs[, 1]) + sum(log(sim[, 2]) * obs[, 2])) * ess - adj
}

#' Population experiencing homelessness in 2010
#'
#' Normal, with standard deviation set by the estimate's interval.
#'
#' @param sim_homeless simulated homeless population, on the scale of the target
#' @noRd
llik_homeless_pop <- function(sim_homeless, targets) {
  pop <- targets$homeless_pop
  llik_normal(pop[1], sim_homeless, diff(pop[2:3]) / 2 / 1.96, 1)
}

#' All-cause deaths in the latest weighted year
#'
#' Normal, with a 95% interval of +/- 10%.
#'
#' @param sim_deaths simulated deaths by model age group, in millions
#' @noRd
llik_deaths_total <- function(sim_deaths, targets, wts) {
  d <- targets$deaths_by_age
  years <- as.integer(d[, 1])
  year <- max(years[years %in% as.integer(names(wts))])
  obs <- d[match(year, d[, 1]), ncol(d)]
  sim <- rowSums(sim_deaths[year_to_idx(year), , drop = FALSE]) * 1e6
  llik_normal(obs, sim, obs * 0.1 / 1.96, year_wts(wts, year))
}

#' Age distribution of all-cause deaths in the latest two weighted years
#'
#' Dirichlet-multinomial.
#'
#' @param sim_deaths simulated deaths by model age group
#' @param rho correlation parameter
#' @noRd
llik_deaths_age <- function(sim_deaths, targets, wts, rho = 0.01) {
  d <- targets$deaths_by_age
  years <- as.integer(d[, 1])
  years <- utils::tail(sort(years[years %in% as.integer(names(wts))]), 2)
  obs <- d[match(years, d[, 1]), -c(1, ncol(d)), drop = FALSE]
  sim <- merge_85_plus(sim_deaths[year_to_idx(years), , drop = FALSE])
  llik_dirmult(obs, sim, rho, year_wts(wts, years), pad = 0.1)
}

#' Distribution of the population across mortality risk groups, 2015-2016
#'
#' Dirichlet-multinomial in each age group, against the model's assumed
#' distribution.
#'
#' @param sim_by_group simulated population with four mortality risk groups for
#'   each model age group, age group by age group
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
