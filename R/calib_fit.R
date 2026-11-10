# The fit of a simulation to one location's calibration targets: each target
# observation with the simulated value it is compared with. The likelihood
# terms are computed from this table, and the fit plots draw it.
#
# The fit table has one row per observation, with columns target, start, end,
# nativity, group, observed, simulated and total. Observed and simulated values
# are on one scale: counts of people for counts and fractions for fractions.
# Distributions are the exception: their simulated values are on the model's
# own scale, and only the proportions within each period are compared. total
# is the observed total an observation is part of: the people tested, the
# people at risk, or a bucket's cases where age groups are suppressed.
# Published studies and the Tiemersma estimates are not tied to calendar years;
# their start and end are NA, and group names the time or quantity.

#' Fit rows from observation rows and their simulated values
#' @noRd
fit_rows <- function(rows, simulated) {
  data.frame(
    target = rows$target, start = rows$start, end = rows$end, nativity = rows$nativity,
    group = rows$group, observed = rows$value, simulated = as.numeric(simulated),
    total = rows$total, stringsAsFactors = FALSE
  )
}

#' Fit rows for observations not tied to a location's data or to calendar years
#' @noRd
fixed_fit_rows <- function(target, group, observed, simulated, total = NA_real_,
                           start = NA_integer_, end = start) {
  data.frame(
    target = target, start = as.integer(start), end = as.integer(end), nativity = "all",
    group = group, observed = as.numeric(observed), simulated = as.numeric(simulated),
    total = as.numeric(total), stringsAsFactors = FALSE
  )
}

#' Fit rows for estimates
#' @noRd
estimate_fit_rows <- function(est, simulated) {
  fixed_fit_rows(est$target, est$group, est$estimate, simulated, start = est$start,
                 end = est$end)
}

#' Observation rows that fall within the simulated years
#'
#' @param n_years number of simulated years, from the model's base year
#' @noRd
simulated_rows <- function(d, n_years) d[year_to_idx(d$end) <= n_years, , drop = FALSE]

#' Fit rows for a target whose simulated values are an annual series summed
#' over each row's years
#'
#' @param sim annual simulated values, in millions
#' @param scale factor taking the simulated values to the target's scale
#' @noRd
period_fit <- function(targets, target, sim, nativity = NULL, scale = 1e6) {
  d <- simulated_rows(target_obs(targets, target, nativity), length(sim))
  fit_rows(d, sum_periods(sim, d$start, d$end) * scale)
}

#' Fit rows for a target observed by group, whose simulated values have one
#' column per group, summed over each row's years
#'
#' @param sim_by_group annual simulated values, one column per group named as in
#'   the target rows
#' @noRd
group_fit <- function(targets, target, sim_by_group, nativity = NULL, scale = 1) {
  d <- simulated_rows(target_obs(targets, target, nativity), nrow(sim_by_group))
  simulated <- vapply(seq_len(nrow(d)), function(i) {
    sum(sim_by_group[year_to_idx(d$start[i]:d$end[i]), d$group[i]])
  }, numeric(1))
  fit_rows(d, simulated * scale)
}

#' Simulated values by the data's ten age groups, with columns named by them
#' @noRd
by_data_age <- function(by_model_age) {
  merged <- merge_85_plus(by_model_age)
  colnames(merged) <- data_age_groups
  merged
}

#' Simulated outputs by model age group
#'
#' @param out annual simulation outputs with named columns
#' @param prefix output name before the age group, e.g. \code{"NOTIF_"}
#' @noRd
by_age <- function(out, prefix) out[, paste0(prefix, model_age_groups), drop = FALSE]

#' The fit of a simulation to a location's calibration targets
#'
#' @param out annual simulation outputs with named columns
#' @param prms model parameters from \code{param_init}
#' @param p model parameter vector
#' @param targets calibration targets, from \code{calib_targets}
#' @return the fit table described at the top of this file
#' @noRd
calib_fit <- function(out, prms, p, targets) {
  rbind(case_fit(out, targets), tb_fit(out, prms, p, targets), demographic_fit(out, targets))
}

#' Fit rows of the case targets
#' @noRd
case_fit <- function(out, targets) {
  cases <- out[, "NOTIF_ALL"] + out[, "NOTIF_MORT_ALL"]
  cases_usb_age <- by_age(out, "NOTIF_US_") + by_age(out, "NOTIF_US_MORT_")
  cases_nusb_age <- (by_age(out, "NOTIF_") + by_age(out, "NOTIF_MORT_")) - cases_usb_age
  cases_usb <- out[, "NOTIF_US"] + out[, "NOTIF_MORT_US"]
  cases_nusb <- out[, "NOTIF_F1"] + out[, "NOTIF_F2"] +
    (out[, "NOTIF_MORT_F1"] + out[, "NOTIF_MORT_F2"])
  nativity_usb <- period_fit(targets, "cases_nativity", cases_usb, "usb")
  nativity_nusb <- period_fit(targets, "cases_nativity", cases_nusb, "nusb")
  rbind(
    period_fit(targets, "cases_total", cases),
    period_fit(targets, "cases_decline", cases),
    nativity_usb,
    nativity_nusb,
    with_bucket_totals(
      group_fit(targets, "cases_age", by_data_age(cases_usb_age), "usb", scale = 1e6),
      nativity_usb
    ),
    with_bucket_totals(
      group_fit(targets, "cases_age", by_data_age(cases_nusb_age), "nusb", scale = 1e6),
      nativity_nusb
    ),
    period_fit(targets, "cases_homeless", out[, "NOTIF_HR"] + out[, "NOTIF_MORT_HR"]),
    period_fit(targets, "cases_recent_entry", out[, "NOTIF_F1"] + out[, "NOTIF_MORT_F1"])
  )
}

#' Fit rows with each bucket's observed total from another target's fit rows
#' for the same buckets
#' @noRd
with_bucket_totals <- function(f, totals) {
  f$total <- totals$observed[match(f$start, totals$start)]
  f
}

#' Fit rows of the TB targets other than cases
#' @noRd
tb_fit <- function(out, prms, p, targets) {
  window <- target_estimate(targets, "recent_transmission")
  window_rows <- year_to_idx(window$start:window$end)
  recent <- sum(out[window_rows, c("INCID_REC_US", "INCID_REC_FB")]) /
    sum(out[window_rows, c("INCID_ALL_US", "INCID_ALL_FB")])
  ltbi_tx_total <- target_estimate(targets, "ltbi_tx_total")
  ltbi_tx_share <- target_estimate(targets, "ltbi_tx_share", c("nusb", "homeless"))
  ltbi_tx_row <- year_to_idx(ltbi_tx_total$start)
  starts <- out[ltbi_tx_row, "TLTBI_INITS"]
  share <- out[ltbi_tx_row, c("TLTBI_INITS_FB", "TLTBI_INITS_HR")] / starts
  outcomes <- out[, c("TBTX_COMPLT", "TBTX_DISCONT", "TBTX_DIED")]
  colnames(outcomes) <- tx_outcome_groups
  tb_deaths <- by_age(out, "TOTMORT_W_TB_")
  progression <- prms[c("Mpfast", "Mrslow", "rfast", "rRecov")]
  rbind(
    estimate_fit_rows(window, recent),
    group_fit(targets, "tx_outcomes", outcomes),
    estimate_fit_rows(ltbi_tx_total, starts * 1e6),
    estimate_fit_rows(ltbi_tx_share, share),
    ltbi_prev_fit(out, targets, "usb"),
    ltbi_prev_fit(out, targets, "nusb"),
    period_fit(targets, "tb_deaths_total", rowSums(tb_deaths)),
    group_fit(targets, "tb_deaths_age", by_data_age(tb_deaths)),
    borgdorff_fit(progression, targets$studies$borgdorff_data),
    cohort_fit("ferebee", progression, targets$studies$ferebee_data),
    cohort_fit("sutherland", progression, targets$studies$sutherland_data),
    tiemersma_fit(p[["rSlfCur"]], p[["muIp"]])
  )
}

#' Fit rows of the demographic targets
#' @noRd
demographic_fit <- function(out, targets) {
  pop <- list(usb = out[, "N_US"], nusb = out[, "N_F1"] + out[, "N_F2"])
  pop_by_age <- list(usb = by_age(out, "N_US_"), nusb = by_age(out, "N_FB_"))
  homeless <- target_estimate(targets, "homeless_pop")
  deaths <- by_age(out, "TOTMORT_")
  per_nativity <- lapply(c("usb", "nusb"), function(nativity) {
    rbind(
      period_fit(targets, "population_census", pop[[nativity]], nativity),
      population_total_fit(targets, pop[[nativity]], nativity),
      group_fit(targets, "population", by_data_age(pop_by_age[[nativity]]), nativity)
    )
  })
  rbind(
    do.call(rbind, per_nativity),
    estimate_fit_rows(homeless, out[year_to_idx(homeless$start), "N_HR"]),
    period_fit(targets, "deaths_total", rowSums(deaths)),
    group_fit(targets, "deaths_age", by_data_age(deaths)),
    mort_dist_fit(out[, paste0("N_ag_", rep(seq_along(model_age_groups), each = 4), "_M", 1:4)])
  )
}

#' Fit rows of a nativity group's total population each year: the annual
#' population summed over age groups, against the simulated total
#' @noRd
population_total_fit <- function(targets, sim_pop, nativity) {
  d <- simulated_rows(target_obs(targets, "population", nativity), length(sim_pop))
  years <- sort(unique(d$start))
  totals <- vapply(years, function(year) sum(d$value[d$start == year]), numeric(1))
  rows <- obs_rows("population_total", years, value = totals, nativity = nativity)
  fit_rows(rows, sim_pop[year_to_idx(years)] * 1e6)
}

#' Fit rows of IGRA survey results for one nativity group
#'
#' The simulated value of each survey age group is the expected number of
#' positive results among the people tested.
#' @noRd
ltbi_prev_fit <- function(out, targets, nativity) {
  d <- target_obs(targets, "ltbi_prev", nativity)
  survey_year <- year_to_idx(unique(d$start))
  prefix <- c(usb = "N_US_", nusb = "N_FB_")[[nativity]]
  ltbi <- by_age(out, paste0(prefix, "LTBI_"))[survey_year, ]
  pop <- by_age(out, prefix)[survey_year, ]
  tested <- if (nativity == "usb") {
    igra_results(ltbi, pop - ltbi, igra_sens[["usb"]], igra_spec[["usb"]])
  } else {
    rbind(
      igra_results(ltbi[1], pop[1] - ltbi[1], igra_sens[["nusb_young"]], igra_spec[["nusb_young"]]),
      igra_results(ltbi[-1], pop[-1] - ltbi[-1], igra_sens[["nusb"]], igra_spec[["nusb"]])
    )
  }
  fit_rows(d, igra_positive_fraction(tested) * d$total)
}

#' Fit rows of Borgdorff's cumulative TB incidence after infection: each
#' interval's share of cases, NaN when the curve cannot be computed
#' @noRd
borgdorff_fit <- function(progression, d) {
  simulated <- tryCatch({
    p_free <- tb_free_probability(progression, d[, "Years"])
    stopifnot(length(p_free) == nrow(d))
    diff(-(1 - (1 - p_free) / (1 - p_free)[nrow(d)]))
  }, error = function(e) rep(NaN, nrow(d) - 1))
  fixed_fit_rows("borgdorff", as.character(d[-1, "Years"]), diff(-d[, "mean"]), simulated)
}

#' Fit rows of a published cohort's annual TB incidence after infection: each
#' year's cases and people at risk, against the cases expected at the
#' simulated rate, NaN when the rates cannot be computed
#' @noRd
cohort_fit <- function(target, progression, d) {
  n_years <- nrow(d)
  simulated <- tryCatch({
    p_free <- tb_free_probability(progression, 0:n_years)
    stopifnot(length(p_free) == n_years + 1)
    -log(1 - diff(-p_free) / p_free[-(n_years + 1)])
  }, error = function(e) rep(NaN, n_years))
  fixed_fit_rows(target, as.character(d[, 1]), d[, 3], simulated * d[, "At_Risk"],
                 total = d[, "At_Risk"])
}

#' Fit rows of Tiemersma's duration of untreated TB and its case fatality
#'
#' @param self_cure,mortality monthly rates of self-cure and death from
#'   untreated TB
#' @noRd
tiemersma_fit <- function(self_cure, mortality) {
  self_cure <- self_cure * 12
  mortality <- mortality * 12
  fixed_fit_rows(
    "tiemersma", c("duration", "fatality"), c(3, 0.9),
    c(1 / (self_cure + mortality), mortality / (self_cure + mortality))
  )
}

# Years over which the population's mortality risk distribution is compared
# with the model's assumed distribution.
mort_dist_years <- 2015:2016

#' Fit rows of the population's distribution across mortality risk groups, in
#' each model age group, against the distribution in the package constant
#' \code{dist_gen}
#'
#' @param sim_by_group simulated population in millions, four mortality risk
#'   group columns for each model age group in turn
#' @noRd
mort_dist_fit <- function(sim_by_group) {
  rows <- year_to_idx(mort_dist_years)
  shares <- rowSums(dist_gen)
  n_groups <- length(shares)
  n_years <- length(mort_dist_years)
  fixed_fit_rows(
    "mort_dist",
    group = rep(paste0(rep(model_age_groups, each = n_groups), "_M", seq_len(n_groups)),
                each = n_years),
    observed = rep(rep(shares, length(model_age_groups)), each = n_years),
    simulated = as.vector(sim_by_group[rows, ]) * 1e6,
    start = rep(mort_dist_years, times = ncol(sim_by_group))
  )
}
