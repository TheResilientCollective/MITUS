# Calibration setup for three calibrated parameters, one from each prior
# family, and one background parameter.
fake_setup <- function() {
  prior <- matrix(0, 3, 8)
  prior[, 6:7] <- rbind(c(2, 3), c(2, 4), c(1, 2))
  prior[, 8] <- c(0.5, 0.5, 0.5)
  list(
    P = c(a = 0, b = 0, c = 0, background = 5),
    ii = c(TRUE, TRUE, TRUE, FALSE),
    ParamInitZ = prior,
    idZ0 = c(TRUE, FALSE, FALSE),
    idZ1 = c(FALSE, TRUE, FALSE),
    idZ2 = c(FALSE, FALSE, TRUE),
    wts = flat_wts()
  )
}

test_that("calib_params maps each prior family from the normal scale", {
  p <- calib_params(c(0, 1, -1), fake_setup())
  expect_equal(p[["a"]], stats::qbeta(0.5, 2, 3))
  expect_equal(p[["b"]], stats::qgamma(stats::pnorm(1), shape = 2, rate = 4))
  expect_equal(p[["c"]], -1)
  expect_identical(p[["background"]], 5)
})

test_that("the 2020 pandemic adjusts diagnosis, transmission and TB mortality", {
  month <- function(year, m) year_month_to_idx(year, m)
  prms <- list(rDxt = matrix(2, month(2100, 12), 3), NixTrans = rep(1, month(2100, 12)))
  adjusted <- apply_pandemic_2020(prms)
  rdxt <- adjusted$rDxt[, 1]
  trans <- adjusted$NixTrans
  rr <- adjusted$RRmuTBPand
  lockdown_rdxt <- 2 * (1 - 0.3707595)
  lockdown_trans <- 1 - 0.1984619

  expect_identical(c(rdxt[month(2020, 2)], trans[month(2020, 2)], rr[month(2020, 2)]), c(2, 1, 1))
  expect_equal(c(rdxt[month(2020, 3)], rdxt[month(2021, 12)]), rep(lockdown_rdxt, 2))
  expect_equal(c(trans[month(2020, 3)], trans[month(2021, 12)]), rep(lockdown_trans, 2))
  expect_equal(c(rr[month(2020, 3)], rr[month(2021, 12)]), rep(1.1158255, 2))

  recovery <- month(2022, 1):month(2023, 12)
  expect_equal(rdxt[recovery], seq(lockdown_rdxt, 2, length.out = 24))
  expect_equal(trans[recovery], seq(lockdown_trans, 1, length.out = 24))
  expect_equal(rr[recovery], seq(1.1158255, 1, length.out = 24))
  expect_identical(c(rdxt[month(2024, 1)], trans[month(2024, 1)], rr[month(2024, 1)]), c(2, 1, 1))
  expect_identical(adjusted$rDxt[, 3], adjusted$rDxt[, 1])
})

test_that("a valid simulation gives the sum of its terms", {
  local_mocked_bindings(calib_llik_terms = function(...) c(a = -1, b = -2))
  expect_identical(calib_llikelihood_one(c(0, 0, 0), "CA", list(), fake_setup(), TRUE, 2019), -3)
})

test_that("an invalid simulation gives the invalid likelihood", {
  local_mocked_bindings(calib_llik_terms = function(...) NULL)
  expect_identical(
    calib_llikelihood_one(c(0, 0, 0), "CA", list(), fake_setup(), TRUE, 2019),
    invalid_llik
  )
})

test_that("a likelihood that cannot be evaluated gives the distance penalty", {
  par <- c(1, 0.5, 0)
  penalty <- invalid_llik - sum((c(0.5, 0.5, 0.5) - par)^2)
  local_mocked_bindings(calib_llik_terms = function(...) stop("simulation failed"))
  expect_message(
    expect_identical(calib_llikelihood_one(par, "CA", list(), fake_setup(), TRUE, 2019), penalty),
    "could not be evaluated"
  )
  local_mocked_bindings(calib_llik_terms = function(...) c(a = -1, b = NaN))
  expect_identical(
    suppressMessages(calib_llikelihood_one(par, "CA", list(), fake_setup(), TRUE, 2019)),
    penalty
  )
  local_mocked_bindings(calib_llik_terms = function(...) c(a = -Inf))
  expect_identical(
    suppressMessages(calib_llikelihood_one(par, "CA", list(), fake_setup(), TRUE, 2019)),
    penalty
  )
})

# Run calib_llikelihood with a fake setup, scoring each parameter set by the sum
# of its values.
summed_llik <- function(par_matrix, n_cores = 1) {
  local_mocked_bindings(
    resolve_calib_end_year = function(calib_end_year) 2019,
    check_target_horizon = function(targets, wts) NULL,
    calib_setup = function() fake_setup(),
    calib_llikelihood_one = function(par, ...) sum(par)
  )
  calib_llikelihood(par_matrix, "CA", list(), n_cores = n_cores)
}

test_that("calib_llikelihood takes one parameter set in any shape, or a matrix of them", {
  expect_identical(summed_llik(c(a = 1, b = 2, c = 3)), 6)
  expect_identical(summed_llik(matrix(c(1, 2, 3), ncol = 1)), 6)
  expect_identical(summed_llik(data.frame(a = 1, b = 2, c = 3)), 6)
  expect_identical(summed_llik(rbind(c(1, 2, 3), c(4, 5, 6))), c(6, 15))
})

test_that("calib_llikelihood needs one value per calibrated parameter", {
  expect_error(summed_llik(c(1, 2)), "Each parameter set needs 3 values, not 2")
})

test_that("a parameter set that fails in parallel stops the calibration", {
  local_mocked_bindings(
    resolve_calib_end_year = function(calib_end_year) 2019,
    check_target_horizon = function(targets, wts) NULL,
    calib_setup = function() fake_setup(),
    calib_llikelihood_one = function(par, ...) if (par[1] > 1) stop("worker failed") else sum(par)
  )
  expect_error(
    suppressWarnings(
      calib_llikelihood(rbind(c(1, 2, 3), c(4, 5, 6)), "CA", list(), n_cores = 2)
    ),
    "failed for parameter set 2.*worker failed"
  )
})

test_that("calib_llikelihood passes the TB switch and end year to every set", {
  seen <- new.env()
  local_mocked_bindings(
    resolve_calib_end_year = function(calib_end_year) calib_end_year,
    check_target_horizon = function(targets, wts) NULL,
    calib_setup = function() fake_setup(),
    calib_llik_terms = function(p, loc, targets, wts, tb, calib_end_year) {
      seen$tb <- tb
      seen$calib_end_year <- calib_end_year
      c(a = 1)
    }
  )
  calib_llikelihood(c(0, 0, 0), "CA", list(), tb = FALSE, calib_end_year = 2021)
  expect_identical(seen$tb, FALSE)
  expect_identical(seen$calib_end_year, 2021)
})

test_that("calib_llik_terms keeps only the demographic terms without TB", {
  local_globals(list(wts = flat_wts()))
  seen <- new.env()
  local_mocked_bindings(
    simulate_calibration = function(p, loc, calib_end_year) {
      seen$calib_end_year <- calib_end_year
      list(out = fake_outputs(), prms = fake_progression())
    }
  )
  terms <- calib_llik_terms(c(rSlfCur = 0.01, muIp = 0.02), "CA", ca_targets(), flat_wts(),
                            tb = FALSE, calib_end_year = 2021)
  expect_named(terms, c("pop_nusb", "pop_usb", "pop_age", "homeless_pop", "deaths_total",
                        "deaths_age", "mort_dist"))
  expect_identical(seen$calib_end_year, 2021L)
})

test_that("weights ending before 2019 are an error", {
  expect_error(check_weight_years(flat_wts()[as.character(1950:2018)]),
               "importance weights end in 2018")
  expect_silent(check_weight_years(flat_wts()))
})

test_that("calib_setup reads the globals model_load sets, and needs all of them", {
  setup <- fake_setup()
  local_globals(setup)
  expect_identical(calib_setup(), setup)
  rm("idZ2", envir = globalenv())
  expect_error(calib_setup(), "Run model_load\\(\\) before calibrating; missing idZ2$")
})

# Simulation run from cSim with every annual output set to 1, after applying
# `change` to the outputs matrix.
local_simulation <- function(change = identity, v1 = 1, env = parent.frame()) {
  outputs <- func_ResNam()
  local_mocked_bindings(
    param_init = function(...) list(RRmuHR = c(1, 1), ResNam = outputs),
    def_prgchng = function(...) NULL,
    def_ttt = function(...) NULL,
    apply_pandemic_2020 = function(prms) prms,
    reblncd = function(...) matrix(0),
    cSim = function(...) {
      list(Outputs = change(matrix(1, 70, length(outputs))), V1 = array(v1, c(2, 2)))
    },
    .env = env
  )
}

test_that("a simulation with negative or missing 2017 outputs or state is invalid", {
  set_2017 <- function(value) function(out) replace(out, cbind(year_to_idx(2017), 5), value)
  local_simulation(set_2017(-1))
  expect_null(suppressMessages(simulate_calibration(1, "CA", 2019)))
  local_simulation(set_2017(NA))
  expect_null(suppressMessages(simulate_calibration(1, "CA", 2019)))
  local_simulation(v1 = -1)
  expect_null(suppressMessages(simulate_calibration(1, "CA", 2019)))
})

test_that("a valid simulation returns named outputs, and only 2017 is checked", {
  local_simulation(function(out) replace(out, cbind(year_to_idx(2016), 5), -1))
  sim <- simulate_calibration(1, "CA", 2019)
  expect_identical(colnames(sim$out), func_ResNam())
  expect_identical(unname(sim$out[year_to_idx(2016), 5]), -1)
})

tb_terms <- function(loc, out = fake_outputs(), targets = ca_targets()) {
  tb_llik_terms(out, fake_progression(), c(rSlfCur = 0.01, muIp = 0.02), loc, targets,
                flat_wts())
}

all_tb_terms <- c(
  "cases_total", "cases_decline", "cases_age_usb", "cases_age_nusb", "cases_nusb", "cases_usb",
  "cases_nativity_trend", "cases_homeless", "recent_transmission", "cases_recent_entry",
  "tx_outcomes", "ltbi_tx_total", "ltbi_tx_dist", "ltbi_prev_usb", "ltbi_prev_nusb",
  "tb_deaths_total", "tb_deaths_age", "borgdorff", "ferebee", "sutherland", "tiemersma"
)

test_that("TB terms are all present and finite for a typical state", {
  terms <- tb_terms("CA")
  expect_named(terms, all_tb_terms)
  expect_true(all(is.finite(terms)))
})

test_that("TB terms leave out the targets each listed state drops", {
  expect_named(tb_terms("MT"), setdiff(all_tb_terms, "cases_nativity_trend"))
  for (loc in c("RI", "ID", "NH", "VT", "WY")) {
    expect_named(tb_terms(loc), setdiff(all_tb_terms, "cases_homeless"), info = loc)
  }
})

test_that("each TB term gets the simulated values it compares", {
  out <- fake_outputs()
  targets <- ca_targets()
  terms <- tb_terms("CA", out)
  nusb_age <- by_age(out, "NOTIF_") + by_age(out, "NOTIF_MORT_") -
    (by_age(out, "NOTIF_US_") + by_age(out, "NOTIF_US_MORT_"))
  expect_identical(
    terms[["cases_age_nusb"]], llik_cases_age_5yr(nusb_age, "nusb", targets, flat_wts())
  )
  rows <- year_to_idx(2015:2018)
  fraction <- sum(out[rows, c("INCID_REC_US", "INCID_REC_FB")]) /
    sum(out[rows, c("INCID_ALL_US", "INCID_ALL_FB")])
  expect_identical(terms[["recent_transmission"]], llik_recent_transmission(fraction, targets))
  y2011 <- year_to_idx(2011)
  expect_identical(
    terms[["ltbi_prev_usb"]],
    llik_ltbi_prev_usb(by_age(out, "N_US_LTBI_")[y2011, ], by_age(out, "N_US_")[y2011, ], targets)
  )
  expect_identical(terms[["tiemersma"]], llik_tiemersma(self_cure = 0.01, mortality = 0.02))
})

test_that("single-year and fixed-window targets take their years from the targets", {
  out <- fake_outputs()
  targets <- ca_targets()
  moved <- targets
  est <- moved$estimates
  est[est$target == "recent_transmission", c("start", "end")] <- c(2012L, 2016L)
  est[est$target %in% c("ltbi_tx_total", "ltbi_tx_share"), c("start", "end")] <- 2005L
  est[est$target == "homeless_pop", c("start", "end")] <- 2013L
  moved$estimates <- est
  moved$obs[moved$obs$target == "ltbi_prev", c("start", "end")] <- 2014L
  terms <- tb_terms("CA", out, moved)
  rows <- year_to_idx(2012:2016)
  fraction <- sum(out[rows, c("INCID_REC_US", "INCID_REC_FB")]) /
    sum(out[rows, c("INCID_ALL_US", "INCID_ALL_FB")])
  expect_identical(terms[["recent_transmission"]], llik_recent_transmission(fraction, moved))
  y2005 <- year_to_idx(2005)
  expect_identical(terms[["ltbi_tx_total"]], llik_ltbi_tx_total(out[y2005, "TLTBI_INITS"], moved))
  y2014 <- year_to_idx(2014)
  expect_identical(
    terms[["ltbi_prev_nusb"]],
    llik_ltbi_prev_nusb(by_age(out, "N_FB_LTBI_")[y2014, ], by_age(out, "N_FB_")[y2014, ], moved)
  )
  expect_identical(
    demographic_llik_terms(out, moved, flat_wts())[["homeless_pop"]],
    llik_homeless_pop(out[year_to_idx(2013), "N_HR"], moved)
  )
})

test_that("North Dakota's US-born case trend target is fixed", {
  out <- fake_outputs()
  usb <- out[, "NOTIF_US"] + out[, "NOTIF_MORT_US"]
  nusb <- out[, "NOTIF_F1"] + out[, "NOTIF_F2"] + (out[, "NOTIF_MORT_F1"] + out[, "NOTIF_MORT_F2"])
  expect_identical(
    tb_terms("ND", out)[["cases_nativity_trend"]],
    llik_cases_nativity_trend(usb, nusb, ca_targets(), usb_trend = -0.1)
  )
  expect_identical(
    tb_terms("CA", out)[["cases_nativity_trend"]],
    llik_cases_nativity_trend(usb, nusb, ca_targets())
  )
})

test_that("demographic terms are named and finite", {
  terms <- demographic_llik_terms(fake_outputs(), ca_targets(), flat_wts())
  expect_named(terms, c("pop_nusb", "pop_usb", "pop_age", "homeless_pop", "deaths_total",
                        "deaths_age", "mort_dist"))
  expect_true(all(is.finite(terms)))
})

test_that("the homeless population term compares the 2010 population", {
  out <- fake_outputs()
  targets <- ca_targets()
  expect_identical(
    demographic_llik_terms(out, targets, flat_wts())[["homeless_pop"]],
    llik_homeless_pop(out[year_to_idx(2010), "N_HR"], targets)
  )
})

# The expected values are what the earlier state likelihood, llikelihood_st,
# gives for the same parameters and targets.
test_that("California's calibrated parameters give the state likelihood", {
  skip_on_cran()
  local_model_load("CA")
  optim <- readRDS(system.file("CA", "CA_Optim_all_15_0708.rds", package = "MITUS"))
  par <- optim[which.min(optim[, ncol(optim)]), seq_len(nrow(calib_setup()$ParamInitZ))]
  targets <- ca_targets()
  expect_equal(calib_llikelihood(par, "CA", targets), -335.3962908092, tolerance = 1e-10)
  expect_equal(calib_llikelihood(par, "CA", targets, tb = FALSE), -18.8659781051,
               tolerance = 1e-10)
})
