# Calibration setup for three calibrated parameters, one from each prior
# family, and one background parameter.
fake_setup <- function() {
  prior <- matrix(0, 3, 8)
  prior[, 6:7] <- rbind(c(2, 3), c(2, 1), c(1, 2))
  prior[, 8] <- c(0.5, 0.5, 0.5)
  list(
    P = c(a = 0, b = 0, c = 0, background = 5),
    ii = c(TRUE, TRUE, TRUE, FALSE),
    ParamInitZ = prior,
    idZ0 = c(TRUE, FALSE, FALSE),
    idZ1 = c(FALSE, TRUE, FALSE),
    idZ2 = c(FALSE, FALSE, TRUE),
    wts = setNames(rep(1, 70), 1950:2019)
  )
}

test_that("calib_params maps each prior family from the normal scale", {
  p <- calib_params(c(0, 1, -1), fake_setup())
  expect_equal(p[["a"]], stats::qbeta(0.5, 2, 3))
  expect_equal(p[["b"]], stats::qgamma(stats::pnorm(1), shape = 2, rate = 1))
  expect_equal(p[["c"]], -1)
  expect_identical(p[["background"]], 5)
})

test_that("the 2020 pandemic lowers diagnosis and transmission, then restores them", {
  months <- year_month_to_idx(2100, 12)
  prms <- list(rDxt = matrix(2, months, 3), NixTrans = rep(1, months))
  adjusted <- apply_pandemic_2020(prms)
  lockdown <- year_month_to_idx(2021, 6)
  expect_equal(adjusted$rDxt[lockdown, ], rep(2 * (1 - pandemic_2020[["Dxt"]]), 3))
  expect_equal(adjusted$NixTrans[lockdown], 1 - pandemic_2020[["Trans"]])
  expect_equal(adjusted$rDxt[year_month_to_idx(2023, 12), ], rep(2, 3))
  expect_identical(adjusted$rDxt[year_month_to_idx(2024, 1), ], rep(2, 3))
  expect_identical(adjusted$rDxt[year_month_to_idx(2020, 2), ], rep(2, 3))
  rr <- adjusted$RRmuTBPand
  expect_identical(rr[year_month_to_idx(2020, 2)], 1)
  expect_equal(rr[lockdown], pandemic_2020[["CaseFat"]])
  expect_equal(rr[year_month_to_idx(2023, 12)], 1)
  expect_length(rr, months)
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

test_that("calib_llikelihood takes one parameter set or a matrix of them", {
  local_mocked_bindings(
    resolve_calib_end_year = function(calib_end_year) 2019,
    calib_setup = function() fake_setup(),
    calib_llikelihood_one = function(par, ...) sum(par)
  )
  expect_identical(calib_llikelihood(c(a = 1, b = 2, c = 3), "CA", list()), 6)
  expect_identical(calib_llikelihood(rbind(c(1, 2, 3), c(4, 5, 6)), "CA", list()), c(6, 15))
})

test_that("calib_setup reads the globals model_load sets, and needs all of them", {
  g <- globalenv()
  setup <- fake_setup()
  saved <- mget(names(setup), envir = g, ifnotfound = list(NULL))
  withr::defer({
    for (name in names(saved)) {
      if (is.null(saved[[name]])) {
        rm(list = intersect(name, ls(g)), envir = g)
      } else {
        assign(name, saved[[name]], envir = g)
      }
    }
  })
  list2env(setup, envir = g)
  expect_identical(calib_setup(), setup)
  rm("idZ2", envir = g)
  expect_error(calib_setup(), "Run model_load\\(\\) before calibrating; missing idZ2$")
})

# Annual outputs with every simulated quantity, filled with positive values.
fake_outputs <- function() {
  set.seed(3)
  names <- func_ResNam()
  matrix(stats::runif(70 * length(names), 1, 2) * 1e-3, 70, length(names),
         dimnames = list(NULL, names))
}

fake_prms <- list(
  Mpfast = matrix(0.1, 4, 4),
  Mrslow = matrix(0.0002, 4, 4),
  rfast = 0.1,
  rRecov = 0.002
)

tb_terms <- function(loc, out, targets) {
  suppressWarnings(tb_llik_terms(
    out, fake_prms, c(rSlfCur = 0.01, muIp = 0.02), loc, targets, setNames(rep(1, 70), 1950:2019)
  ))
}

test_that("TB terms leave out the targets some states drop", {
  targets <- suppressWarnings(calib_targets("CA"))
  out <- fake_outputs()
  all_terms <- c(
    "cases_total", "cases_decline", "cases_age_usb", "cases_age_nusb", "cases_nusb", "cases_usb",
    "cases_nat_trend", "cases_homeless", "recent_transmission", "cases_recent_entry",
    "tx_outcomes", "ltbi_tx_total", "ltbi_tx_dist", "ltbi_prev_usb", "ltbi_prev_nusb",
    "tb_deaths_total", "tb_deaths_age", "borgdorff", "ferebee", "sutherland", "tiemersma"
  )
  expect_named(tb_terms("CA", out, targets), all_terms)
  expect_named(tb_terms("MT", out, targets), setdiff(all_terms, "cases_nat_trend"))
  expect_named(tb_terms("VT", out, targets), setdiff(all_terms, "cases_homeless"))
})

test_that("North Dakota's US-born case trend target is fixed", {
  targets <- suppressWarnings(calib_targets("CA"))
  out <- fake_outputs()
  nusb <- out[, "NOTIF_F1"] + out[, "NOTIF_F2"] + (out[, "NOTIF_MORT_F1"] + out[, "NOTIF_MORT_F2"])
  usb <- out[, "NOTIF_US"] + out[, "NOTIF_MORT_US"]
  expect_identical(
    tb_terms("ND", out, targets)[["cases_nat_trend"]],
    llik_cases_nat_trend(nusb, usb, targets, usb_trend = -0.1)
  )
  expect_identical(
    tb_terms("CA", out, targets)[["cases_nat_trend"]],
    llik_cases_nat_trend(nusb, usb, targets)
  )
})

test_that("demographic terms are named and finite", {
  targets <- suppressWarnings(calib_targets("CA"))
  terms <- demographic_llik_terms(fake_outputs(), targets, setNames(rep(1, 70), 1950:2019))
  expect_named(terms, c("pop_nusb", "pop_usb", "pop_age", "homeless_pop", "deaths_total",
                        "deaths_age", "mort_dist"))
  expect_true(all(is.finite(terms)))
})

test_that("California's calibrated parameters give the state likelihood", {
  skip_on_cran()
  suppressWarnings(model_load("CA"))
  optim <- readRDS(system.file("CA", "CA_Optim_all_15_0708.rds", package = "MITUS"))
  par <- optim[which.min(optim[, ncol(optim)]), seq_len(nrow(calib_setup()$ParamInitZ))]
  targets <- suppressWarnings(calib_targets("CA"))
  expect_equal(calib_llikelihood(par, "CA", targets), -352.4801406337, tolerance = 1e-10)
})
