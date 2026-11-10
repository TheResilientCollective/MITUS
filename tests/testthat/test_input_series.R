test_that("mortality is forecast from 2017 by fixed yearly declines", {
  mortality <- cbind(1950:2100, matrix(1, 151, 11))
  mortality[year_to_idx(2016), -1] <- 2
  forecast <- forecast_mortality(mortality)
  expect_identical(dim(forecast), c(151L, 11L))
  expect_identical(forecast[year_to_idx(2016), ], rep(2, 11))
  expect_equal(forecast[year_to_idx(2017), 1:10], 2 * (1 - mortality_declines[2:11]))
  expect_equal(forecast[year_to_idx(2018), 1], 2 * (1 - 0.0159)^2)
  expect_identical(forecast[year_to_idx(2017):151, 11], rep(1, 84))
})

test_that("treatment inputs hold their last value, named by year", {
  held <- held_to_end(cbind(1993:1995, c(0.1, 0.2, 0.3)))
  expect_identical(names(held), as.character(1993:2100))
  expect_identical(unname(held[c("1995", "2100")]), c(0.3, 0.3))
})

test_that("the immigrant TB trend declines after its data and doubles for ND", {
  inputs <- list(ImmigInputs = list(TBBurdenImmig = seq(2, 1, length.out = 69)))
  trend <- immigrant_tb_trend(inputs, "CA")
  expect_identical(names(trend), as.character(1950:2100))
  expect_equal(trend[["2019"]], trend[["2018"]] * 0.985)
  nd <- immigrant_tb_trend(inputs, "ND")
  expect_equal(nd[["2020"]], 2 * nd[["2010"]])
  expect_identical(nd[["2009"]], trend[["2009"]])
})

test_that("model_load builds the location's input series", {
  local_model_load("CA")
  series <- input_series()
  expect_identical(series$loc, "CA")
  expect_identical(nrow(series$mortality), 151L)
  expect_gte(nrow(series$mortality_monthly), year_month_to_idx(2100, 1))
  expect_identical(names(series$tx_default)[c(1, length(series$tx_default))], c("1993", "2100"))
  expect_error(param_init(P, "TX", prg_chng = def_prgchng(P), ttt_list = def_ttt()),
               "param_init\\(\\) for 'TX' needs model_load\\(\"TX\"\\); 'CA' is loaded")
})

test_that("param_init needs model_load first", {
  had <- exists("input_series", envir = .mitus_state, inherits = FALSE)
  saved <- if (had) get("input_series", envir = .mitus_state)
  if (had) rm("input_series", envir = .mitus_state)
  withr::defer(if (had) assign("input_series", saved, envir = .mitus_state))
  expect_error(input_series(), "Run model_load\\(\\) before building model parameters")
})
