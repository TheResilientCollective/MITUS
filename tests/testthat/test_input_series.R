test_that("mortality is forecast by each age group's yearly decline", {
  mortality <- matrix(1, 151, 11)
  mortality[year_to_idx(2016), ] <- 2
  forecast <- forecast_mortality(mortality, 2017)
  expect_identical(forecast[year_to_idx(2016), ], rep(2, 11))
  expect_equal(forecast[year_to_idx(2017), ], 2 * (1 - mortality_declines))
  expect_equal(forecast[year_to_idx(2018), c(1, 3, 11)], 2 * (1 - c(0.0159, 0.0090, 0.0069))^2)
  expect_identical(forecast[year_to_idx(2015), ], rep(1, 11))
})

test_that("yearly inputs are read by their years, or from 1950 without them", {
  expect_identical(annual_values(c(`1950` = 1, `1951` = 2)), c(`1950` = 1, `1951` = 2))
  expect_identical(annual_values(c(5, 6, 7)), c(`1950` = 5, `1951` = 6, `1952` = 7))
  expect_error(annual_values(c(`1951` = 1)), "must run year by year from 1950")
  expect_error(annual_values(c(`1950` = 1, `1952` = 2)), "must run year by year from 1950")
})

test_that("an input is used through a given year, or all of it", {
  x <- annual_values(rep(1, 151))
  expect_identical(names(through_year(x, 2019))[c(1, 70)], c("1950", "2019"))
  expect_length(through_year(x, 2019), 70)
  expect_identical(through_year(x, NULL), x)
})

test_that("immigration grows and births hold after their data", {
  immigration <- forecast_immigration(annual_values(c(rep(1, 73), 2)))
  expect_identical(names(immigration), as.character(1950:2100))
  expect_identical(immigration[["2023"]], 2)
  expect_equal(immigration[c("2024", "2033")], c(`2024` = 2 * 1.005, `2033` = 2 * 1.005^10))
  births <- forecast_births(annual_values(c(rep(1, 73), 3)))
  expect_identical(unname(births[c("2022", "2023", "2024", "2100")]), c(1, 3, 3, 3))
  full <- annual_values(seq_len(151))
  expect_identical(forecast_immigration(full), full)
  expect_identical(forecast_births(full), full)
})

test_that("the immigrant age distribution holds its last year after its data", {
  dist <- immigration_age_dist(matrix(c(rep(1, 11 * 73), 2:12), 11, 74))
  expect_identical(dim(dist), c(11L, 151L))
  expect_identical(dist[, 151], 2:12 + 0)
  full <- matrix(seq_len(11 * 151), 11, 151, dimnames = list(NULL, 1950:2100))
  expect_identical(immigration_age_dist(full[, 151:1]), full)
})

test_that("national counts are read by year, including a labelled 1959+", {
  counts <- national_by_age("US/US_PopCountsByAge.rds")
  expect_identical(colnames(counts)[c(1, 10, ncol(counts))], c("1950", "1959", "2018"))
  deaths <- national_by_age("US/US_MortalityCountsByAge.rds")
  expect_identical(colnames(deaths)[ncol(deaths)], "2017")
})

test_that("a life table with the national rates gives the national mortality", {
  deaths <- national_by_age("US/US_MortalityCountsByAge.rds")
  pop <- national_by_age("US/US_PopCountsByAge.rds")
  years <- c("1990", "2000", "2016")
  local_mocked_bindings(location_life_table = function(loc) {
    as.data.frame(deaths[, years] / pop[, years])
  })
  expect_equal(weight_mort("XX"), weight_mort("US")[years, ])
})

test_that("each age group's mortality comes only from its own ages", {
  ages <- 0:110
  local_mocked_bindings(location_life_table = function(loc) {
    data.frame(`2000` = as.numeric(ages >= 5 & ages <= 14),
               `2016` = as.numeric(ages >= 5 & ages <= 14), check.names = FALSE)
  })
  expected <- matrix(0, 2, 11, dimnames = list(c("2000", "2016"), NULL))
  expected[, 2] <- 1
  expect_equal(weight_mort("XX"), expected)
})

test_that("treatment inputs hold their last value, named by year", {
  held <- held_to_end(cbind(1993:1995, c(0.1, 0.2, 0.3)))
  expect_identical(names(held), as.character(1993:2100))
  expect_identical(unname(held[c("1995", "2100")]), c(0.3, 0.3))
})

test_that("the immigrant TB trend declines after its data and doubles for ND", {
  burden <- annual_values(seq(2, 1, length.out = 74))
  trend <- immigrant_tb_trend(burden, "CA")
  expect_identical(names(trend), as.character(1950:2100))
  expect_equal(trend[["2023"]], 90 / 1e5)
  expect_equal(trend[["2024"]], trend[["2023"]] * (1 - immigrant_tb_decline))
  nd <- immigrant_tb_trend(burden, "ND")
  expect_equal(nd[["2020"]], 2 * nd[["2010"]])
  expect_identical(nd[["2009"]], trend[["2009"]])
})

test_that("model_load builds the location's input series", {
  local_model_load("CA")
  series <- input_series()
  expect_identical(series$loc, "CA")
  expect_identical(series$input_through,
                   c(mortality = 2016, immigration = 2019, immigrant_tb_trend = 2018L))
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

test_that("the pandemic immigration reduction applies only to forecast immigration", {
  local_model_load("CA")
  run <- function(immig) {
    p <- param_init(P, "CA", prg_chng = def_prgchng(P), ttt_list = def_ttt(), immig = immig)
    p$ImmNon
  }
  month <- year_month_to_idx(2021, 1)
  expect_lt(run(0.5)[month, 5], run(99)[month, 5])
  series <- input_series()
  series$input_through[["immigration"]] <- 2023
  assign("input_series", series, envir = .mitus_state)
  expect_identical(run(0.5), run(99))
})

test_that("mortality through 2100 is used as given beyond the life table years", {
  local_model_load("CA")
  file_rates <- as.matrix(Inputs$BgMort)[, -1]
  given <- location_mortality("CA", Inputs$BgMort)
  expect_identical(attr(given, "input_through"), 2100L)
  later <- year_to_idx(2017):151
  expect_identical(unname(given[later, ]), unname(file_rates[later, ]))
  legacy <- location_mortality("CA", Inputs$BgMort, through = 2016)
  expect_identical(attr(legacy, "input_through"), 2016)
  expect_identical(legacy[year_to_idx(2016), ], given[year_to_idx(2016), ])
  expect_equal(legacy[year_to_idx(2017), ], given[year_to_idx(2016), ] * (1 - mortality_declines))
})

test_that("a location's own model input file is used as given", {
  local_model_load("CA")
  own <- build_input_series("CA", Inputs, legacy = FALSE)
  expect_identical(own$input_through,
                   c(mortality = 2100L, immigration = 2100L, immigrant_tb_trend = 2018L))
  expect_identical(own$immigration, annual_values(Inputs$ImmigInputs$TotByYear))
})

test_that("the package's own data files are recognised", {
  expect_true(is_package_file(system.file("ST", "ST_PopCountsByAge.rds", package = "MITUS")))
  expect_false(is_package_file(withr::local_tempfile(fileext = ".rds", lines = "")))
})
