#' Aging denominators by age group
#'
#' For each age group but the oldest, the inverse of the share of its
#' population in its oldest single year of age, smoothed to months and held at
#' the end of the smoothed curve through \code{month}. The US uses its yearly
#' national population counts. A state or county uses its census population
#' counts, which the smoothing spaces ten years apart, and sets missing or zero
#' shares to 0.0198.
#'
#' State census counts give the oldest ages as one open-ended group, such as 90
#' and over, recorded at a single age that differs between census years. They
#' cannot give the share of the 85-94 group at age 94, so a state or county
#' takes that share from the national counts in each of its census years.
#'
#' @name age_denom
#' @param loc location code
#' @param month number of months
#' @return matrix of monthly aging denominators, one row per month and one
#'   column per age group but the oldest
#' @export
age_denom <- function(loc, month) {
  national <- national_by_age("US/US_PopCountsByAge.rds")
  if (loc == "US") {
    ltd <- oldest_age_share(national)
    smooth <- SmoCurve
  } else {
    popdist <- location_census_population(loc)
    if (!all(colnames(popdist) %in% colnames(national))) {
      stop("The census population for ", loc, " has years without national population counts",
           call. = FALSE)
    }
    ltd <- oldest_age_share(popdist)
    ltd[10, ] <- oldest_age_share(national[, colnames(popdist)])[10, ]
    ltd[is.na(ltd)] <- 0
    ltd[ltd == 0] <- 0.01980337
    smooth <- SmoCurve_decade_month
  }
  ltd <- 1 / ltd
  td <- matrix(NA, 10, month)
  for (i in 1:10) {
    smoothed <- smooth(as.numeric(ltd[i, ]))
    td[i, seq_along(smoothed)] <- smoothed
    td[i, (length(smoothed) + 1):month] <- smoothed[length(smoothed)]
  }
  t(td) * 12
}

#' Share of each age group but the oldest in its oldest single year of age
#'
#' @param popdist population counts, one row per single year of age from 0 and
#'   one column per year
#' @return matrix with one row per age group but the oldest and one column per
#'   year
#' @noRd
oldest_age_share <- function(popdist) {
  t(sapply(age_group_rows[1:10], function(rows) {
    popdist[max(rows), ] / colSums(popdist[rows, ])
  }))
}

#' Census population by single year of age for a state or county
#'
#' Read from the location's own file if it has one, otherwise from the
#' all-state file. For AK and HI the first census year is estimated by a linear
#' trend over the later census years.
#'
#' @return matrix with one row per single year of age and one column per census
#'   year
#' @noRd
location_census_population <- function(loc) {
  loc_file <- find_loc_file(loc, "PopCountsByAge", data_dir = get_data_dir(), required = FALSE)
  if (!is.null(loc_file)) {
    pop <- readRDS(loc_file)
  } else {
    state_id <- as.data.frame(get_stateID())
    fips <- as.numeric(as.matrix(state_id[which(state_id$USPS == loc), 2]))
    pop <- readRDS(system.file("ST/ST_PopCountsByAge.rds", package = "MITUS"))
    pop <- pop[pop$fips == fips, ]
  }
  popdist <- pop[, grep("^[0-9]{4}$", colnames(pop))]
  rownames(popdist) <- pop[, 2]
  popdist <- as.matrix(popdist)
  if (loc %in% c("AK", "HI")) {
    later_years <- rev(as.integer(colnames(popdist)[-1]))
    for (i in seq_len(nrow(popdist))) {
      trend <- data.frame(count = rev(popdist[i, -1]), year = later_years)
      coeff <- stats::lm(count ~ year, data = trend)$coefficients
      estimate <- round((coeff[2] * as.integer(colnames(popdist)[1])) + coeff[1])
      popdist[i, 1] <- max(estimate, 0)
    }
  }
  popdist
}
