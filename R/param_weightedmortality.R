# Single years of age, 0 to 110, in each of the model's 11 age groups.
age_group_rows <- list(1:5, 6:15, 16:25, 26:35, 36:45, 46:55, 56:65, 66:75, 76:85, 86:95, 96:111)

#' Background mortality rates by age group for a location
#'
#' For the US, the national death counts over the national population in each
#' age group. For a state or county, the location's life table rates averaged
#' over the single years of age in each age group, weighted by the national
#' population at each age, for the years the life table covers.
#'
#' @name weight_mort
#' @param loc location code
#' @return matrix of yearly mortality rates, one row per year, named by year,
#'   and one column per age group
#' @export
weight_mort <- function(loc) {
  popdist <- national_by_age("US/US_PopCountsByAge.rds")
  if (loc == "US") {
    death_age <- national_by_age("US/US_MortalityCountsByAge.rds")
    years <- intersect(colnames(death_age), colnames(popdist))
    return(sapply(age_group_rows, function(rows) {
      colSums(death_age[rows, years]) / colSums(popdist[rows, years])
    }))
  }
  life_table <- as.matrix(location_life_table(loc))
  life_years <- colnames(life_table)
  if (!all(life_years %in% colnames(popdist))) {
    stop("The life table for ", loc, " covers years without national population data",
         call. = FALSE)
  }
  popdist <- popdist[, life_years]
  sapply(age_group_rows, function(rows) {
    colSums(life_table[rows, ] * popdist[rows, ]) / colSums(popdist[rows, ])
  })
}

#' Mortality rates by single year of age and year for a state or county
#'
#' Read from the location's own life table file if it has one, otherwise from
#' the all-state life tables.
#'
#' @return data frame with one column per year
#' @noRd
location_life_table <- function(loc) {
  loc_file <- find_loc_file(loc, "SingleYearLifeTable", data_dir = get_data_dir(),
                            required = FALSE)
  if (!is.null(loc_file)) {
    return(readRDS(loc_file)[, -1])
  }
  index <- which(as.data.frame(get_stateID())$USPS == loc)
  readRDS(system.file("ST/ST_SingleYearLifeTable.rds", package = "MITUS"))[[index]][, -1]
}

#' National counts by single year of age, one column per year, from a package
#' file with an age column followed by yearly columns
#'
#' Year columns are read by their leading four digits: the national population
#' counts label 1959 as "1959+".
#' @noRd
national_by_age <- function(file) {
  table <- readRDS(system.file(file, package = "MITUS"))
  year_cols <- grep("^[0-9]{4}", colnames(table))
  years <- substr(colnames(table)[year_cols], 1, 4)
  if (any(diff(as.integer(years)) != 1)) {
    stop(file, " must have one column per year", call. = FALSE)
  }
  counts <- as.matrix(table[, year_cols])
  dimnames(counts) <- list(table[, 1], years)
  counts
}
