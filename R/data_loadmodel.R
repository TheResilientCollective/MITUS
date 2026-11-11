#' Load calibration data and inputs for a location
#'
#' Reads the location's model inputs and parameter files into the global
#' environment, builds the location's input series for \code{param_init}, and
#' sets the importance weights \code{wts}. For the US it also reads the
#' calibration data its likelihood uses. State and county calibration reads one
#' location's targets with \code{calib_targets}.
#'
#' Each input is read by its years and taken as observed through its last
#' year. The years after it, through 2100, are forecast by fixed rules, all set
#' at the top of \code{R/input_series.R}:
#' \itemize{
#'   \item background mortality declines by a fixed yearly rate for each age group
#'   \item immigration volume grows 0.5\% a year
#'   \item the immigrant age distribution holds its last observed year
#'   \item births hold their last observed value
#'   \item the immigrant TB prevalence trend declines 1.5\% a year
#'   \item treatment default and quality hold their last values
#'   \item the aging denominators hold the end of their smoothed curve
#' }
#' An input that runs through 2100 holds its own forecast after its data, so it
#' is taken as observed through the year the previous code assumed: 2016 for
#' background mortality and 2019 for immigration. The 2020 pandemic immigration
#' reduction applies only where immigration is observed through a year before
#' 2020.
#'
#' @name model_load
#' @param loc location code: \code{"US"}, a state USPS code, or the code of
#'   another location whose files are in \code{data_dir}
#' @param data_dir optional external data directory to search before package
#'   inst/
#' @param last_weight_year last calendar year covered by the calibration
#'   importance weights \code{wts}. Likelihoods that take their year window from
#'   the calibration data use only the data years present in \code{names(wts)},
#'   so this sets the most recent year of data those targets can use.
#'   Calibration simulates through this year unless a later
#'   \code{calib_end_year} is given.
#' @return NULL, invisibly
model_load <- function(loc = "US", data_dir = NULL, last_weight_year = 2019) {
  library(mnormt)
  library(parallel)
  library(lhs)
  library(Rcpp)
  library(MCMCpack)
  library(MASS)

  # Store data_dir, or NULL, so a directory from an earlier call is not reused.
  assign("data_dir", data_dir, envir = .mitus_state)

# load necessary datasets
if (loc=="US"){
  CalibDat<<-readRDS(find_loc_file("US","CalibDat",data_dir=data_dir))
  CalibDatCases<<-CalibDat
  ## ParamInit and StartVal were replaced on 8/1 to clean up, but 07-07 are last calibrated version
  ParamInit<<-as.data.frame(readRDS(find_loc_file("US","ParamInit",data_dir=data_dir)))
  StartVal<<-readRDS(find_loc_file("US","StartVal",data_dir=data_dir))
  Inputs<<-readRDS(find_loc_file("US","Inputs",data_dir=data_dir))
  Opt <<- readRDS(find_loc_file("US","Optim",data_dir=data_dir,required=FALSE))
  Par <<- readRDS(find_loc_file("US","Param",data_dir=data_dir,required=FALSE))
} else {
  # State and county calibration reads one location's targets with
  # calib_targets(), so any US target globals left from an earlier load are
  # removed.
  target_globals <- c("CalibDat", "CalibDatCases")
  rm(list = intersect(target_globals, ls(globalenv())), envir = globalenv())
  ParamInit_st<<-ParamInit<<-readRDS(find_loc_file("ST","ParamInit",fallback_loc=loc,data_dir=data_dir))
  StartVal_st<<-StartVal<<-readRDS(find_loc_file("ST","StartVal",fallback_loc=loc,data_dir=data_dir))
  Inputs<<-readRDS(find_loc_file(loc,"ModelInputs",data_dir=data_dir))
  par_file<-find_loc_file(loc,"Param",data_dir=data_dir,required=FALSE)
  if(!is.null(par_file)){
    Par <<- readRDS(par_file)
  }
  #last input change was to update the RR active TB by age in immigrants
}
if (loc=="US"){
  wts <<- make_impt_weights(last_year=last_weight_year)
  P  <<- ParamInit[,1]
  names(P) <<- rownames(ParamInit)

  ii <<-  ParamInit[,5]==1
  ParamInitZ <<- ParamInit[ParamInit$Calib==1,]
  idZ0 <<- ParamInitZ[,4]==0
  idZ1 <<- ParamInitZ[,4]==1
  idZ2 <<- ParamInitZ[,4]==2
} else {
  wts <<- make_impt_weights(last_year=last_weight_year)

  #creation of background parameters
  #elements of P will be replaced from either the StartVals in the case
  #of optimization or the user inputted dataset

  P  <<- ParamInit_st[,1]
  names(P) <<- rownames(ParamInit_st)
  ii <<-  ParamInit_st[,5]==1
  ParamInitZ <<- ParamInit_st[ParamInit_st$Calib==1,]
  idZ0 <<- ParamInitZ[,4]==0
  idZ1 <<- ParamInitZ[,4]==1
  idZ2 <<- ParamInitZ[,4]==2
  ParamInit<<-ParamInit_st
}
  prgchng<<-def_prgchng(P)
  assign("input_series", build_input_series(loc, Inputs), envir = .mitus_state)

return(invisible(NULL))
}

#'loads data for the geography of interest
#'@name model_load_demo
#'@param loc \code{"US"}, the only location the demography calibration covers
#'@return void
model_load_demo<-function(loc="US"){
  library(mnormt)
  library(parallel)
  library(lhs)
  library(Rcpp)
  library(MCMCpack)
  library(MASS)

  if (loc != "US") {
    stop("model_load_demo() covers the US only", call. = FALSE)
  }

  #'lazy load necessary datasets
  #'Model Input
  CalibDat <<- readRDS(find_loc_file("US", "CalibDat"))
  ParamInit_demo <<- readRDS(find_loc_file("US", "ParamInitdemo"))
  StartVal_demo <<- readRDS(find_loc_file("US", "StartValdemo"))
  Inputs <<- readRDS(find_loc_file("US", "Inputs"))
  wts <<- CalibDat[["ImptWeights"]]
  #'creation of background parameters
  #'elements of P will be replaced from either the StartVals in the case
  #'of optimization or the user inputted dataset
  #'
  P  <<- ParamInit_demo[,1]
  names(P) <<- rownames(ParamInit_demo)

  ii <<-  ParamInit_demo[,5]==1
  ParamInitZ <<- ParamInit_demo[ParamInit_demo$Calib==1,]
  idZ0 <<- ParamInitZ[,4]==0
  idZ1 <<- ParamInitZ[,4]==1
  idZ2 <<- ParamInitZ[,4]==2

  #'format the calibration data
  # targets[["ParamInitZ"]]<-ParamInitZ
  # targets[["idZ0"]]<-idZ0
  # targets[["idZ1"]]<-idZ1
  # targets[["idZ2"]]<-idZ2
  # targets[["P"]]<-P
  ## ALSO LOAD IN THE BASELINE PROGRAM CHANGE VECTOR
  prgchng<<-def_prgchng(P)
  return(invisible(NULL))
}
