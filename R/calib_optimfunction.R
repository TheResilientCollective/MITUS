#'This function is used to run a optimization on the input parameters.
#'The only input to the function is a "batch number" which will determine
#'which row of the starting values dataframe the optimization will use.
#'It will return 8 datasets of the optimized parameters --one from each
#'optimization step. The final (8th) optimization step is an univariate
#'optimization.
#'The function requires the use of the ParamInit Rdata file and the
#'StartValues Rdata file.
#'@name optim_b
#'@param df dataframe or matrix of starting values data frame
#'@param samp_i which rows of the data frame to use
#'@param TB boolean for TB likelihoods
#'@param n_cores parallelization
#'@return 8 datasets from optimization loop
#'@export

optim_b <- function(df, samp_i=1,n_cores=2, TB=1){
# data("StartVal_2018-08-06", package = "MITUS")


posterior = function(theta) {
   -lprior(theta) - llikelihood(theta,n_cores, TB=TB)
 }

if(min(dim(as.data.frame(df)))==1) {
  df1 <- as.numeric(df)
  names(df1) <- names(df)
} else{
  df1 <- as.numeric(df[samp_i,])
}

b<-samp_i
# for (i in min(b, nrow(StartVal))){
  o1  <- optim(df1, posterior, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o1$value
  save(o1,file=paste("Opt_US_r1_", b,"_", Sys.Date(),".rda",sep=""))
  o2  <- optim(o1$par,posterior, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o2$value
  save(o2,file=paste("Opt_US_r2_", b,"_", Sys.Date(),".rda",sep=""))
  o3  <- optim(o2$par, posterior, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o3$value
  save(o3,file=paste("Opt_US_r3_", b,"_", Sys.Date(),".rda",sep=""))
  o4  <- optim(o3$par ,posterior, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o4$value
  save(o4,file=paste("Opt_US_r4_", b,"_", Sys.Date(),".rda",sep=""))
  o5  <- optim(o4$par, posterior, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o5$value
  save(o5,file=paste("Opt_US_r5_", b,"_", Sys.Date(),".rda",sep=""))
  o6  <- optim(o5$par ,posterior, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o6$value
  save(o6,file=paste("Opt_US_r6_", b,"_", Sys.Date(),".rda",sep=""))
  o7  <- optim(o6$par, posterior, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o7$value
  save(o7,file=paste("Opt_US_r7_", b,"_", Sys.Date(),".rda",sep=""))
  # o8 <- UnivOptim(o7$par) ; o8$value
  # save(o8,file=paste("Opt_US_r8_", b,"_", Sys.Date(),".rda",sep=""))
# }
}

#' Optimize sub-national parameters from one starting point
#'
#' Runs nine rounds of \code{optim} on the negative log posterior, alternating
#' BFGS and Nelder-Mead, each round starting from the previous round's result.
#' Requires the globals set by \code{model_load}.
#'
#' @name optim_b_st
#' @param df starting values on the standard normal scale, as a named vector or
#'   as a matrix with one row per starting point
#' @param samp_i row of \code{df} to start from, also used in the saved file
#'   names
#' @param n_cores number of cores passed to \code{llikelihood_st}
#' @param loc location code: a state USPS code or a county code loaded by
#'   \code{model_load}
#' @param TB 1 to include the TB likelihood terms, 0 for demographic terms only
#' @param calib_end_year last simulated year; NULL simulates through the last
#'   year of the importance weights set by \code{model_load}
#' @return NULL, invisibly. Each round's \code{optim} result is saved to the
#'   working directory as \code{Opt_<loc>_r<round>_<samp_i>_<date>.rda}.
#' @export


optim_b_st <- function(df,
                       samp_i = 1,
                       n_cores = 2,
                       loc = "MA",
                       TB = 1,
                       calib_end_year = NULL) {
  calib_end_year <- resolve_calib_end_year(calib_end_year)
  posterior_st <- function(theta) {
    -lprior(theta) - llikelihood_st(theta, loc, n_cores, TB = TB,
                                    calib_end_year = calib_end_year)
  }

  if(min(dim(as.data.frame(df)))==1) {
    df1 <- as.numeric(df)
    names(df1) <- names(df)
  } else{
    df1 <- as.numeric(df[samp_i,])
  }

  b<-samp_i
  o1  <- optim(df1, posterior_st, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o1$value
  save(o1,file=paste("Opt_",loc,"_r1_", b,"_", Sys.Date(),".rda",sep=""))
  o2  <- optim(o1$par,posterior_st, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o2$value
  save(o2,file=paste("Opt_",loc,"_r2_", b,"_", Sys.Date(),".rda",sep=""))
  o3  <- optim(o2$par, posterior_st, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o3$value
  save(o3,file=paste("Opt_",loc,"_r3_", b,"_", Sys.Date(),".rda",sep=""))
  o4  <- optim(o3$par ,posterior_st, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o4$value
  save(o4,file=paste("Opt_",loc,"_r4_", b,"_", Sys.Date(),".rda",sep=""))
  o5  <- optim(o4$par, posterior_st, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o5$value
  save(o5,file=paste("Opt_",loc,"_r5_", b,"_", Sys.Date(),".rda",sep=""))
  o6  <- optim(o5$par ,posterior_st, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o6$value
  save(o6,file=paste("Opt_",loc,"_r6_", b,"_", Sys.Date(),".rda",sep=""))
  o7  <- optim(o6$par, posterior_st, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o7$value
  save(o7,file=paste("Opt_",loc,"_r7_", b,"_", Sys.Date(),".rda",sep=""))
  o8  <- optim(o7$par ,posterior_st, method ="Nelder-Mead", control=list(maxit=1000,trace=1,reltol=sqrt(.Machine$double.eps)/5));  o8$value
  save(o8,file=paste("Opt_",loc,"_r8_", b,"_", Sys.Date(),".rda",sep=""))
  o9  <- optim(o8$par, posterior_st, method ="BFGS", control=list(maxit=400,trace=5,reltol=sqrt(.Machine$double.eps)/5)) ; o9$value
  save(o9,file=paste("Opt_",loc,"_r9_", b,"_", Sys.Date(),".rda",sep=""))
}
