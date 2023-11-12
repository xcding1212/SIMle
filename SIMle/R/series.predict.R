

source(paste(getwd(), "/R/SIMle.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.db1-20.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.original_code.v1.R", sep = ""))

# source(paste(getwd(), "/R/SIMle.Legen.v1.R", sep = ""))

# for validation, t != 1 


# prediction step, what exact formula used in this step? 

#' Predicting time series with 1 step
#' 
#' @description This function predicts the time series data basis on the estimation.
#' 
#' @param ts ts is the data set which is a time series data typically
#' @param c number of basis for time input
#' @param d number of basis for variate input 
#' @param b_time type of basis for time input 
#' @param b_timese type of basis for variate input
#' @param mp_type select type of mapping function, "algeb" indicates algebraic mapping on the real line. "logari" represents logarithmic mapping on the real line
#' @param r indicates number of variate
#' @param s  s is a positive scaling factor, the default is 1
#' @param n_esti number of points for estimation, the default is 2000 
#'
#' @return predictive values for time series 
#' @export
#' 

series.predict = function(ts, c, d, b_time, b_timese, mp_type, r=1, s=1, n_esti = 2000){
  res = esti_ts(ts, c, d, b_time, b_timese, mp_type, r=r, s=s, n_esti = n_esti)
  return(res[length(res)])
}