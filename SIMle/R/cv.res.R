
# cv results 
source(paste(getwd(), "/R/SIMle.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.db1-20.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.plot.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.original_code.v1.R", sep = ""))




#' Cross validation result by specific criteria
#' @description this function gets the cross validation result by specific criteria.
#'
#' @param ts ts is the data set which is a time series data typically
#' @param c the maximum value of number of basis for time input
#' @param d the maximum value of number of basis for variate input
#' @param b_time type of basis for time input 
#' @param b_timese type of basis for variate input
#' @param mp_type select type of mapping function, "algeb" indicates algebraic mapping on the real line. "logari" represents logarithmic mapping on the real line
#' @param ops Criteria for choosing the number of bases are provided by the package, offering four options: "AIC," "BIC," "CV," and "Kfold," each corresponding to a specific Criteria 
#' @param r indicates number of variate 
#' @param s s is a positive scaling factor, the default is 1
#' @param per the percentage for test set used in "CV" option 
#' @param k the number of fold used in "Kfold" option
#' 
#' @return A data frame containing the criterion values corresponding to "c" and "d". The first element refers to the optimal number of basis for time input, and the second 
#'         element refers to the optimal number of basis for variate. 
#' @export

cv.res <- function(ts,c, d, b_time, b_timese, mp_type, ops, r = 1, s = 1, per = 0, k = 0){
  if(ops == "AIC"){
    ops_cd = cd_selection(ts, c, d, b_time, b_timese, mp_type, ops, r = r, s = s)
    ops_cd = matrix(as.numeric(ops_cd),ncol = ncol(ops_cd))
    return(ops_cd)
    
  } else if(ops == "BIC"){
    ops_cd = cd_selection(ts, c, d, b_time, b_timese, mp_type, ops, r = r, s = s)
    ops_cd = matrix(as.numeric(ops_cd),ncol = ncol(ops_cd))
    return(ops_cd)
    
  } else if(ops == "CV"){
    ops_cd = cd_selection(ts, c, d, b_time, b_timese, mp_type, ops, per = per, r = r, s = s)
    ops_cd = matrix(as.numeric(ops_cd),ncol = ncol(ops_cd))
    return(ops_cd)
    
  } else if(ops == "Kfold"){
    ops_cd = cd_selection(ts, c, d, b_time, b_timese, mp_type, ops, k = k, r = r, s = s)
    ops_cd = matrix(as.numeric(ops_cd),ncol = ncol(ops_cd))
    return(ops_cd)
  } else{
    return(stop("Invalid option!"))
  }
}


