source(paste(getwd(), "/R/SIMle.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.db1-20.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.plot.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.original_code.v1.R", sep = ""))




#' Automated creation of a Simultaneous Confidence Region (SCR) for the estimated function
#' @description This function generates a Simultaneous Confidence Region (SCR) for the estimated function with chosen bases.
#' @param ts ts is the data set which is a time series data typically
#' @param c the maximum value of number of basis for time input 
#' @param d the maximum value of number of basis for variate input
#' @param b_time type of basis for time input 
#' @param b_timese type of basis for variate input
#' @param mp_type select type of mapping function, "algeb" indicates algebraic mapping on the real line. "logari" represents logarithmic mapping on the real line
#' @param type select type of estimation."fixt" indicates fix time t estimation. "fixx" represents fixed variate estimation
#' @param ops Criteria for choosing the number of bases are provided by the package, offering four options: "AIC," "BIC," "CV," and "Kfold," each corresponding to a specific Criteria
#' @param m the window size for the simultaneous confidence region procedure, with the default being 'MV,' which stands for the Minimum Volatility method
#' @param fix_num fix_num indicates the use of fixed-value nonlinear time series regression. If "fixt" is chosen, it represents a fixed time value. Otherwise, 
#'        if not selected, it pertains to a fixed variate value
#' @param r indicates number of variate 
#' @param s s is a positive scaling factor, the default is 1
#' @param per the percentage for test set used in "CV" option 
#' @param k the number of fold used in "Kfold" option
#' @param upper upper The upper bound for the variate basis domain. The default value is 10. When "algeb" or "logari" is chosen, the domain is automatically set from -upper to upper
#'        
#' 
#' @return A list containing dataframes with three columns each. The first column corresponds to input values. The second column contains values of the estimated function, along with 
#'         their upper and lower bounds. The third column is a factor that indicates the types associated with the values in the second column.
#' @export




auto.SCR <- function(ts, c, d, b_time, b_timese, mp_type, type, ops, m = "MV", fix_num = 0, r = 1, s = 1, per = 0, k = 0, upper = 10){ # ops: criteria, type: fixx or fixt
  best_cd = best_cd.auto.fit(ts, c, d, b_time, b_timese, mp_type, ops, r = r, s = s, per = per, k = k)
  res = fix.SCR(ts, as.numeric(best_cd[1]), as.numeric(best_cd[2]), m, b_time, b_timese, mp_type, type, fix_num = fix_num, r = r, s = s, upper = upper)
  return(res)
}
