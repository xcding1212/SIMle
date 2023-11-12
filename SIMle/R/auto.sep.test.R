source(paste(getwd(), "/R/SIMle.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.db1-20.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.plot.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.original_code.v1.R", sep = ""))
source(paste(getwd(), "/R/sep.test.R", sep = ""))


#' Automated separability test 
#' @description This function utilizes Simultaneous Confidence Regions (SCR) for the automated execution of separability tests with with chosen bases.
#'
#' @param ts ts is the data set which is a time series data typically
#' @param c the maximum value of number of basis for time input
#' @param d the maximum value of number of basis for variate input
#' @param b_time type of basis for time input 
#' @param b_timese type of basis for variate input
#' @param mp_type select type of mapping function, "algeb" indicates algebraic mapping on the real line. "logari" represents logarithmic mapping on the real line
#' @param type select type of estimation."fixt" indicates fix time t. "fixx" represents fix variate. The default is "not" which indicates no fix value 
#' @param ops Criteria for choosing the number of bases are provided by the package, offering four options: "AIC," "BIC," "CV," and "Kfold," each corresponding to a specific Criteria
#' @param m the window size for the simultaneous confidence region procedure, with the default being 'MV,' which stands for the Minimum Volatility method
#' @param option Choose either the fixed value separability test or the for all values separability test. The default is "each"
#' @param fix_num fix_num indicates the use of fixed-value nonlinear time series regression. If "fixt" is chosen, it represents a fixed time value. Otherwise, 
#'        if not selected, it pertains to a fixed variate value
#' @param r indicates number of variate 
#' @param s s is a positive scaling factor, the default is 1
#' @param per the percentage for test set used in "CV" option 
#' @param k the number of fold used in "Kfold" option
#' @param upper upper The upper bound for the variate basis domain. The default value is 10. When "algeb" or "logari" is chosen, the domain is automatically set from -upper to upper 
#'        
#' 
#' @return If option equals to "each". A list containing dataframes with three columns each. The first column represents input values. The second column contains values of the estimated function, 
#'         along with their upper and lower bounds, which are used for separability testing. The third column is a factor indicating the types corresponding to the values 
#'         in the second column. If option equals to "all". It returns a list containing 2 lists. The first one contains SCR dataframe  and the second one is estimated function by using separating 
#'         function of t and function of x.
#' @export



auto.sep.test <- function(ts, c, d, b_time, b_timese, mp_type, type = "not", ops, m = "MV", option ="each", fix_num = 0, r = 1, s = 1, per = 0, k = 0, upper = 10){
  if(option == "each"){
    best_cd = best_cd.auto.fit(ts, c, d, b_time, b_timese, mp_type, ops, r = r, s = s, per = per, k = k)
    res = sep.test(ts, as.numeric(best_cd[1]), as.numeric(best_cd[2]), m, b_time, b_timese, mp_type, type,  fix_num = fix_num, option = option, r = r, s = s, upper = upper)
    return(res)
    
  } else if(option == "all"){
    res = sepera_3d(ts, as.numeric(best_cd[1]), as.numeric(best_cd[2]), m, b_time, b_timese, mp_type, r = r, s = s,  upper = upper)
    return(res)
  } else{
    return(stop("Invalid option!"))
  }
}
