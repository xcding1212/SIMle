
source(paste(getwd(), "/R/SIMle.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.db1-20.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.plot.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.original_code.v1.R", sep = ""))




#' Exact form test 
#' @description This function employs the L2 test for the user-specific execution of exact form tests.
#'
#' @param ts ts is the data set which is a time series data typically
#' @param c number of basis for time input
#' @param d number of basis for variate input
#' @param m the window size for the simultaneous confidence region procedure, with the default being 'MV,' which stands for the Minimum Volatility method 
#' @param b_time type of basis for time input 
#' @param b_timese type of basis for variate input
#' @param mp_type select type of mapping function, "algeb" indicates algebraic mapping on the real line. "logari" represents logarithmic mapping on the real line
#' @param exact_func A list contains elements that are matrix contain exact functions, which are desired to be tested. The k-th element represents the k-th variable. 
#'        The matrix contains values of the exact function within its domain  
#' @param r indicates number of variate  
#' @param s s is a positive scaling factor, the default is 1
#' @param upper The upper bound for the variate basis domain. The default value is 10. When "algeb" or "logari" is chosen, the domain is automatically set from -upper to upper. 
#'        
#' 
#' @return A list whose elements are p value of exact form test. Each element in the list represents p-values in the order of variates.
#' @export




exact.test <- function(ts, c, d, m = "MV", b_time, b_timese, mp_type, exact_func, r = 1, s = 1, upper = 10){
  basis_candi = c("Legen","Cheby","tri", "cos", "sin", "Cspli", "db1", "db2", "db3", "db4", "db5",
                  "db6", "db7", "db8", "db9", "db10",
                  "db11", "db12", "db13", "db14", "db15",
                  "db16", "db17", "db18", "db19", "db20",
                  "cf1", "cf2", "cf3", "cf4", "cf5"
  )
  
  #if(m == "MV"){ # refer to minimum volatility (MV) method
  #  m =  16 # min_vola()
  #} 
  
  if((b_time %in% basis_candi) && (b_timese %in% basis_candi)){
      res = exact_form_test(ts, c, d, m, b_time, b_timese, mp_type, exact_func, r = r, s = s, upper = upper)
      return(res)
    } else{
    return(stop("Invalid option!"))
  }
}
