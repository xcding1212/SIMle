source(paste(getwd(), "/R/SIMle.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.db1-20.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.plot.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.original_code.v1.R", sep = ""))




#' User-specified time-homogeneity test
#' @description This function utilizes Simultaneous Confidence Regions (SCR) for the automated execution of time-homogeneity tests 
#' @param ts ts is the data set which is a time series data typically
#' @param c number of basis for time input
#' @param d number of basis for variate input
#' @param m the window size for the simultaneous confidence region procedure, with the default being 'MV,' which stands for the Minimum Volatility method 
#' @param b_time type of basis for time input 
#' @param b_timese type of basis for variate input
#' @param mp_type select type of mapping function, "algeb" indicates algebraic mapping on the real line. "logari" represents logarithmic mapping on the real line
#' @param fix_num fix_num indicates fixed value for time
#' @param r indicates number of variate 
#' @param s s is a positive scaling factor, the default is 1
#' @param n_point number of points for SCR, the default is 2000
#' @param upper upper The upper bound for the variate basis domain. The default value is 10. When "algeb" or "logari" is chosen, the domain is automatically set from -upper to upper 
#'        
#' 
#' @return A list is returned, containing dataframes with three columns each. The first column pertains to input values, the second column contains values of the estimated function 
#'         along with their upper and lower bounds, which are used for time-homogeneity testing. The third column serves as a factor indicating the types corresponding to the values in 
#'         the second column.
#' @export


homo.test <- function(ts, c, d, m = "MV", b_time, b_timese, mp_type,  fix_num = 0, r = 1, s = 1, n_point = 4000, upper = 10){
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
    res = homot_scr(ts, fix_num, c, d, m, b_time, b_timese, mp_type, r = r, s = s, n_point = n_point, upper = upper)
    return(res)
  } else{
    return(stop("Invalid option!"))
  }
}
