
source(paste(getwd(), "/R/SIMle.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.db1-20.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.plot.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.original_code.v1.R", sep = ""))





#' User-specified creation of a Simultaneous Confidence Region (SCR) for the estimated function
#' @description This function generates a Simultaneous Confidence Region (SCR) for the estimated function
#' @param ts ts is the data set which is a time series data typically
#' @param c the maximum value of number of basis for time input 
#' @param d the maximum value of number of basis for variate input
#' @param m the window size for the simultaneous confidence region procedure, with the default being 'MV,' which stands for the Minimum Volatility method
#' @param b_time type of basis for time input 
#' @param b_timese type of basis for variate input
#' @param mp_type select type of mapping function, "algeb" indicates algebraic mapping on the real line. "logari" represents logarithmic mapping on the real line
#' @param type select type of estimation."fixt" indicates fixed time t value. "fixx" represents fix variate value 
#' @param fix_num fix_num indicates the use of fixed-value nonlinear time series regression. If "fixt" is chosen, it represents a fixed time value. Otherwise, 
#'        if not selected, it pertains to a fixed variate value
#' @param r indicates number of variate 
#' @param s s is a positive scaling factor, the default is 1
#' @param n_point number of points for SCR, the default is 4000
#' @param upper upper The upper bound for the variate basis domain. The default value is 10. When "algeb" or "logari" is chosen, the domain is automatically set from -upper to upper 
#'
#'
#' @return A list containing dataframes with three columns each. The first column corresponds to input values. The second column contains values of the estimated function, along with 
#'         their upper and lower bounds. The third column is a factor that indicates the types associated with the values in the second column.
#' @export


fix.SCR <- function(ts, c, d, m = "MV", b_time, b_timese, mp_type, type, fix_num = 0, r = 1, s = 1, n_point = 4000, upper = 10){
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
    if(type == "fixt"){
      res = fix_t_scr(ts, fix_num, c, d, m, b_time, b_timese, mp_type, r = r, s = s, n_point = n_point, upper = upper)
      return(res)
    } else if(type == "fixx"){
      res = fix_x_scr(ts, fix_num, c, d, m, b_time, b_timese, mp_type, r = r, s = s, n_point = n_point, upper = upper)
      return(res)
    } else{
      return(stop("Invalid option!"))
    }
    
  } else{
    return(stop("Invalid option!"))
  }
}



