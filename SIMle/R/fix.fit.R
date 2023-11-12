# Estimate coefficients for dynamic systems


source(paste(getwd(), "/R/SIMle.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.db1-20.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.plot.v1.R", sep = ""))
# source(paste(getwd(), "/R/SIMle.original_code.v1.R", sep = ""))




#' User-specified estimation of nonlinear time series regression
#' @description This function estimates nonlinear time series regression by sieve methods 
#'
#' @param ts ts is the data set which is a time series data typically
#' @param c number of basis for time input
#' @param d number of basis for variate input 
#' @param b_time type of basis for time input 
#' @param b_timese type of basis for variate input
#' @param mp_type select type of mapping function, "algeb" indicates algebraic mapping on the real line. "logari" represents logarithmic mapping on the real line
#' @param type select type of estimation."nfix" refers to no fix estimation. "fixt" indicates fix time t estimation. 
#'             "fixx" represents fix variate estimation
#' @param fix_num fix_num indicates the use of fixed-value nonlinear time series regression. The default value is 0, which is employed for non-fixed estimation. 
#'        If "fixt" is chosen, it represents a fixed time value. Otherwise, if not selected, it pertains to a fixed variate value
#' @param r indicates number of variate 
#' @param s s is a positive scaling factor, the default is 1
#' @param n_esti number of points for estimation, the default is 2000
#' @param upper upper The upper bound for the variate basis domain. The default value is 10. When "algeb" or "logari" is chosen, the domain is automatically set from -upper to upper 
#'
#'
#' @return If "nfix" is selected, the function returns a list where each element is a matrix representing the estimation function in two dimensions. Otherwise, 
#'         if "nfix" is not selected, the function returns a list where each element is a vector representing the estimation function.
#' @export



fix.fit <- function(ts, c, d, b_time, b_timese, mp_type, type, fix_num = 0, r = 1, s = 1, n_esti = 2000, upper = 10){
  basis_candi = c("Legen","Cheby","tri", "cos", "sin", "Cspli", "db1", "db2", "db3", "db4", "db5",
                    "db6", "db7", "db8", "db9", "db10",
                    "db11", "db12", "db13", "db14", "db15",
                    "db16", "db17", "db18", "db19", "db20",
                    "cf1", "cf2", "cf3", "cf4", "cf5"
  )

  if((b_time %in% basis_candi) && (b_timese %in% basis_candi)){
    if(type == "nfix"){
      res = general_esti(ts, c, d, b_time, b_timese, mp_type, r, s, n_esti = n_esti, upper = upper)
      return(res)
    } else if(type == "fixt"){
      res = fix_t_esti(ts, c, d, fix_num, b_time, b_timese, mp_type, r = r, s = s, upper = upper)
      return(res)
    } else if(type == "fixx"){
      res = fix_x_esti(ts, c, d, fix_num, b_time, b_timese, mp_type, r = r, s = s, upper = upper)
        return(res)
    } else{
      return(stop("Invalid option!"))
    }
    
  } else{
    return(stop("Invalid option!"))
  }
}
