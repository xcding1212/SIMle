source(paste(getwd(), "/R/SIMle.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.db1-20.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.plot.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.original_code.v1.R", sep = ""))




#' Visualization of simultaneous confidence region (SCR)
#'
#' @param scr_df the result of estimation 
#' @param ops select type of estimation."nfix" refers to no fix estimation. "fixt" indicates fix time t estimation. 
#'             "fixx" represents fix variate estimation
#' @param title give the title for plot 
#' @param lower give the lower bound for scale limits, the default is -1.3
#' @param upper give the upper bound for scale limits, the default is 1.3
#' @return the plot shows estimated function and its simultaneous confidence region (SCR)
#' @export
#'
#' @examples
#' generate_nAR1 = function(n, v){
#'  ts = c()
#'  w = rnorm(n, 0, 1/v)
#'  x_ini = runif(1,0,1)
#'  for(i in 1:n){
#'    if(i == 1){
#'      ts[i] = sin(2*pi*(i/n))*exp(-x_ini^2)  + w[i] #
#'    } else{
#'      ts[i] = sin(2*pi*(i/n))*exp(-ts[i-1]^2) + w[i]
#'    }
#'  }
#'  return(ts)
#' }
#' ts  = generate_nAR1(27, 1) #change sample size in real case. 
#' res_esti = fix.SCR(ts, 1, 1, m = "MV", "Legen", "Legen", "algeb", "fixt", 0.6, r = 1)
#' scr.plot(res_esti[[1]], "fixt")

scr.plot = function(scr_df, ops, title = "", lower = -1.3, upper = 1.3){
  res = plot_scr(scr_df, ops, title = title, lower = lower, upper = upper)
  return(res)
}
