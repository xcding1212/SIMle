# visualized result

source(paste(getwd(), "/R/SIMle.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.db1-20.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.plot.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.original_code.v1.R", sep = ""))




#' Visualization of estimation
#'
#' @param res_esti the result of estimation 
#' @param ops select type of estimation."nfix" refers to no fix estimation. "fixt" indicates fix time t estimation. 
#'             "fixx" represents fix variate estimation
#' @param mp_type select type of mapping function, "algeb" indicates algebraic mapping on the real line. "logari" represents logarithmic mapping on the real line
#' @param title give the title for plot 
#' @param lower give the lower bound for scale limits, the default is -1.3
#' @param upper give the upper bound for scale limits, the default is 1.3
#' @param domain The upper bound for the variate basis domain. The default value is 10. When "algeb" or "logari" is chosen, the domain is automatically set from -domain to domain. 
#'
#' @return the plot shows estimated function 
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
#' ts  = generate_nAR1(200, 1) # change sample size in real case
#' res_esti = fix.fit(ts, 5, 2, "Legen", "Legen", "algeb", "fixt", 0.1)
#' fit.plot(res_esti[[1]], "fixt", "algeb")

fit.plot <- function(res_esti, ops, mp_type, title="", lower = -1.3, upper = 1.3, domain = 10){
  res = plot_esti(res_esti, ops, mp_type = mp_type, title = title, lower = lower, upper = upper, domain = domain)
  return(res)
}
