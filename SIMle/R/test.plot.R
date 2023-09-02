
source(paste(getwd(), "/R/SIMle.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.db1-20.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.plot.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.original_code.v1.R", sep = ""))



#' the plot for simultaneous confidence region(SCR)
#'
#' @param res_esti the result of estimation and simultaneous confidence region
#' @param ops select type of estimation."nfix" refers to no fix estimation. "fixt" indicates fix time t estimation. 
#'             "fixx" represents fix time series estimation
#' @param title give the title for plot 
#'
#' @return the plot shows estimated function and its simultaneous confidence region(SCR)
#' @export
#'
#' @examples

#' Visulization of Simultaneous Confidence Region(SCR) for test result 
#'
#' @param df the result of test (estimated function under null and Simultaneous Confidence Region (SCR) )
#' @param type specify type of test, "homot" represents time-homogeneity test. "separa" is separability test
#' @param ops select type of estimation."nfix" refers to no fix estimation. "fixt" indicates fix time t estimation. 
#'             "fixx" represents fix variate estimation
#' @param title give the title for plot 
#' @param lower give the lower bound for scale limits, the default is -1.3
#' @param upper give the upper bound for scale limits, the default is 1.3
#' @return the plot shows test estimated function and simultaneous confidence region (SCR)
#' @export
#'
#' @examples
#' res = auto.sep.test(ts, 2, 2, "tri", "db4", "algeb","fixx", "Kfold", k = 4, fix_num = 0.1, r = 2)
#' test.plot(res[[1]], "separa", "fixx")
#' test.plot(res[[2]], "separa", "fixx")


test.plot = function(df, type, ops = "", title = "", lower = -1.3, upper = 1.3){
  if(type == "homot"){
    res = homot_plot(df, title = title, lower = lower, upper = upper)
    return(res)
  } else if (type == "separa"){
    res = sep_plot(df, ops, title = title, lower = lower, upper = upper)
    return(res)
  } else{
    return(stop("Invalid option!"))
  }
  
}




