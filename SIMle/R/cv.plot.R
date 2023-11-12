# visualized result

source(paste(getwd(), "/R/SIMle.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.db1-20.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.plot.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.original_code.v1.R", sep = ""))




#' Visualization of the cross-validation results
#' 
#' @param cv_m give the cross validation data frame
#' @param title give the title for plot 
#' 
#' @return the plot shows cross validation result (3D)
#' @export
 
cv.plot = function(cv_m, title=""){
  res = fit.plot.cvm(cv_m, title = title)
  return(res)
}
