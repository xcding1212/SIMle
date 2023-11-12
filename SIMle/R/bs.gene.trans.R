
source(paste(getwd(), "/R/SIMle.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.db1-20.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.plot.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.original_code.v1.R", sep = ""))




#' Generate Mapping Basis
#'
#' @description this function generates the value of k-th basis function. (The wavelet basis options return the full table)
#' @param type type indicates which type of basis is used
#' @param mp_type select type of mapping function, "algeb" indicates algebraic mapping on the real line. "logari" represents logarithmic mapping on the real line
#' @param k k-th basis function
#' @param upper the upper bound for basis domain, the default is 10 
#' @param s s is a positive scaling factor, the default is 1
#' @param n_esti the number of values got from k-th basis function, the default is 500
#' @param c c only used in Cspli which indicates the total number of knots to generate, the default is 10, c should not be less than k.(for splines, the true
#'        number of basis is c-2+or)
#' @param or indicates the order of spline and only used in Cspli type, default is 4 which indicates cubic spline
#' @return A matrix in which the k-th column corresponds to the values of the k-th mapped basis function
#' @references \[1] Chen, Xiaohong. “Large Sample Sieve Estimation of Semi-Nonparametric Models.” Handbook of Econometrics, 6(B): 5549–5632,2007.
#' @export
#'
#' @examples
#' bs.gene.trans("Legen", "algeb", 5)


bs.gene.trans <- function(type, mp_type, k, upper = 10, s = 1, n_esti = 500, c = 10, or = 4){
  # typically, scalar function is s. 
  d = k
  
  pos_map = c( "algebp", "logarip") 
  if(mp_type %in% pos_map){
    x_i = seq(0, upper, length.out = n_esti)
  }else{
    x_i = seq(-upper, upper, length.out = n_esti) 
  }
  
  
  df_basis = 0
  if(type == "Cspli"){
    df_basis = Cspline_table(d) # true is d - 2 + or 
    d = dim(df_basis)[2]
  }
  
  # change for wavelet
  wavelet_basis = c("db1", "db2", "db3", "db4", "db5",
                    "db6", "db7", "db8", "db9", "db10",
                    "db11", "db12", "db13", "db14", "db15",
                    "db16", "db17", "db18", "db19", "db20",
                    "cf1", "cf2", "cf3", "cf4", "cf5"
  )  
  if(type %in% wavelet_basis){
    df_basis = db_table(d, type)
    d = dim(df_basis)[2]
  }
  
  if(type == "tri"){
    d_tri = 2*d - 1
  } else{
    d_tri = d
  }
  
  basis_xi = matrix(ncol = d_tri) 
  
  for(i in 1:n_esti){
    basis_xi = rbind(basis_xi, sqrt(mp_selection(mp_type, x_i[i], s)[2])*select_basis_timese(d, mp_selection(mp_type, x_i[i], s)[1], type, df_basis = df_basis))
  }
  basis_xi = basis_xi[-1, ]
  basis_xi = as.matrix(basis_xi, ncol = d_tri)
  
  return(basis_xi)
}


