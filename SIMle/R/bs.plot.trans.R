
source(paste(getwd(), "/R/SIMle.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.db1-20.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.plot.v1.R", sep = ""))
source(paste(getwd(), "/R/SIMle.original_code.v1.R", sep = ""))

#' Plots of mapping basis
#' @description This function generates the plot of first k basis function.
#' @param type type indicates which type of basis is used
#' @param mp_type select type of mapping function, "algeb" indicates algebraic mapping on the real line. "logari" represents logarithmic mapping on the real line
#' @param k The k is the number of basis functions represented (If wavelet are chosen, the real number of basis is 2^k. If
#'     Cspli is chosen, the real number of basis is k - 2 + or)
#' @param upper the upper bound for basis domain, the default is 10      
#' @param s s is a positive scaling factor, the default is 1   
#' @param or indicates the order of spline and only used in Cspli type, default is 4 which indicates cubic spline
#' @param title give the title for the basis plot
#' @return The plot of 1 to k basis functions
#' @export
#'
#' @examples
#' bs.plot.trans("Legen", "algeb", 2)


bs.plot.trans <- function(type, mp_type, k, upper = 10, s = 1, or = 4, title = ""){
  # library(ggplot2)
  # coeffi = Chebyshev_coeff(n, kind = 1)
  
  pos_map = c( "algebp", "logarip")  
  aux_basis = bs.gene.trans(type, mp_type, k)
  
  n_esti = dim(aux_basis)[1]
  
  if(mp_type %in% pos_map){
    aux_x = seq(0, upper, length.out = n_esti)
  }else{
    aux_x = seq(-upper, upper, length.out = n_esti)
  }
  
  if(type == "Cspli"){
    df_basis = Cspline_table(k) # true is d - 2 + or 
    k = dim(df_basis)[2]
  }
  
  wavelet_basis = c("db1", "db2", "db3", "db4", "db5",
                    "db6", "db7", "db8", "db9", "db10",
                    "db11", "db12", "db13", "db14", "db15",
                    "db16", "db17", "db18", "db19", "db20",
                    "cf1", "cf2", "cf3", "cf4", "cf5"
  )  
  if(type %in% wavelet_basis){
    k = 2^k
  }
  
  if(type == "tri"){
    k = 2*k - 1
  } else{
    k = k
  }
  
  x = rep(aux_x, k)
  
  aux = c()
  for(i in 1:k){
    aux = c(aux, aux_basis[,i])
  }
  aux_order = as.factor(rep(1:k, each = n_esti))
  
  f.df = data.frame(x = x, new = aux, order = aux_order)

  theme_update(plot.title = element_text(hjust = 0.5))
  p1 <- ggplot(f.df, aes(x=x, y = new, group = order, colour = order))+ geom_line()  + ggtitle(title) +
    xlab("") + ylab("") + scale_colour_discrete(name  ="order")+theme(plot.title = element_text(size=18, face="bold"),
                                                                      legend.text=element_text(size=24, face = "bold"),
                                                                      axis.text.x = element_text(face="bold", color="#993333",                                                                                                                              size=22, angle=0),
                                                                      axis.text.y = element_text(face="bold", color="#993333",size=22, angle=0),
                                                                      axis.title.x=element_text(size=22,face='bold'),
                                                                      axis.title.y=element_text(angle=90, face='bold', size=22),
                                                                      legend.title = element_text(face = "bold"))
  return(p1)

}

