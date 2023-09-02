# Cubic-spline Demo
# Write function for C-spline, # c is the number of the cubic spline.In detail, splines library is called.

## ord options default 4   number of basis c-2+r
# ord = 5 we need add 4 additional points

# c is at least 2, which is the number of the basis function - 2, the true number of basis is c+2

phi_csp <- function(cspline, beta, b){
  c = length(cspline)
  b_res = list()
  for(i in 0:b){
    B.aux = matrix(c(rep(0, c*i), cspline, rep(0, c*(b-i))), ncol = 1)
    b_res[[i+1]] = as.numeric(t(beta)%*%B.aux)
  }
  
  return(b_res)
}

csp_basis.f <- function(b.table, x){
  return(as.numeric(b.table[which(abs(b.table[,1]-x) == min(abs(b.table[,1]-x)))[1], -1]))
}

Cspline_table = function(c, n=1000, or = 4){
  #library(splines)
  x = seq(0, 1, length=n)
  knots=c(rep(0, or - 1), seq(0,1, length.out = c), rep(1, or - 1)) ## need to add three additional knots at the two ends;
  B=splineDesign(knots, x, ord = or)                ## default order: ord =4, corresponds to cubic splines
  csptable = cbind(as.matrix(x, ncol = 1), B)
  inte = csp_basis.f(csptable, 1/10000)^2*(1/10000)
  for(i in 2:10000){
    inte = inte + csp_basis.f(csptable, i/10000)^2*(1/10000)
  }
  
  for(i in 2:dim(csptable)[2]){
    csptable[, i] = sqrt(1/inte)[i-1]*csptable[, i]
  }
  B = csptable[, -1]
  return(B)
}



# c is number of basis function, x is the inputs value. 
Spline_basis <- function(k, x_val, df){ # true k is c-2+or.    k = c-2+or   
  n = dim(df)[1]
  x = seq(0,1, length.out = n)
  return(df[unique(which(abs(x - x_val) == min(abs(x -x_val)))),])
}


