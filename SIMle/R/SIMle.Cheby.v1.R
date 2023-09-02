# Chebyshev basis function (kind 1)
# input is the coefficients of polynomial order by (1, x, x^2, x^3,...) until the highest order. The out put is this polynomial multiple by x.
move_order = function(poly){ # get the coefficient when polynomial multiple 1 times of x
  return(c(0,poly))
}

# we can get the plot from the coefficients
# input is n (the number of basis function)
# output is the coefficients of polynomial order by (1, x, x^2, x^3,...) until the highest order.


# No normalized, kind only can be taken 1 and 2.
# The Chebyshev polynomials of the first kind are a special case of the Jacobi polynomials where alpha=beta=-1/2


Chebyshev_coeff = function(n, kind = 1){

  p_c = list()
  p_c[1] = c(1)
  if(n == 1){
    return(p_c[[1]])
  } else if(n==2){
    if(kind == 1){
      p_c[[2]] = c(0,1)
      return(p_c)
    } else{
      p_c[[2]] = c(0,2)
      return(p_c)
    }

  } else{
    if (kind == 1){
      for(i in 3:n){
        aux.i = i - 1
        p_c[[2]] = c(0,1)
        p_c[[i]] = 2*move_order(p_c[[i-1]]) - c(p_c[[i-2]], 0, 0)
      }
    } else{
      for(i in 3:n){
        aux.i = i - 1
        p_c[[2]] = c(0,2)
        p_c[[i]] = 2*move_order(p_c[[i-1]]) - c(p_c[[i-2]], 0, 0)
      }
    }
  }
  return(p_c)

}


# input is n (the number of basis function). coeffi is the list and contain the coefficients of polynomial order by (1, x, x^2, x^3,...) until the highest order.
poly_val = function(coeffi, x){
  for(i in 1:length(coeffi)){
    aux = 0
    for(j in 1:length(coeffi[[i]])){
      aux = aux + coeffi[[i]][j]*(2*x-1)^(j - 1)
    }
    coeffi[[i]] = aux
  }
  return(coeffi)
}

# c is number of basis function, x is the inputs value.
Chebyshev_basis = function(c, x, kind = 1){

  aux_li = Chebyshev_coeff(c, kind)
  res = c()
  res.aux = unlist(poly_val(aux_li, x))
  normcos = c()

  for(i in 1:c){
    fx = function(x){
      aux = 0
      len = length(Chebyshev_coeff(c)[[i]])
      for(j in 1:len){
        aux = aux + Chebyshev_coeff(c)[[i]][j]*((2*x-1)^(j-1))
      }
      return(aux^2)
    }
    normcos[i] = sqrt(1/integrate(fx, 0, 1)[[1]])
  }

  for(n in 1:c){
    res[n] = res.aux[n]*normcos[n]
  }
  return(res)
}

