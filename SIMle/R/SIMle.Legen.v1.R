
# This is the first function in the SIMle package and here are several tasks remain. 

# 
# 2. Created guideline 
# 3. Finish function to SIMle

# Remaining: 
# 1. Addiing Cspline and wavelet for basis 
# 2. Change matrix structure, now only allow r =1 but r can be any number. 
# 4. Not including real data.
# 1. There are 4 mapped functions, x should close to 0 to apply this approach.




# Legendre
move_order <- function(poly){ # get the coefficient when polynomial multiple 1 times of x
  return(c(0,poly))
}
# we can get the plot from the coefficients
# input is n (the number of basis function)
# output is the coefficients of polynomial order by (1, x, x^2, x^3,...) until the highest order.
# always normalized
legendre_coeff <- function(n){
  p_c = list()
  p_c[1] = c(1)
  if(n == 1){
    return(p_c[[1]])
  } else if(n==2){
    p_c[[2]] = c(0,sqrt(3))
    return(p_c)
    
  } else{
    for(i in 3:n){
      aux.i = i - 1
      p_c[[2]] = c(0, 1)
      p_c[[i]] = ((2*aux.i-1)/aux.i)*move_order(p_c[[i-1]]) - ((aux.i-1)/aux.i)*c(p_c[[i-2]], 0, 0)
    }
    p_n = list(p_c[[1]])
    for (i in 2:length(p_c)){
      p_n[[i]] = sqrt((2*i-1))*p_c[[i]]
    }
    # p_n[[1]] = 1/sqrt(2)
    return(p_n)
  }
}
# x \in [0,1]
# input is n (the number of basis function). coeffi is the list and contain the coefficients of polynomial order by (1, x, x^2, x^3,...) until the highest order.
poly_val <- function(coeffi, x){
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
Legendre_basis <- function(c, x){
  aux_li = legendre_coeff(c)
  return(unlist(poly_val(aux_li, x)))
}



