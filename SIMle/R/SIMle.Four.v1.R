# c indicate the number of tri basis function. The true number of basis functions are 2*c - 1. For the rest basis, c expresses the number of basis function.
# Compared with version1, F_general() is removed

F_tri = function(c, x){
  aux = c()
  if(c == 1){
    aux[c] = 1
    return(aux)
  }
  ind = 1
  aux[ind] = 1
  ind = ind + 1
  for(i in 1:(c-1)){
    aux[ind] = sqrt(2)*sin(2*i*pi*x)
    aux[ind+1] = sqrt(2)*cos(2*i*pi*x)
    ind = ind + 2
  }
  return(aux)

}

F_cosPol = function(c, x){
  aux = c()
  if(c == 1){
    aux[c] = 1
    return(aux)
  }
  for(i in 1:c){
    if(i == 1){
      aux[i] = 1
    } else{
      aux[i] = sqrt(2)*cos((i-1)*pi*x)
    }
  }
  return(aux)

}

F_sinPol = function(c, x){
  aux = c()
  for(i in 1:c){
    aux[i] = sqrt(2)*sin(i*pi*x)

  }
  return(aux)
}

# F_general = function(c, x){
#   aux = c()
#   for(i in 1:c){
#     if (i == 1){
#       aux[i] = 1
#     } else if (i %% 2 == 0){
#       aux[i] = sqrt(2)*sin(i*pi*x)
#     } else{
#       aux[i] = sqrt(2)*cos((i-1)*pi*x)
#     }
#   }
#   return(aux)
# }

# If the option parameter equals tri, it means we choose trigometric basis, cos means cospol, sin means sinpol.
select.basis = function(c,x, ops = "tri"){
  if(ops == "tri"){
    return(F_tri(c, x))
  } else if (ops == "cos"){
    return(F_cosPol(c,x))
  } else{
    return(F_sinPol(c,x))
  }
}


# tri --> it gives you (2*d - 1) basis 

