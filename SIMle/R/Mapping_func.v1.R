# algebraic mapping R to [0,1]
algeb_u <- function(x, s){
  return(x/(2*sqrt(x^2+s^2)) + 1/2) 
}
algeb_up <- function(x, s){
  return((s^2)/(2*((x^2+s^2)^(3/2))))
}



# algebraic mapping R+ to [0,1]
algeb_u_posi <- function(x, s){
  return(x/(x+s)) 
}
algeb_up_posi <- function(x, s){
  return(s/((x+s)^2))
}


# The logarithmic mapping from R to [0, 1]:
logari_u <- function(x, s){
  return((tanh(x/s) + 1)/2) 
}
logari_up <- function(x, s){
  return((1-(tanh(x/s)^2))/(2*s))
}



# The logarithmic mapping from R+ to [0, 1]:
logari_u_posi <- function(x, s){
  return(tanh(x/s)) 
}
logari_up_posi <- function(x, s){
  return((1-(tanh(x/s)^2))/(s))
}






