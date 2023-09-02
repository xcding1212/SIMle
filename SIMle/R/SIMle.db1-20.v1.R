# Daubechies Wavelet 1 and some general functions to be used for Daubechies, the basis of wavelet is formulated by the method of Meyer (S.2 in the paper)
# However, we also try the S.1 for estimating coefficients and choose the J_0 = 0.
#library(ggplot2)


# find the minimal value correspond, inter is the interval created, valtable is the db2 value table which is generated in advance.upper is the upper bound of basis. val express the input x.
w_find = function(valtable, inter, upper, val){
  if(val < 0 | val >= upper){
    return(0)
  } else{
    return(valtable[which(abs(inter-val) == min(abs(inter-val)))])
  }
}

db1_f = function(t){
  return(ifelse(t<1 & t>=0, 1,0))
}


wavelet_kth_b = function(k, point, ops){ # the k-th basis in 2^k total basis, n refers number of points
  w = k
  if(ops == "db1"){
    dbt = valdb(w,psi.f = 0, 1, point)
    n = dim(dbt)[1]
    x = seq(0, 1, length = n)
    
    df = data.frame(x)
    for (i in 2:(2^w + 1)){
      res = as.data.frame(dbt[, i]) #unlist(poly_val(coeffi, i))
      df = cbind(df, res)
    }
    return(data.matrix(df[,-1]))
  } else{
    #library(stringr)
    filename = paste(ops,"_fa_table", sep = "")
    aux_str = str_split(ops, "")[[1]]
    if(aux_str[1] == "d"){
      if(length(aux_str) == 3){
        c = as.numeric(aux_str[3])
      }else
      {
        c = as.numeric(paste(aux_str[3], aux_str[length(aux_str)], sep = ""))
      }
      
    } else{
      c = as.numeric(aux_str[3])*3
    }
    
    #library(RCurl)
    x <- getURL(paste("https://raw.githubusercontent.com/xcding1212/Sie2nts/main/db_table/", filename, sep = ""))
    dbt = read.csv(text = x)
    dbt = dbt[,2]
    t=seq(0, 2*c - 1, length.out = length(dbt))
    n = length(t)
    x =seq(0,1, length.out = n)
    df = data.frame(x)
    for(h in 0:(2^w-1)){
      p=rep(0, n)
      for (k in 1:n){
        for (l in -70:70){
          p[k]=p[k]+ w_find(dbt, t, 2*c - 1, (2^(w)*(x[k]+l)-h))
        }
        p[k]=2^(w/2)*p[k]
      }
      
      df=cbind(df, data.frame(value = p))
      p=rep(0, n)
    }
    db.leng = seq(0, 1, length.out = point)
    new_x = c()
    for (i in 1:point){
      new_x[i] = unique(which(abs(df[,1] - db.leng[i]) == min(abs(df[,1] - db.leng[i]))))
    }
    df = df[new_x, -1]
    return(data.matrix(df))
  }
}


# res = wavelet_kth_b(1, 2000, "db3")



# This basis is formulated by Meyer for Daubechies1
valdb1 = function(w, n){
  x =seq(0,1, length.out = n)
  df.db1 = data.frame(x = x)

  for(h in 0:(2^w-1)){
    p=rep(0, n);
    for (k in 1:n){
      for (l in 0:70){
        p[k]=p[k]+ db1_f((2^(w)*(x[k]+l)-h))
      }
      p[k]=2^(w/2)*p[k]
    }
    df.db1[,h+2] = p
  }
  return(df.db1)
}


# get the res with different w.

# db_number represent which order of Daubechies to be used. 1-20 right now could be chosen. w indicates the number of basis functions.

valdb = function(w, psi.f=0, db_number, len.n){ # len.n indicates the point chosen
  if(db_number == 1){
    return(valdb1(w, len.n))
  } else{
    n = length(psi.f)
    x =seq(0,1, length.out = n)
    df.db = data.frame(x)
    t=seq(0, 2*db_number - 1, length.out = n)

    for(k in 0:(2^w-1)){
      p=rep(0, n);
      for (ind in 1:n){
        for (l in 0:70){
          p[ind]=p[ind]+ w_find(psi.f, t, 2*db_number - 1, (2^(w)*(x[ind]+l)-k))
        }
        p[ind]=(2^(w/2))*p[ind]
      }
      df.db[,k+2] = p
    }
    return(df.db)
  }

}


library(stringr)
library(RCurl)

db_table = function(w, ops){ # c indicate the order of db and w indicate the number of basis(the true basis is 2^w).
  if(ops == "db1"){
    point = 10000
    dbt = valdb(w,psi.f = 0, 1, point)
    n = dim(dbt)[1]
    x = seq(0, 1, length=n)
    
    df = matrix(nrow = n) 
    for (i in 2:(2^w + 1)){
      res = dbt[, i] #unlist(poly_val(coeffi, i))
      df = cbind(df, res)
    }
    df = df[,-1]
    return(df)
  } else{
    #library(stringr)
    filename = paste(ops,"_fa_table", sep = "")
    aux_str = str_split(ops, "")[[1]]
    if(aux_str[1] == "d"){
      if(length(aux_str) == 3){
        c = as.numeric(aux_str[3])
      }else
      {
        c = as.numeric(paste(aux_str[3], aux_str[length(aux_str)], sep = ""))
      }
      
    } else{
      c = as.numeric(aux_str[3])*3
    }
    
    # library(RCurl)
    
    x <- getURL(paste("https://raw.githubusercontent.com/xcding1212/Sie2nts/main/db_table/", filename, sep = ""))
    dbt = read.csv(text = x)
    
    dbt = dbt[,2]
    n = length(dbt)
    df = matrix(nrow = n) 
    t=seq(0, 2*c - 1, length.out = length(dbt))
    n = length(t)
    x =seq(0,1, length.out = n)
    for(h in 0:(2^w-1)){
      p=rep(0, n)
      for (k in 1:n){
        for (l in -70:70){
          p[k]=p[k]+ w_find(dbt, t, 2*c - 1, (2^(w)*(x[k]+l)-h))
        }
        p[k]=2^(w/2)*p[k]
      }
      
      df= cbind(df, p)
      p=rep(0, n)
    }
    df = df[,-1]
    return(df)
  }
}



# c is number of basis function, x is the inputs value. 
Wavelet_basis <- function(k, x_val, df){ # true k is c-2+or.    k = c-2+or   
  n = dim(df)[1]
  x = seq(0,1, length.out = n)
  return(df[unique(which(abs(x - x_val) == min(abs(x -x_val)))),])
}

