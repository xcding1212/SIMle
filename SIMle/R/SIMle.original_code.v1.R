
source(paste(getwd(), "/R/Mapping_func.v1.R", sep = ""))


# Adding Cspline and wavelets
select_basis_timese <- function(d, x, b_timese, df_basis){ # d: number-th of basis
                                                # x: value corresponding 
                                                # b_timese: type of basis used
  # we should confirm that b_timese is a valid option in this function 
  wavelet_basis = c("db1", "db2", "db3", "db4", "db5",
                    "db6", "db7", "db8", "db9", "db10",
                    "db11", "db12", "db13", "db14", "db15",
                    "db16", "db17", "db18", "db19", "db20",
                    "cf1", "cf2", "cf3", "cf4", "cf5"
  )  
  
  if(b_timese == "Legen"){
    return(Legendre_basis(d, x))
    
  }else if(b_timese == "Cheby"){
    return(Chebyshev_basis(d, x))
    
  }else if(b_timese %in% c("tri", "cos", "sin")){
    return(select.basis(d, x, ops = b_timese))
    
  }else if(b_timese == "Cspli"){
    return(Spline_basis(d, x, df_basis)) # this d is c for spline and the true number of basis is c-2+or
    
  }else{ # for wavelet
    return(Wavelet_basis(d, x, df_basis))
  }
}

# Mapping selection 

mp_selection <- function(mp_type, x, s = 1){ # 
  if(mp_type == "algeb"){
    return(c(algeb_u(x, s), algeb_up(x, s))) # first is mapping function, second is derivative of mapping function 
    
  }else if(mp_type == "logari"){
    return(c(logari_u(x, s), logari_up(x, s)))
    
  #else if(mp_type == "algebp"){
    #return(c(algeb_u_posi(x, s), algeb_up_posi(x, s)))
  #}else if(mp_type == "logarip"){
  #  return(c(logari_u_posi(x, s), logari_up_posi(x, s)))
    
  }else{
    return(stop("Invalid option!"))
  }
  
}

# selecting best cv results
best_cd.auto.fit <- function(ts, c, d, b_time, b_timese, mp_type, ops, r = 1, s = 1, per = 0, k = 0){
  if(ops == "AIC"){
    res = cd_selection(ts, c, d, b_time, b_timese, mp_type, ops, r = r, s = s)
    #res = res[which(res[,3] != NA), ]
    return(res[which(res[,3] == min(res[,3])), c(1,2)])
    
  } else if(ops == "BIC"){
    res = cd_selection(ts, c, d, b_time, b_timese, mp_type, ops, r = r, s = s)
    #res = res[which(res[,3] != NA), ]
    return(res[which(res[,3] == min(res[,3])), c(1,2)])
    
  } else if(ops == "CV"){
    res = cd_selection(ts, c, d, b_time, b_timese, mp_type, ops, per = per, r = r, s = s)
    #res = res[which(res[,3] != NA), ]
    return(res[which(res[,3] == min(res[,3])), c(1,2)])
    
  } else if(ops == "Kfold"){
    res = cd_selection(ts, c, d, b_time, b_timese, mp_type, ops, k = k, r = r, s = s)
    #res = res[which(res[,3] != NA), ]
    return(res[which(res[,3] == min(res[,3])), c(1,2)])
  } else{
    return(stop("Invalid option!"))
  }
}


cd_selection <- function(ts, c, d, b_time, b_timese, mp_type, ops, per = 0, k = 0, r = 1, s = 1){ 
  # choosing the minimum value for c and d
  # ops for selecting method, AIC, BIC, cross validation(train and test set), k-fold
  # per only for cross-validation 
  # k only for k-folds
  if(ops == "AIC"){
    res = aic(ts, c, d, b_time, b_timese, mp_type, r, s)
    # automatically select c and d
    # res[which(res[,3] == min(res[,3])), c(1,2)] to get the best c and d
    res = res[which(res[,3] != "NA"), ]
    return(res)
  } else if(ops == "BIC"){
    res = bic(ts, c, d, b_time, b_timese, mp_type, r, s)
    res = res[which(res[,3] != "NA"), ]
    # automatically select c and d
    # res[which(res[,3] == min(res[,3])), c(1,2)]
    return(res)
  } else if(ops == "CV"){
    res = cv(ts, c, d, per, b_time, b_timese, mp_type, r, s)
    res = res[which(res[,3] != "NA"), ]
    # automatically select c and d
    # res[which(res[,3] == min(res[,3])), c(1,2)]
    return(res)
  } else if(ops == "Kfold"){
    
    res = kcv(ts, c, d, k, b_time, b_timese, mp_type, r, s)
    res = res[which(res[,3] != "NA"), ]
    # automatically select c and d
    # res[which(res[,3] == min(res[,3])), c(1,2)]
    return(res)
    
  } else{
    return(stop("Invalid option!"))
  }
  
}

esti_beta <- function(timese, c, d, b_time, b_timese, mp_type, r = 1, s = 1){  # b_time indicates type of basis for time t 
                                                     # b_timese indicates name of basis for time series
  n = length(timese)
  Y = matrix(timese[(r+1):n], ncol = 1)
  W = matrix(nrow = n-r, ncol = 1)
  
  # change for Cspline
  df_basis = 0
  if(b_timese == "Cspli"){
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
  if(b_timese %in% wavelet_basis){
    df_basis = db_table(d, b_timese)
    d = dim(df_basis)[2]
  }
  
  if(b_timese == "tri"){
    d_tri = 2*d - 1
  } else{
    d_tri = d
  }
  #}
  
  phi_x = matrix(nrow = n-r, ncol = 1) 
  aux_phi_x = matrix(ncol = d_tri)
  
  
  
  for(j in 1:r){ # j for lag
    for(i in (r+1):n){
      aux_phi_x = rbind(aux_phi_x, sqrt(mp_selection(mp_type, timese[i-j], s)[2])*select_basis_timese(d, mp_selection(mp_type, timese[i-j], s)[1], b_timese, df_basis))
    }
    aux_phi_x = aux_phi_x[-1, ]
    aux_phi_x = as.matrix(aux_phi_x, ncol = d_tri)
    phi_x = cbind(phi_x, aux_phi_x)
    aux_phi_x = matrix(ncol = d_tri)
  }
  
  
  phi_x = matrix(phi_x[,-1], ncol = r*d_tri)
  
  if(b_time == "Cspli"){
    df_basis = Cspline_table(c, n) # true is d - 2 + or 
    c = dim(df_basis)[2]
  }
  
  if(b_time %in% wavelet_basis){
    df_basis_1 = wavelet_kth_b(c, n, b_time)
    c = dim(df_basis_1)[2]
  }
  
  if(b_time == "tri"){
    c = c*2 -1 
  }
  
  for(j in 1:r){
    for(l1 in 1:c){
      if(b_time == "Cspli"){
        aux_basis_t = as.matrix(df_basis[(r+1):n, l1])
      }else if(b_time %in% wavelet_basis){
        aux_basis_t = as.matrix(df_basis_1[(r+1):n, l1])
      }else{
        aux_basis_t = as.matrix(bs.gene(b_time, l1, n)[(r+1):n,])
      }
      
      
      for(l2 in 1:d_tri){
        W = cbind(W, aux_basis_t*matrix(phi_x[, (j-1)*d_tri + l2], ncol = 1))
      }
    }
  }
  
  
  
  W = W[,-1]
  W = as.matrix(W, ncol = r*c*d_tri)
  colnames(W) = NULL
  beta_hat = solve(t(W)%*%W, tol = 1e-40)%*%t(W)%*%Y
  return(list(fit_beta = beta_hat, d_m = W)) # estimation need to change
}


# for Cspli, at least 2, the default is or = 4 represents Cspline,  
# this d is c for spline and the true number of basis is c-2+or. For db and cf, the true number of basis are 2^c and 
# 2^d

# Estimation: overall (3D plot), fix_t(2D plot), fix_x(2D plot) (auto version, automatically choosing c and d)
general_esti <- function(ts, c, d, b_time, b_timese, mp_type, r = 1, s = 1, n_esti = 2000, upper = 10){ # ts for train beta
  
  # typically, scalar function is s. 
  basis_ti = matrix(nrow = n_esti) 
  
  pos_map = c( "algebp", "logarip") 
  if(mp_type %in% pos_map){
    x_i = seq(0, upper, length.out = n_esti)
  }else{
    x_i = seq(-upper, upper, length.out = n_esti) 
  }
  
  
  esti_df = esti_beta(ts, c, d, b_time, b_timese,mp_type, r, s)
  beta_hat = esti_df[[1]]
  W =  esti_df[[2]]
  # n = length(ts)
  
  df_basis = 0
  if(b_timese == "Cspli"){
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
  if(b_timese %in% wavelet_basis){
    df_basis = db_table(d, b_timese)
    d = dim(df_basis)[2]
  }
  
  if(b_timese == "tri"){
    d_tri = 2*d - 1
  } else{
    d_tri = d
  }
  
  basis_xi = matrix(ncol = d_tri) 
  
  for(i in 1:n_esti){
    basis_xi = rbind(basis_xi, sqrt(mp_selection(mp_type, x_i[i], s)[2])*select_basis_timese(d, mp_selection(mp_type, x_i[i], s)[1], b_timese, df_basis = df_basis))
  }
  basis_xi = basis_xi[-1, ]
  basis_xi = as.matrix(basis_xi, ncol = d_tri)
  
  if(b_time == "Cspli"){
    df_basis = Cspline_table(c, n_esti) # true is d - 2 + or 
    c = dim(df_basis)[2]
  }
  
  if(b_time %in% wavelet_basis){
    df_basis_1 = wavelet_kth_b(c, n_esti, b_time)
    c = dim(df_basis_1)[2]
  }
  
  if(b_time == "tri"){
    c = c*2 -1 
  }
  
  for(l1 in 1:c){
    if(b_time == "Cspli"){
      aux_bti = as.matrix(df_basis[, l1])
    }else if(b_time %in% wavelet_basis){
      aux_bti = as.matrix(df_basis_1[, l1])
    }else{
      aux_bti = bs.gene(b_time, l1, n_esti)
    }
    basis_ti = cbind(basis_ti, aux_bti)
  }
  
  basis_ti = as.matrix(basis_ti[,-1], ncol = c)
  colnames(basis_ti) = NULL
  
  # 2 dimensinal integrate
  
  t=matrix(seq(0,1,length.out = n_esti))
  x=matrix(seq(-10, 10,length.out = n_esti))
  m_hat_ij = matrix(rep(0, (n_esti)^2), nrow = n_esti, ncol = n_esti)
  
  B_j = matrix(rep(0, (c*d_tri)^2), nrow = c*d_tri, ncol = c*d_tri)
  
  res_m_hat_ij = list()
  
  for(k in 1:r){
    for(i in 1:n_esti){
      for(j in 1:n_esti){
        ti = matrix(as.numeric(basis_ti[i,]), ncol = 1)
        xj = matrix(as.numeric(basis_xi[j,]), ncol = 1)
        B_j = B_j + kronecker(ti, xj)%*%t(kronecker(ti, xj))
        aux_beta_hat = matrix(beta_hat[((k-1)*c*d_tri+1):(k*c*d_tri), 1], ncol = 1)
        aux_m_hat_ij = t(aux_beta_hat)%*%kronecker(ti, xj)
        m_hat_ij[i, j] = aux_m_hat_ij
      }
    }
    res_m_hat_ij[[k]] = m_hat_ij
  }
  
  
  return(res_m_hat_ij)
  
}

# res = general_esti(timese[1:100], 2,3, "tri", "Legen", "algeb", n_esti = 50)

# Get estimated time series (for prediction)
esti_ts <- function(ts, c, d, b_time, b_timese, mp_type, r=1, s=1, n_esti = 2000, upper = 10){ # in r=1 case 
  n = length(ts)
  fited_res = c()
  
  pos_map = c( "algebp", "logarip") 
  if(mp_type %in% pos_map){
    x_i = seq(0, upper, length.out = n_esti)
  }else{
    x_i = seq(-upper, upper, length.out = n_esti) 
  }
  
  esti_func = general_esti(ts, c, d, b_time, b_timese, mp_type, r, s, n_esti) # estimated function using ts
  aux_esti_func = 0
  for(i in (r+1):n_esti){
    for(j in 1:r){
      aux_fix_inde = unique(which(abs(x_i-ts[i-r]) == min(abs(x_i-ts[i-r]))))[1]
      aux_esti_func = aux_esti_func + esti_func[[j]][i, aux_fix_inde]
    }
    fited_res[i-1] = aux_esti_func
    aux_esti_func = 0 
  }
  
  return(fited_res)
}

# fix t or time series to estimate 
fix_t_esti <- function(timese, c, d, t, b_time, b_timese, mp_type, r=1, s = 1, n_esti = 2000, upper = 10){ # default n_esti = 4000
  # t choose from [0, 1], default of number of estimation is 2000, so t/2000 is the number chosen from [0, 1]
  esti_li = esti_beta(timese, c, d, b_time, b_timese, mp_type, r, s)
  beta_hat = esti_li[[1]]
  df_basis = 0
  
  if(b_timese == "Cspli"){
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
  if(b_timese %in% wavelet_basis){
    df_basis = db_table(d, b_timese)
    d = dim(df_basis)[2]
  }
  
  if(b_timese == "tri"){
    d_tri = 2*d - 1
  } else{
    d_tri = d
  }
  
  basis_x = matrix(ncol = d_tri) 
  basis_t = matrix(nrow = n_esti) 
  b = matrix(nrow = n_esti, ncol = 1)
  
  
  pos_map = c( "algebp", "logarip") 
  if(mp_type %in% pos_map){
    x = seq(0, upper, length.out = n_esti)
  }else{
    x = seq(-upper, upper, length.out = n_esti) 
  }
  
  func_t = seq(0, 1, length.out = n_esti)
  index_t = unique(which(abs(func_t-t) == min(abs(func_t-t))))[1]
  
  
  for(i in 1:n_esti){
    basis_x = rbind(basis_x, sqrt(mp_selection(mp_type, x[i], s)[2])*select_basis_timese(d, mp_selection(mp_type, x[i], s)[1], b_timese, df_basis))
  }
  basis_x = basis_x[-1, ]
  basis_x = matrix(basis_x, ncol = d_tri)
  
  if(b_time == "Cspli"){
    df_basis = Cspline_table(c, n_esti) # true is d - 2 + or 
    c = dim(df_basis)[2]
  }
  
  if(b_time %in% wavelet_basis){
    df_basis_1 = wavelet_kth_b(c, n_esti, b_time)
    c = dim(df_basis_1)[2]
  }
  
  if(b_time == "tri"){
    c = c*2 -1 
  }
  
  for(l1 in 1:c){
    if(b_time == "Cspli"){
      basis_t = as.matrix(df_basis[,l1])
    }else if(b_time %in% wavelet_basis){
      basis_t = as.matrix(df_basis_1[,l1])
    }else{
      basis_t = as.matrix(bs.gene(b_time, l1, n_esti))
    }

    for(l2 in 1:d_tri){
      b = cbind(b, basis_t[index_t]*matrix(basis_x[, l2], ncol = 1))
    }
  }
  
  
  colnames(b) = NULL
  # I_r = diag(r)
  b = as.matrix(b[,-1], ncol = c*d_tri)
  
  # I_r = diag(1)
  # l_j = diag(c*d) 
  # esti = list()
  # for(k in 1:r){
  #  aux_beta_hat = matrix(beta_hat[((k-1)*c*d+1):(k*c*d), 1], ncol = 1)
    # print(dim(aux_beta_hat))
  #  esti[[k]] = t(t(l_j%*%aux_beta_hat)%*%(t(kronecker(b, I_r))))[,1]
  # }
 
  # l_j = diag(c*d) 
  # esti = t(t(l_j%*%beta_hat)%*%(t(kronecker(b, I_r))))[,1]
  # return(esti)
  # return(esti)
  
  esti = list()
  for(k in 1:r){
    aux_beta_hat = matrix(beta_hat[((k-1)*c*d_tri+1):(k*c*d_tri), 1], ncol = 1)
    # esti[[k]] = t(t(l_j%*%aux_beta_hat)%*%(t(kronecker(b, I_r))))[,1]
    esti[[k]] = (b%*%aux_beta_hat)[,1]
  }
  return(esti)
  
}


# library(cIRT)

fix_x_esti <- function(timese, c, d, fix_x, b_time, b_timese, mp_type, r=1, s=1, n_esti = 2000, upper = 10){ 
  # # fix_x choose from [-10, 10], default of number of estimation is 2000, so t/2000 is the number chosen from [-10, 10]
  
  basis_t = matrix(nrow = n_esti) 
  esti_li = esti_beta(timese, c, d, b_time, b_timese, mp_type, r = r, s)
  beta_hat = esti_li[[1]]
  b = matrix(nrow = n_esti, ncol = 1)
  
  pos_map = c( "algebp", "logarip") 
  if(mp_type %in% pos_map){
    x = seq(0, upper, length.out = n_esti)
  }else{
    x = seq(-upper, upper, length.out = n_esti) 
  }
  
  df_basis = 0
  if(b_timese == "Cspli"){
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
  if(b_timese %in% wavelet_basis){
    df_basis = db_table(d, b_timese)
    d = dim(df_basis)[2]
  }
  
  if(b_timese == "tri"){
    d_tri = 2*d - 1
  } else{
    d_tri = d
  }
  
  basis_x = matrix(ncol = d_tri) 
  
  
  for(i in 1:n_esti){
    basis_x = rbind(basis_x, sqrt(mp_selection(mp_type, x[i], s)[2])*select_basis_timese(d, mp_selection(mp_type, x[i], s)[1], b_timese, df_basis))
  }
  basis_x = basis_x[-1, ]
  basis_x = as.matrix(basis_x, ncol = d_tri)
  
  
  if(b_time == "Cspli"){
    df_basis = Cspline_table(c, n_esti) # true is d - 2 + or 
    c = dim(df_basis)[2]
  }
  
  if(b_time %in% wavelet_basis){
    df_basis_1 = wavelet_kth_b(c, n_esti, b_time)
    c = dim(df_basis_1)[2]
  }
  
  if(b_time == "tri"){
    c = c*2 -1 
  }
  
  index_x = unique(which(abs(x-fix_x) == min(abs(x-fix_x))))[1]
  
  for(l1 in 1:c){
    if(b_time == "Cspli"){
      basis_t = as.matrix(df_basis[,l1])
    }else if(b_time %in% wavelet_basis){
      basis_t = as.matrix(df_basis_1[,l1])
    }else{
      basis_t = as.matrix(bs.gene(b_time, l1, n_esti))
    }
    for(l2 in 1:d_tri){
      b = cbind(b, basis_t*basis_x[index_x, l2])
    }
  }
  
  colnames(b) = NULL
  b = as.matrix(b[,-1], ncol = c*d_tri)
  
  # I_r = diag(r)
  # aux_L_j = list()
  # esti = list()
  # for(j in 1:r){
  #   for(k in 1:r){
  #    if(j == k){
  #      aux_L_j[[k]] = diag(c*d) 
  #    }else{
  #      aux_L_j[[k]] = 0*diag(c*d) 
  #    }
      
  #  }
    #L_j = direct_sum(aux_L_j)
    #esti[[j]] = t(t(beta_hat)%*%L_j %*%(t(kronecker(b, I_r))))[,1]
    #aux_L_j = list()
  # }
  
  # I_r = diag(1)
  # l_j = diag(c*d) 
  esti = list()
  for(k in 1:r){
    aux_beta_hat = matrix(beta_hat[((k-1)*c*d_tri+1):(k*c*d_tri), 1], ncol = 1)
    # esti[[k]] = t(t(l_j%*%aux_beta_hat)%*%(t(kronecker(b, I_r))))[,1]
    esti[[k]] = (b%*%aux_beta_hat)[,1]
    }
  return(esti)
}

# Prediction and choosing c and d for auto.fit ((auto version prediction, automatically choosing c and d))
# cross_validation (separate in percentage) 

cross_validation <- function(ts, c, d, per, b_time, b_timese, mp_type, r = 1, s = 1){ # per: percentage for testing set, typically is 0.1
  n = length(ts)
  gcv = 0
  aux.true = ts[(floor(n*(1-per))+1):n]
  esti_df = esti_beta(ts[1:floor(n*(1-per))], c, d, b_time, b_timese,mp_type, r, s)
  
  W = esti_df[[2]]
  beta_hat = esti_df[[1]]
  
  df_basis = 0
  if(b_timese == "Cspli"){
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
  if(b_timese %in% wavelet_basis){
    df_basis = db_table(d, b_timese)
    d = dim(df_basis)[2]
  }
  
  if(b_timese == "tri"){
    d_tri = 2*d - 1
  } else{
    d_tri = d
  }
  
  n_esti = n-r
  b = matrix(nrow = n_esti, ncol = 1)
  
  basis_x = matrix(nrow = n-r, ncol = 1) 
  aux_phi_x = matrix(ncol = d_tri)
  
  for(j in 1:r){ # j for lag
    for(i in (r+1):n){
      aux_phi_x = rbind(aux_phi_x, sqrt(mp_selection(mp_type, ts[i-j], s)[2])*select_basis_timese(d, mp_selection(mp_type, ts[i-j], s)[1], b_timese, df_basis))
    }
    aux_phi_x = aux_phi_x[-1, ]
    aux_phi_x = as.matrix(aux_phi_x, ncol = d_tri)
    basis_x = cbind(basis_x, aux_phi_x)
    aux_phi_x = matrix(ncol = d_tri)
  }
  
  basis_x = matrix(basis_x[,-1], ncol = r*d_tri)
  
  # print(dim(basis_x))
  
  #basis_x = matrix(ncol = d) 

  #for(i in (r+1):n){
  #  basis_x = rbind(basis_x, sqrt(mp_selection(mp_type, ts[i-1], s)[2])*select_basis_timese(d, mp_selection(mp_type, ts[i-1], s)[1], b_timese, df_basis))
  #}
  
  #basis_x = basis_x[-1, ]
  #basis_x = as.matrix(basis_x, ncol = d)
  
  
  if(b_time == "Cspli"){
    df_basis = Cspline_table(c, n) # true is d - 2 + or 
    c = dim(df_basis)[2]
  }
  
  if(b_time %in% wavelet_basis){
    df_basis_1 = wavelet_kth_b(c, n, b_time)
    c = dim(df_basis_1)[2]
  }
  
  if(b_time == "tri"){
    c = c*2 -1 
  }
  
  for(k in 1:r){
    for(l1 in 1:c){
      if(b_time == "Cspli"){
        basis_t = as.matrix(df_basis[(r+1):n, l1])
      }else if(b_time %in% wavelet_basis){
        basis_t = as.matrix(df_basis_1[(r+1):n, l1])
      }else{
        basis_t = bs.gene(b_time, l1, n)[(r+1):n, 1]
      }
      for(l2 in 1:d_tri){
        b = cbind(b, matrix(basis_t*basis_x[, (k-1)*d_tri + l2], ncol = 1)) #  matrix(basis_t*basis_x[, l2], ncol = 1)
      }
    }
  }
  
  colnames(b) = NULL
  b = as.matrix(b[,-1], ncol = r*c*d_tri)
  #I_r = diag(r)
  #l_j = diag(c*d)  
  #aux.esti = t(t(l_j%*%beta_hat)%*%(t(kronecker(b, I_r))))[,1]
  
  esti = list()
  for(k in 1:r){
    aux_beta_hat = matrix(beta_hat[((k-1)*c*d_tri+1):(k*c*d_tri), 1], ncol = 1)
    # esti[[k]] = t(t(l_j%*%aux_beta_hat)%*%(t(kronecker(b, I_r))))[,1]
    
    aux_b = matrix(b[,((k-1)*c*d_tri + 1):(k*c*d_tri)], ncol = c*d_tri)
    esti[[k]] = (aux_b%*%aux_beta_hat)[,1]
  }
  
  aux.esti = rep(0, n_esti)
  for(i in 1:n_esti){
    for(k in 1:r){
      aux.esti[i] = aux.esti[i] + esti[[k]][i] 
    }
  }
  
  MSE = sum((aux.esti[(floor(n*(1-per))+1-r): (n-r)] - aux.true)^2)/length(aux.true)
  
  return(list(CV = MSE, y_hat = aux.esti[(floor(n*(1-per))+1-r): (n-r)]))
}

cv = function(ts, c, d, per, b_time, b_timese, mp_type, r = 1, s = 1){
  res = matrix(nrow = 1, ncol = 3)
  for(i in 1:c){
    for(j in 1:d){
      if(i == 1 & b_time == "Cspli"){
        res = rbind(res, c(i, j, "NA"))
      }else if(j == 1 & b_timese == "Cspli"){
        res = rbind(res, c(i, j, "NA"))
      } else{
        cv_score = cross_validation(ts, i, j, per,b_time, b_timese, mp_type, r = r, s = s)[[1]]
        res = rbind(res, c(i, j, cv_score))
      }
    }
  }
  res = res[-1,]
  return(res)
}

# cross_validation(k-fold for time series) need change

k_fold <- function(ts, c, d, k, b_time, b_timese, mp_type, r = 1, s = 1){
  # k at least 2. 
  aux_c = c
  aux_d = d
  n = length(ts)
  n_esti = n - r
  kcv = 0
  sep_p = (n - n%%k)/k
  b = matrix(nrow = n_esti, ncol = 1)
  n_esti = n
  # basis_x = matrix(ncol = d)
  
  df_basis = 0
  if(b_timese == "Cspli"){
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
  if(b_timese %in% wavelet_basis){
    df_basis = db_table(d, b_timese)
    d = dim(df_basis)[2]
  }
  
  if(b_timese == "tri"){
    d_tri = 2*d - 1
  } else{
    d_tri = d
  }
  
  #for(i in (r+1):n){
  #  basis_x = rbind(basis_x, sqrt(mp_selection(mp_type, ts[i-1], s)[2])*select_basis_timese(d, mp_selection(mp_type, ts[i-1], s)[1], b_timese, df_basis))
  #}
  
  #basis_x = basis_x[-1, ]
  #basis_x = as.matrix(basis_x, ncol = d)
  
  basis_x = matrix(nrow = n-r, ncol = 1) 
  aux_phi_x = matrix(ncol = d_tri)
  for(j in 1:r){ # j for lag
    for(i in (r+1):n){
      aux_phi_x = rbind(aux_phi_x, sqrt(mp_selection(mp_type, ts[i-j], s)[2])*select_basis_timese(d, mp_selection(mp_type, ts[i-j], s)[1], b_timese, df_basis))
    }
    aux_phi_x = aux_phi_x[-1, ]
    aux_phi_x = as.matrix(aux_phi_x, ncol = d_tri)
    basis_x = cbind(basis_x, aux_phi_x)
    aux_phi_x = matrix(ncol = d_tri)
  }
  
  basis_x = matrix(basis_x[,-1], ncol = r*d_tri)
  
  if(b_time == "Cspli"){
    df_basis = Cspline_table(c, n) # true is d - 2 + or 
    c = dim(df_basis)[2]
  }
  
  if(b_time %in% wavelet_basis){
    df_basis_1 = wavelet_kth_b(c, n, b_time)
    c = dim(df_basis_1)[2]
  }
  
  if(b_time == "tri"){
    c = c*2 -1 
  }
  
  for(k_lag in 1:r){
    for(l1 in 1:c){
      if(b_time == "Cspli"){
        basis_t = as.matrix(df_basis[(r+1):n, l1])
      }else if(b_time %in% wavelet_basis){
        basis_t = as.matrix(df_basis_1[(r+1):n, l1])
      }else{
        basis_t = bs.gene(b_time, l1, n)[(r+1):n, 1]
      }
      for(l2 in 1:d_tri){
        b = cbind(b, matrix(basis_t*basis_x[, (k_lag-1)*d_tri +  l2], ncol = 1)) #  matrix(basis_t*basis_x[, l2], ncol = 1)
      }
    }
  }
  
  colnames(b) = NULL
  b = as.matrix(b[,-1], ncol = r*c*d_tri)
  
  #I_r = diag(r)
  #l_j = diag(c*d) 
  
  for(i_sep in 1:(k-1)){
    if(i_sep != (k-1)){
      aux.true = ts[(sep_p*i_sep+1):(sep_p*(i_sep+1))]
      
      esti_df = esti_beta(ts[1:(sep_p*i_sep)], aux_c, aux_d, b_time, b_timese, mp_type, r, s)
      W = esti_df[[2]]
      beta_hat = esti_df[[1]]
      
      
      esti = list()
      for(k_ind in 1:r){
        aux_beta_hat = matrix(beta_hat[((k_ind-1)*c*d_tri+1):(k_ind*c*d_tri), 1], ncol = 1)
        # esti[[k]] = t(t(l_j%*%aux_beta_hat)%*%(t(kronecker(b, I_r))))[,1]
        
        aux_b = matrix(b[,((k_ind-1)*c*d_tri + 1):(k_ind*c*d_tri)], ncol = c*d_tri)
        esti[[k_ind]] = (aux_b%*%aux_beta_hat)[,1]
      }
      
      aux.esti = rep(0, n_esti - r)
      for(i in 1:(n_esti-r)){
        for(k_ind in 1:r){
          aux.esti[i] = aux.esti[i] + esti[[k_ind]][i] 
        }
      }
      
      
      # aux.esti = t(t(l_j%*%beta_hat)%*%(t(kronecker(b, I_r))))[,1]
      
      MSE = sum((aux.esti[(sep_p*i_sep + 1 - r):(sep_p*(i_sep+1) - r)] - aux.true)^2)/length(aux.true)
      kcv = kcv + MSE 
    } else{
      
      aux.true = ts[(sep_p*i_sep + 1):n]
      
      esti_df = esti_beta(ts[1:(sep_p*i_sep)], aux_c, aux_d, b_time, b_timese, mp_type, r, s)
      
      W = esti_df[[2]]
      beta_hat = esti_df[[1]]
      
      # aux.esti = t(t(l_j%*%beta_hat)%*%(t(kronecker(b, I_r))))[,1]
      esti = list()
      for(k_ind in 1:r){
        aux_beta_hat = matrix(beta_hat[((k_ind-1)*c*d_tri+1):(k_ind*c*d_tri), 1], ncol = 1)
        # esti[[k]] = t(t(l_j%*%aux_beta_hat)%*%(t(kronecker(b, I_r))))[,1]
        
        aux_b = matrix(b[,((k_ind-1)*c*d_tri + 1):(k_ind*c*d_tri)], ncol = c*d_tri)
        esti[[k_ind]] = (aux_b%*%aux_beta_hat)[,1]
      }
      
      aux.esti = rep(0, n_esti - r)
      for(i in 1:(n_esti-r)){
        for(k_ind in 1:r){
          aux.esti[i] = aux.esti[i] + esti[[k_ind]][i] 
        }
      }
      # print(aux.esti) from x_(r+1) to x_n
      MSE = sum((aux.esti[(sep_p*i_sep+1-r):(n_esti-r)] - aux.true)^2)/length(aux.true)
      kcv = kcv + MSE 
    }
    
    
  }
  
  return(kcv/k)
}

kcv <- function(ts, c, d, k, b_time, b_timese, mp_type, r = 1, s = 1){
  res = matrix(nrow = 1, ncol = 3)
  for(i in 1:c){
    for(j in 1:d){
      if(i == 1 & b_time == "Cspli"){
        res = rbind(res, c(i, j, "NA"))
      }else if(j == 1 & b_timese == "Cspli"){
        res = rbind(res, c(i, j, "NA"))
      } else{
        cv_score = k_fold(ts, i, j, k, b_time, b_timese, mp_type, r, s)[[1]]
        res = rbind(res, c(i, j, cv_score))
      }
    }
  }
  res = res[-1,]
  return(res)
}
# AIC
AIC_score <- function(ts, c, d, b_time, b_timese, mp_type, r = 1, s = 1){
  n = length(ts)
  esti_df = esti_beta(ts, c, d, b_time, b_timese, mp_type, r, s)

  Y = matrix(ts[(r+1):n], ncol = 1)
  
  W = esti_df[[2]]
  p = dim(W)[2]
  hat = W%*%solve(t(W)%*%W, tol = 1e-40)%*%t(W)
  SSE = t(Y)%*%(diag(n-r)-hat)%*%Y
  #BIC = n*log(SSE/n) + log(n)*p
  AIC =  n*log(SSE/n) + 2*p
  return(AIC)
}

aic <- function(ts, c, d, b_time, b_timese, mp_type, r = 1, s = 1){ # c is the largest 
  res = matrix(nrow = 1, ncol = 3)
  for(i in 1:c){
    for(j in 1:d){
      if(i == 1 & b_time == "Cspli"){
        res = rbind(res, c(i, j, "NA"))
      }else if(j == 1 & b_timese == "Cspli"){
        res = rbind(res, c(i, j, "NA"))
      } else{
        AIC = AIC_score(ts, i, j, b_time, b_timese,mp_type, r, s)
        res = rbind(res, c(i, j, AIC))
      }
    }
  }
  res = res[-1,]
  return(res)
}

# BIC 
BIC_score <- function(ts, c, d, b_time, b_timese, mp_type, r = 1, s = 1){
  n = length(ts)
  esti_df = esti_beta(ts, c, d, b_time, b_timese, mp_type, r, s)
  Y = matrix(ts[(r+1):n], ncol = 1)
  
  W = esti_df[[2]]
  p = dim(W)[2]
  hat = W%*%solve(t(W)%*%W, tol = 1e-40)%*%t(W)
  SSE = t(Y)%*%(diag(n-r)-hat)%*%Y
  BIC = n*log(SSE/n) + log(n)*p
  #AIC =  n*log(SSE/n) + 2*p
  return(BIC)
}

bic <- function(ts, c, d, b_time, b_timese, mp_type, r = 1, s = 1){
  res = matrix(nrow = 1, ncol = 3)
  for(i in 1:c){
    for(j in 1:d){
      if(i == 1 & b_time == "Cspli"){
        res = rbind(res, c(i, j, "NA"))
      }else if(j == 1 & b_timese == "Cspli"){
        res = rbind(res, c(i, j, "NA"))
      } else{
        BIC = BIC_score(ts, i, j, b_time, b_timese, mp_type, r, s)
        res = rbind(res, c(i, j, BIC))
      }
    }
  }
  res = res[-1,]
  return(res)
}





# Simulation for Simultaneous Confidence Region(SCR)  (auto version, automatically choosing c and d), we do not perform 3D plot. 
# fix_t_scr (fix time t)

SCR_fix_t <- function(ts, m, t, c, d, b_time, b_timese, mp_type, r=1, s = 1, n_point = 2000, upper = 10){
  
  esti_df = esti_beta(ts, c, d, b_time, b_timese, mp_type, r, s)
  W = esti_df[[2]]
  n = length(ts)
  beta_hat = esti_df[[1]]
  
  df_basis = 0
  df_basis_d = 0
  if(b_timese == "Cspli"){
    df_basis_d = Cspline_table(d) # true is d - 2 + or 
    d = dim(df_basis_d)[2]
  }
  
  # print(d)
  
  # change for wavelet
  wavelet_basis = c("db1", "db2", "db3", "db4", "db5",
                    "db6", "db7", "db8", "db9", "db10",
                    "db11", "db12", "db13", "db14", "db15",
                    "db16", "db17", "db18", "db19", "db20",
                    "cf1", "cf2", "cf3", "cf4", "cf5"
  )
  
  if(b_timese %in% wavelet_basis){
    df_basis_d = db_table(d, b_timese)
    d = dim(df_basis_d)[2]
  }
  
  if(b_timese == "tri"){
    d_tri = 2*d - 1
  } else{
    d_tri = d
  }
  
  n_esti = n-r
  basis_x = matrix(nrow = n-r, ncol = 1) 
  aux_phi_x = matrix(ncol = d_tri)
  
  basis_t = matrix(nrow = n_esti) 
  aux_bt  = matrix(nrow = n_esti)
  
  b = matrix(nrow = n_esti, ncol = 1)
  
  # print(select_basis_timese(d, mp_selection(mp_type, ts[3], s)[1], b_timese, df_basis))
  
  for(j in 1:r){ # j for lag
    for(i in (r+1):n){
      aux_phi_x = rbind(aux_phi_x, sqrt(mp_selection(mp_type, ts[i-j], s)[2])*select_basis_timese(d, mp_selection(mp_type, ts[i-j], s)[1], b_timese, df_basis_d))
    }
    aux_phi_x = aux_phi_x[-1, ]
    aux_phi_x = as.matrix(aux_phi_x, ncol = d_tri)
    basis_x = cbind(basis_x, aux_phi_x)
    aux_phi_x = matrix(ncol = d_tri)
  }
  
  basis_x = matrix(basis_x[,-1], ncol = r*d_tri)
  
  aux_c = c
  
  if(b_time == "Cspli"){
    df_basis = Cspline_table(aux_c, n) # true is d - 2 + or 
    c = dim(df_basis)[2]
  }
  
  if(b_time %in% wavelet_basis){
    df_basis_1 = wavelet_kth_b(aux_c, n, b_time)
    c = dim(df_basis_1)[2]
  }
  
  if(b_time == "tri"){
    c = aux_c*2 - 1 
  }
  
  # print(c(c, d_tri))
  
  
  
  # for(i in (r+1):n){
  #  basis_x = rbind(basis_x, sqrt(mp_selection(mp_type, ts[i-1], s)[2])*select_basis_timese(d, mp_selection(mp_type, ts[i-1], s)[1], b_timese, df_basis))
  #}
  
 # basis_x = basis_x[-1, ]
 # basis_x = as.matrix(basis_x, ncol = d)
  
  for(k in 1:r){
    for(l1 in 1:c){
      if(b_time == "Cspli"){
        basis_t = as.matrix(df_basis[(r+1):n, l1])
      }else if(b_time %in% wavelet_basis){
        basis_t = as.matrix(df_basis_1[(r+1):n, l1])
      }else{
        basis_t = bs.gene(b_time, l1, n)[(r+1):n, 1]
      }
      aux_bt = cbind(aux_bt, basis_t)
      
      for(l2 in 1:d_tri){
        b = cbind(b, matrix(basis_t*basis_x[, (k-1)*d_tri + l2], ncol = 1))
      }
    }
  }
  
 
  colnames(b) = NULL
  colnames(aux_bt) = NULL
  
  b = as.matrix(b[,-1], ncol = r*c*d_tri)
  aux_bt = as.matrix(aux_bt[,-1], ncol = r*c) 
  aux_bt = as.matrix(aux_bt[ , 1:c], ncol = c)
  # I_r = diag(r)
  # l_j = diag(c*d)  
  # aux.esti = t(t(l_j%*%beta_hat)%*%(t(kronecker(b, I_r))))[,1]
  
  esti = list()
  for(k in 1:r){
    aux_beta_hat = matrix(beta_hat[((k-1)*c*d_tri+1):(k*c*d_tri), 1], ncol = 1)
    # esti[[k]] = t(t(l_j%*%aux_beta_hat)%*%(t(kronecker(b, I_r))))[,1]
    
    aux_b = matrix(b[,((k-1)*c*d_tri + 1):(k*c*d_tri)], ncol = c*d_tri)
    esti[[k]] = (aux_b%*%aux_beta_hat)[,1]
  }
  
  
  
  # inference for SCR 
  
  aux_Sigma_hat = list()
  
  for(k in 1:r){
    aux_W  = matrix(W[,((k-1)*c*d_tri + 1):(k*c*d_tri)], ncol = c*d_tri)
    aux_Sigma_hat[[k]] = 1/n*t(aux_W)%*%aux_W
  }
  
  
  
  n_len = n_point
  t_i = seq(0, 1, length.out = n_len)
  
  # mimic x_i in R  or R+
  
  pos_map = c( "algebp", "logarip") 
  if(mp_type %in% pos_map){
    x_i = seq(0, upper, length.out = n_len)
  }else{
    x_i = seq(-upper, upper, length.out = n_len) 
  }
   
  
  
  basis_xi = matrix(ncol = d_tri)
  basis_ti = matrix(nrow = n_len) 
  
  if(b_time == "Cspli"){
    df_basis = Cspline_table(aux_c, n_len) # true is d - 2 + or 
    c = dim(df_basis)[2]
  }
  
  if(b_time %in% wavelet_basis){
    df_basis_1 = wavelet_kth_b(aux_c, n_len, b_time)
    c = dim(df_basis_1)[2]
  }
  
  #print(c(c, d_tri))
  
  for(l1 in 1:c){
    if(b_time == "Cspli"){
      aux_bti = as.matrix(df_basis[, l1])
    }else if(b_time %in% wavelet_basis){
      aux_bti = as.matrix(df_basis_1[, l1])
    }else{
      aux_bti = bs.gene(b_time, l1, n_len)
    }
    basis_ti = cbind(basis_ti, aux_bti)
  }
  
  basis_ti = basis_ti[,-1]
  basis_ti = as.matrix(basis_ti, ncol = c)
  colnames(basis_ti) = NULL
  
  
  for(i in 1:n_len){
    basis_xi = rbind(basis_xi, sqrt(mp_selection(mp_type, x_i[i], s)[2])*select_basis_timese(d, mp_selection(mp_type, x_i[i], s)[1], b_timese, df_basis_d))
  }
  basis_xi = basis_xi[-1, ]
  basis_xi = as.matrix(basis_xi, ncol = d_tri)
  
  
  
  ## Bootstrap
  tao_k = c()
  h_hat = rep(0, 1000)
  T_1k = list()
  b = matrix(nrow = c*d_tri)
  
  func_t = seq(0, 1, length.out = dim(basis_ti)[1])
  index_t = unique(which(abs(func_t - t) == min(abs(func_t - t))))[1]
  
  
  for(M in 1:n_len){
    ti = matrix(as.numeric(basis_ti[index_t,]), ncol = 1) # we fix time t 
    xj = matrix(as.numeric(basis_xi[M,]), ncol = 1)
    
    b = cbind(b, kronecker(ti, xj))
  }
  b = b[,-1]
  
  aux.esti = rep(0, n_esti)
  for(i in 1:n_esti){
    for(k in 1:r){
      aux.esti[i] = aux.esti[i] + esti[[k]][i] 
    }
  }
  
  if(m == "MV"){
    m = min_vola(ts, aux.esti, basis_x, aux_bt, r)
  }
  
  
  interval = list()
  
  #print(c(c, d_tri))
  
  for(index_r in 1:r){
    for(k in 1:1000){
      Xi_i = matrix(rep(0, c*d_tri), ncol = 1)
      for(i in (r+1):(n-m)){
        R_i = rnorm(1, 0, 1)
        z_j = (ts[i]-aux.esti[i-r])*matrix(basis_x[i-r, ((index_r-1)*d_tri + 1):(index_r*d_tri)],ncol = 1)
        for(j in (i+1):(i+m)){
          z_j = z_j + (ts[j]-aux.esti[j-r])*matrix(basis_x[j-r, ((index_r-1)*d_tri + 1):(index_r*d_tri)], ncol = 1)
        }
        a_ti = matrix(aux_bt[i-r,], ncol =1)
        #print(c(dim(Xi_i), dim(kronecker(z_j, a_ti)) ))
        Xi_i = Xi_i + R_i*kronecker(z_j, a_ti)
      }
      Xi_i = 1/(sqrt((n-m-r)*(m)))*Xi_i
      
      T_1k[[k]] = t(Xi_i)%*%solve(aux_Sigma_hat[[index_r]],  tol = 1e-40)%*%diag(c*d_tri)%*%b
    }
    
    df_T_1k = matrix(unlist(T_1k), ncol = 1000)
    h_hat = apply(df_T_1k, 1, sd) # using function sd()
    
    # df_T_1k = matrix(unlist(T_1k), ncol = 1000)
    for(k in 1:1000){
      tao_k[k] = max(abs(df_T_1k[,k]/h_hat))
    }
    s_tao = sort(tao_k)
    c_alpha = s_tao[1000*(1-0.05)]
    interval[[index_r]] = c_alpha/sqrt(n)*h_hat
  }
  
  return(interval)
}

fix_t_scr <- function(timese, t, c, d, m, b_time, b_timese, mp_type, r = 1, s = 1, n_point = 2000, upper = 10){
  # n_esti = length(timese)
  
  width = SCR_fix_t(timese, m, t, c, d, b_time, b_timese, mp_type, r, s, n_point = n_point, upper = upper)
  esti = fix_t_esti(timese, c, d, t, b_time, b_timese, mp_type, r, s, n_esti = n_point, upper = upper)
  f.df = list()
  
  pos_map = c( "algebp", "logarip") 
  
  if(mp_type %in% pos_map){
    x_i = seq(0, upper, length.out = n_point)
  }else{
    x_i = rep(seq(-upper, upper, length.out = n_point), 3)
  }
  
  for(k in 1:r){
    upper = esti[[k]] + width[[k]]
    lower = esti[[k]] - width[[k]]
    f.df[[k]] = data.frame(x=x_i, y = c(esti[[k]], upper, lower), order = as.factor(rep(c("esti", "upper", "lower"), each = n_point)))
  }
  return(f.df)
}


SCR_fix_x <- function(ts, m, fix_x, c, d, b_time, b_timese, mp_type, r = 1, s = 1, n_point = 2000, upper = 10){

  esti_df = esti_beta(ts, c, d, b_time, b_timese, mp_type, r, s)
  W = esti_df[[2]]
  n = length(ts)
  beta_hat = esti_df[[1]]
  
  df_basis = 0
  df_basis_d = 0
  if(b_timese == "Cspli"){
    df_basis_d = Cspline_table(d) # true is d - 2 + or 
    d = dim(df_basis_d)[2]
  }
  
  # change for wavelet
  wavelet_basis = c("db1", "db2", "db3", "db4", "db5",
                    "db6", "db7", "db8", "db9", "db10",
                    "db11", "db12", "db13", "db14", "db15",
                    "db16", "db17", "db18", "db19", "db20",
                    "cf1", "cf2", "cf3", "cf4", "cf5"
  )  
  if(b_timese %in% wavelet_basis){
    df_basis_d = db_table(d, b_timese)
    d = dim(df_basis_d)[2]
  }
  
  if(b_timese == "tri"){
    d_tri = 2*d - 1
  } else{
    d_tri = d
  }
  
  n_esti = n-r
  # basis_x = matrix(ncol = d) 
  basis_x = matrix(nrow = n-r, ncol = 1) 
  aux_phi_x = matrix(ncol = d_tri)
  
  basis_t = matrix(nrow = n_esti) 
  aux_bt  = matrix(nrow = n_esti)
  
  b = matrix(nrow = n_esti, ncol = 1)
  
 
  for(j in 1:r){ # j for lag
    for(i in (r+1):n){
      aux_phi_x = rbind(aux_phi_x, sqrt(mp_selection(mp_type, ts[i-j], s)[2])*select_basis_timese(d, mp_selection(mp_type, ts[i-j], s)[1], b_timese, df_basis_d))
    }
    aux_phi_x = aux_phi_x[-1, ]
    aux_phi_x = as.matrix(aux_phi_x, ncol = d_tri)
    basis_x = cbind(basis_x, aux_phi_x)
    aux_phi_x = matrix(ncol = d_tri)
  }
  
  basis_x = matrix(basis_x[,-1], ncol = r*d_tri)
  
  #for(i in (r+1):n){
  #  basis_x = rbind(basis_x, sqrt(mp_selection(mp_type, ts[i-1], s)[2])*select_basis_timese(d, mp_selection(mp_type, ts[i-1], s)[1], b_timese, df_basis))
  #}
 # basis_x = basis_x[-1, ]
  #basis_x = as.matrix(basis_x, ncol = d)
  
  aux_c = c
  
  if(b_time == "Cspli"){
    df_basis = Cspline_table(aux_c, n) # true is d - 2 + or 
    c = dim(df_basis)[2]
  }
  
  if(b_time %in% wavelet_basis){
    df_basis_1 = wavelet_kth_b(aux_c, n, b_time)
    c = dim(df_basis_1)[2]
  }
  
  if(b_time == "tri"){
    c = c*2 - 1 
  }
  
  for(k in 1:r){
    for(l1 in 1:c){
      if(b_time == "Cspli"){
        basis_t = as.matrix(df_basis[(r+1):n, l1])
      }else if(b_time %in% wavelet_basis){
        basis_t = as.matrix(df_basis_1[(r+1):n, l1])
      }else{
        basis_t = bs.gene(b_time, l1, n)[(r+1):n, 1]
      }
      aux_bt = cbind(aux_bt, basis_t)
      for(l2 in 1:d_tri){
        b = cbind(b, matrix(basis_t*basis_x[, (k-1)*d_tri + l2], ncol = 1))
      }
    }
  }
  
  #for(l1 in 1:c){
  #  if(b_time == "Cspli"){
  #    basis_t = as.matrix(df_basis[(r+1):n, l1])
  #  }else if(b_time %in% wavelet_basis){
  #    basis_t = as.matrix(df_basis_1[(r+1):n, l1])
  #  }else{
  #    basis_t = bs.gene(b_time, l1, n)[(r+1):n, 1]
  #  }
  #  aux_bt = cbind(aux_bt, basis_t)
  #  for(l2 in 1:d){
  #    b = cbind(b, matrix(basis_t*basis_x[, l2], ncol = 1))
  #  }
  # }
  
  colnames(b) = NULL
  colnames(aux_bt) = NULL
  
  b = as.matrix(b[,-1], ncol = r*c*d_tri)
  aux_bt = as.matrix(aux_bt[,-1], ncol = r*c) 
  aux_bt = as.matrix(aux_bt[ ,1:c], ncol = c)
  
  #I_r = diag(r)
  #l_j = diag(c*d)  
  #aux.esti = t(t(l_j%*%beta_hat)%*%(t(kronecker(b, I_r))))[,1]
  
  esti = list()
  for(k in 1:r){
    aux_beta_hat = matrix(beta_hat[((k-1)*c*d_tri + 1):(k*c*d_tri), 1], ncol = 1)
    # esti[[k]] = t(t(l_j%*%aux_beta_hat)%*%(t(kronecker(b, I_r))))[,1]
    
    aux_b = matrix(b[,((k-1)*c*d_tri + 1):(k*c*d_tri)], ncol = c*d_tri)
    esti[[k]] = (aux_b%*%aux_beta_hat)[,1]
  }
  
  
  # inference for SCR 
  
  # Sigma_hat = 1/n*t(W)%*%W
  aux_Sigma_hat = list()
  
  for(k in 1:r){
    aux_W  = matrix(W[,((k-1)*c*d_tri + 1):(k*c*d_tri)], ncol = c*d_tri)
    aux_Sigma_hat[[k]] = 1/n*t(aux_W)%*%aux_W
  }
  
  n_len = n_point
  t_i = seq(0, 1, length.out = n_len)
  
  # mimic x_i in R or R+
  pos_map = c( "algebp", "logarip") 
  if(mp_type %in% pos_map){
    x_i = seq(0, upper, length.out = n_len)
  }else{
    x_i = seq(-upper, upper, length.out = n_len) 
  }
  
  basis_xi = matrix(ncol = d_tri) 
  basis_ti = matrix(nrow = n_len) 
  
  
  if(b_time == "Cspli"){
    df_basis = Cspline_table(aux_c, n_len) # true is d - 2 + or 
    c = dim(df_basis)[2]
  }
  
  if(b_time %in% wavelet_basis){
    df_basis_1 = wavelet_kth_b(aux_c, n_len, b_time)
    c = dim(df_basis_1)[2]
  }
  
  
  for(l1 in 1:c){
    if(b_time == "Cspli"){
      aux_bti = as.matrix(df_basis[, l1])
    }else if(b_time %in% wavelet_basis){
      aux_bti = as.matrix(df_basis_1[, l1])
    }else{
      aux_bti = bs.gene(b_time, l1, n_len)
    }
    basis_ti = cbind(basis_ti, aux_bti)
  }
  
  basis_ti = basis_ti[,-1]
  basis_ti = as.matrix(basis_ti, ncol = c)
  colnames(basis_ti) = NULL
  
  
  
  for(i in 1:n_len){
    basis_xi = rbind(basis_xi, sqrt(mp_selection(mp_type, x_i[i], s)[2])*select_basis_timese(d, mp_selection(mp_type, x_i[i], s)[1], b_timese, df_basis_d))
  }
  basis_xi = basis_xi[-1, ]
  basis_xi = as.matrix(basis_xi, ncol = d_tri)
  
  
  index_x = unique(which(abs(x_i-fix_x) == min(abs(x_i-fix_x))))[1]
  
  ## Bootstrap
  tao_k = c()
  h_hat = rep(0, 1000)
  T_1k = list()
  b = matrix(nrow = c*d_tri)
  
  for(M in 1:n_len){
    ti = matrix(as.numeric(basis_ti[M,]), ncol = 1) # we fix time t 
    xj = matrix(as.numeric(basis_xi[index_x,]), ncol = 1)
    
    b = cbind(b, kronecker(ti, xj))
  }
  b = b[,-1]

  
  aux.esti = rep(0, n_esti)
  for(i in 1:n_esti){
    for(k in 1:r){
      aux.esti[i] = aux.esti[i] + esti[[k]][i] 
    }
  }
  
  if(m == "MV"){
    m = min_vola(ts, aux.esti, basis_x, aux_bt, r)
  }
  
  interval = list()

  for(index_r in 1:r){
    for(k in 1:1000){
      Xi_i = matrix(rep(0, c*d_tri), ncol = 1)
      for(i in (r+1):(n-m)){
        R_i = rnorm(1, 0, 1)
        z_j = (ts[i]-aux.esti[i-r])*matrix(basis_x[i-r, ((index_r-1)*d_tri + 1):(index_r*d_tri)],ncol = 1)
        for(j in (i+1):(i+m)){
          z_j = z_j + (ts[j]-aux.esti[j-r])*matrix(basis_x[j-r, ((index_r-1)*d_tri + 1):(index_r*d_tri)], ncol = 1)
        }
        a_ti = matrix(aux_bt[i-r,], ncol =1)
        
        #print(dim(R_i*kronecker(z_j, a_ti)))
        Xi_i = Xi_i + R_i*kronecker(z_j, a_ti)
      }
      Xi_i = 1/(sqrt((n-m-r)*(m)))*Xi_i
      # print(dim(Xi_i))
      # print(dim(Sigma_hat))
      T_1k[[k]] = t(Xi_i)%*%solve(aux_Sigma_hat[[index_r]],  tol = 1e-40)%*%diag(c*d_tri)%*%b
    }
    
    df_T_1k = matrix(unlist(T_1k), ncol = 1000)
    h_hat = apply(df_T_1k, 1, sd)
    
    # df_T_1k = matrix(unlist(T_1k), ncol = 1000)
    for(k in 1:1000){
      tao_k[k] = max(abs(df_T_1k[,k]/h_hat))
    }
    s_tao = sort(tao_k)
    c_alpha = s_tao[1000*(1-0.05)]
    interval[[index_r]] = c_alpha/sqrt(n)*h_hat
    
  }
  
  return(interval)
}

fix_x_scr <- function(timese, fix_x, c, d, m, b_time, b_timese, mp_type, r=1, s=1, n_point = 2000, upper = 10){
  #n_esti = length(timese)
  width = SCR_fix_x(timese, m, fix_x, c, d, b_time, b_timese, mp_type, r, s, n_point = n_point, upper = upper)
  esti = fix_x_esti(timese, c, d, fix_x, b_time, b_timese, mp_type, r, s, n_esti = n_point, upper = upper)
  f.df = list()
  for(k in 1:r){
    upper = esti[[k]] + width[[k]]
    lower = esti[[k]] - width[[k]]
    f.df[[k]] = data.frame(t=rep(seq(0,1, length.out = n_point),3), y = c(esti[[k]], upper, lower), order = as.factor(rep(c("esti", "upper", "lower"), each = n_point)))
  }
  return(f.df)
}



#Inference (for each test, we have auto version)

# Testing time-homogeneity
homot_scr <- function(timese, t, c, d, m, b_time, b_timese, mp_type, r = 1, s = 1, n_point = 2000,  upper = 10){
  width = SCR_fix_t(timese, m, t, c, d, b_time, b_timese, mp_type, r, s, n_point = n_point, upper = upper)
  n_esti = n_point # length(timese)
  beta_hat = esti_beta(timese,c, d, b_time, b_timese, mp_type, r, s)[[1]]
  
  df_basis = 0
  df_basis_d = 0
  if(b_timese == "Cspli"){
    df_basis_d = Cspline_table(d) # true is d - 2 + or 
    d = dim(df_basis_d)[2]
  }
  
  # change for wavelet
  wavelet_basis = c("db1", "db2", "db3", "db4", "db5",
                    "db6", "db7", "db8", "db9", "db10",
                    "db11", "db12", "db13", "db14", "db15",
                    "db16", "db17", "db18", "db19", "db20",
                    "cf1", "cf2", "cf3", "cf4", "cf5"
  )  
  if(b_timese %in% wavelet_basis){
    df_basis_d = db_table(d, b_timese)
    d = dim(df_basis_d)[2]
  }
  
  if(b_timese == "tri"){
    d_tri = 2*d - 1
  } else{
    d_tri = d
  }
  
  basis_x = matrix(ncol = d_tri) 
  
  # basis_x = matrix(nrow = n_esti, ncol = 1) 
  # aux_phi_x = matrix(ncol = d)

  b = matrix(nrow = n_esti, ncol = 1)
  b_homo = matrix(nrow = n_esti, ncol = 1)
  
  # mimic x_i in R or R+
  pos_map = c( "algebp", "logarip") 
  if(mp_type %in% pos_map){
    x = seq(0, upper, length.out = n_esti)
  }else{
    x = seq(-upper, upper, length.out = n_esti) 
  }
  
  
  func_t = seq(0, 1, length.out = n_esti)
  index_t = unique(which(abs(func_t - t) == min(abs(func_t - t))))[1]
  
  #for(j in 1:r){ # j for lag
  #  for(i in (r+1):n){
  #    aux_phi_x = rbind(aux_phi_x, sqrt(mp_selection(mp_type, timese[i-j], s)[2])*select_basis_timese(d, mp_selection(mp_type, timese[i-j], s)[1], b_timese, df_basis))
  #  }
  #  aux_phi_x = aux_phi_x[-1, ]
  #  aux_phi_x = as.matrix(aux_phi_x, ncol = d)
  #  basis_x = cbind(basis_x, aux_phi_x)
  #  aux_phi_x = matrix(ncol = d)
 # }
  
 # basis_x = matrix(basis_x[,-1], ncol = r*d)
  
  for(i in 1:n_esti){
    basis_x = rbind(basis_x, sqrt(mp_selection(mp_type, x[i], s)[2])*select_basis_timese(d, mp_selection(mp_type, x[i], s)[1], b_timese, df_basis_d))
  }
  
  basis_x = basis_x[-1, ]
  # print(dim(basis_x))
  basis_x = as.matrix(basis_x, ncol = d_tri)
  
  
  if(b_time == "Cspli"){
    df_basis = Cspline_table(c, n_esti) # true is d - 2 + or 
    c = dim(df_basis)[2]
  }
  
  if(b_time %in% wavelet_basis){
    df_basis_1 = wavelet_kth_b(c, n_esti, b_time)
    c = dim(df_basis_1)[2]
  }
  
  if(b_time == "tri"){
    c = c*2 -1 
  }
  
  for(l1 in 1:c){
    if(b_time == "Cspli"){
      basis_t = as.matrix(df_basis[, l1])
      }else if(b_time %in% wavelet_basis){
        basis_t = as.matrix(df_basis_1[, l1])
      }else{
        basis_t = as.matrix(bs.gene(b_time, l1, n_esti))
      }
    for(l2 in 1:d_tri){
      b = cbind(b, basis_t[index_t]*matrix(basis_x[, l2], ncol = 1))
      b_homo = cbind(b_homo, matrix(basis_x[, l2], ncol = 1))
    }
    }

  colnames(b) = NULL
  colnames(b_homo) = NULL
  # I_r = diag(r)
  b = as.matrix(b[,-1], ncol = c*d_tri) 
  b_homo = as.matrix(b_homo[,-1], ncol = c*d_tri) 
  
  # l_j = diag(c*d)  
  
  
  esti = list()
  for(k in 1:r){
    aux_beta_hat = matrix(beta_hat[((k-1)*c*d_tri+1):(k*c*d_tri), 1], ncol = 1)
    # esti[[k]] = t(t(l_j%*%aux_beta_hat)%*%(t(kronecker(b, I_r))))[,1]
    #print(c(dim(b), dim(aux_beta_hat)))
    # aux_b = matrix(b[,((k-1)*c*d + 1):(k*c*d)], ncol = c*d)
    esti[[k]] = (b%*%aux_beta_hat)[,1]
  }
  
  upper_b = list()
  lower = list()
  # esti = t(t(l_j%*%beta_hat)%*%(t(kronecker(b, I_r))))[,1]
  for(k in 1:r){
    upper_b[[k]] = esti[[k]] + width[[k]]
    lower[[k]] = esti[[k]] - width[[k]]
  }
  
  esti = list()
  pos_map = c( "algebp", "logarip") 
  for(k in 1:r){
    aux_beta_hat = matrix(beta_hat[((k-1)*c*d_tri+1):(k*c*d_tri), 1], ncol = 1)
    # esti[[k]] = t(t(l_j%*%aux_beta_hat)%*%(t(kronecker(b, I_r))))[,1]
    
    # aux_b = matrix(b_homo[,((k-1)*c*d + 1):(k*c*d)], ncol = c*d)
    esti[[k]] = (b_homo%*%aux_beta_hat)[,1]
  }
  # esti = t(t(l_j%*%beta_hat)%*%(t(kronecker(b_homo, I_r))))[,1]
  f.df = list()
  
  if(mp_type %in% pos_map){
    x_temp = rep(seq(0, upper, length.out = n_esti),3)
  }else{
    x_temp = rep(seq(-upper, upper, length.out = n_esti),3)
  }
  
  for(k in 1:r){
    f.df[[k]] = data.frame(x = x_temp, y = c(esti[[k]], upper_b[[k]], lower[[k]]), order = as.factor(rep(c("esti", "upper", "lower"), each = n_esti)))
  }
  return(f.df)
}

# Testing separability(product) we need make sure that integration is not 0. 


# fix t

fix_t_sepera_scr <- function(timese, t, c, d, m, b_time, b_timese, mp_type, r = 1, s = 1, n_point = 2000,  upper = 10){
  
  n_esti = n_point # length(timese)
  beta_hat = esti_beta(timese, c, d, b_time, b_timese, mp_type, r, s)[[1]]
  m_hat_ij = matrix(rep(0, n_esti*n_esti), nrow = n_esti, ncol = n_esti)
  width = SCR_fix_t(timese, m, t, c, d, b_time, b_timese, mp_type, r = r, s = s, n_point = n_point,  upper = upper)
  
  
  df_basis = 0
  if(b_timese == "Cspli"){
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
  if(b_timese %in% wavelet_basis){
    df_basis = db_table(d, b_timese)
    d = dim(df_basis)[2]
  }
  
  if(b_timese == "tri"){
    d_tri = 2*d - 1
  } else{
    d_tri = d
  }
  
  m_sep_ij = c() # fix t
  basis_xi = matrix(ncol = d_tri) 
  basis_ti = matrix(nrow = n_esti) 
  
  
  pos_map = c( "algebp", "logarip") 
  if(mp_type %in% pos_map){
    x_i = seq(0, upper, length.out = n_esti)
  }else{
    x_i = seq(-upper, upper, length.out = n_esti) 
  }

  for(i in 1:n_esti){
    basis_xi = rbind(basis_xi, sqrt(mp_selection(mp_type, x_i[i], s)[2])*select_basis_timese(d, mp_selection(mp_type, x_i[i], s)[1], b_timese, df_basis))
  }
  basis_xi = basis_xi[-1, ]
  basis_xi = as.matrix(basis_xi, ncol = d_tri)
  
  
  if(b_time == "Cspli"){
    df_basis = Cspline_table(c, n_esti) # true is d - 2 + or 
    c = dim(df_basis)[2]
  }
  
  if(b_time %in% wavelet_basis){
    df_basis_1 = wavelet_kth_b(c, n_esti, b_time)
    c = dim(df_basis_1)[2]
  }
  
  if(b_time == "tri"){
    c = c*2 -1 
  }
  
  for(l1 in 1:c){
    if(b_time == "Cspli"){
      aux_bti = as.matrix(df_basis[, l1])
    }else if(b_time %in% wavelet_basis){
      aux_bti = as.matrix(df_basis_1[, l1])
    }else{
      aux_bti = bs.gene(b_time, l1, n_esti)
    }
    basis_ti = cbind(basis_ti, aux_bti)
  }
  colnames(basis_ti) = NULL
  basis_ti =  as.matrix(basis_ti[,-1], ncol = c) 
  
  
  # m_hat(i,j)
  esti_sep = list()
  for(k in 1:r){
    for(i in 1:n_esti){
      for(j in 1:n_esti){
        ti = matrix(as.numeric(basis_ti[i,]), ncol = 1)
        xj = matrix(as.numeric(basis_xi[j,]), ncol = 1)
        m_hat_ij[i,j] = t(matrix(beta_hat[((k-1)*c*d+1):(k*c*d),], ncol = 1))%*%kronecker(ti, xj)
      }
    }
    
    
    func_t = seq(0, 1, length.out = n_esti)
    index_t = unique(which(abs(func_t - t) == min(abs(func_t - t))))[1]
    
    f_t = (apply(m_hat_ij,1, sum)/n_esti)*2*upper 
    g_x = (apply(m_hat_ij,2, sum)/n_esti) 
    const = sum((apply(m_hat_ij,1, sum)/n_esti)*2*upper)/n_esti
    
    esti_sep[[k]] = (g_x)*f_t[index_t] / const # there are 2 seperates
  }
  
  # basis_x = matrix(ncol = d_tri) 
  basis_t = matrix(nrow = n_esti) 
  b = matrix(nrow = n_esti, ncol = 1)
  # x = seq(-10, 10, length.out = n_esti)
  
  
  # df_basis = 0
  #if(b_timese == "Cspli"){
  #  df_basis = Cspline_table(d) # true is d - 2 + or 
  #  d = dim(df_basis)[2]
  #}
  
  # for(i in 1:n_esti){
  #  basis_x = rbind(basis_x, sqrt(mp_selection(mp_type, x[i], s)[2])*select_basis_timese(d, mp_selection(mp_type, x[i], s)[1], b_timese, df_basis))
  #}
  #basis_x = basis_x[-1, ]
  
  #if(b_time == "Cspli"){
  #  df_basis = Cspline_table(c, n_esti) # true is d - 2 + or 
  #  c = dim(df_basis)[2]
  #}
  
  for(l1 in 1:c){
    if(b_time == "Cspli"){
      basis_t = as.matrix(df_basis[, l1])
    }else if(b_time %in% wavelet_basis){
      basis_t = as.matrix(df_basis_1[, l1])
    }else{
      basis_t = as.matrix(bs.gene(b_time, l1, n_esti))
    }
    for(l2 in 1:d_tri){
      b = cbind(b, basis_t[index_t]*matrix(basis_xi[, l2], ncol = 1))
    }
  }
  colnames(b) = NULL
  # I_r = diag(r)
  b = as.matrix(b[,-1], ncol = c*d_tri) 
  # l_j = diag(c*d)  
  # esti = t(t(l_j%*%beta_hat)%*%(t(kronecker(b, I_r))))[,1]
  
  esti = list()
  pos_map = c( "algebp", "logarip") 
  for(k in 1:r){
    aux_beta_hat = matrix(beta_hat[((k-1)*c*d_tri+1):(k*c*d_tri), 1], ncol = 1)
    # esti[[k]] = t(t(l_j%*%aux_beta_hat)%*%(t(kronecker(b, I_r))))[,1]
    
    # aux_b = matrix(b[,((k-1)*c*d + 1):(k*c*d)], ncol = c*d)
    esti[[k]] = b%*%aux_beta_hat[,1]
  }
  
  if(mp_type %in% pos_map){
    x_temp = rep(seq(0, upper, length.out = n_esti),3)
  }else{
    x_temp = rep(seq(-upper, upper, length.out = n_esti),3)
  }
  
  f.df = list()
  for(k in 1:r){
    upper_b = esti[[k]] + width[[k]]
    lower = esti[[k]] - width[[k]]
    
    f.df[[k]] = data.frame(x = x_temp, y = c(esti_sep[[k]], upper_b, lower), order = as.factor(rep(c("esti_sep", "upper", "lower"), each = n_esti)))
  }
  
  return(f.df)
}

# fix time series, f_x indicates the index of x_i 
fix_x_sepera_scr <- function(timese, f_x, c, d, m, b_time, b_timese, mp_type, r=1, s=1, n_point = 2000, upper = 10){
  
  n_esti = n_point # length(timese)
  beta_hat = esti_beta(timese, c, d, b_time, b_timese, mp_type, r, s)[[1]]
  m_hat_ij = matrix(rep(0, n_esti*n_esti), nrow = n_esti, ncol = n_esti)
  width = SCR_fix_x(timese, m, f_x, c, d, b_time, b_timese, mp_type, r, s, n_point = n_point)
  
  df_basis = 0
  if(b_timese == "Cspli"){
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
  if(b_timese %in% wavelet_basis){
    df_basis = db_table(d, b_timese)
    d = dim(df_basis)[2]
  }
  
  if(b_timese == "tri"){
    d_tri = 2*d - 1
  } else{
    d_tri = d
  }
  
  m_sep_ij = c() # fix t
  basis_xi = matrix(ncol = d_tri) 
  basis_ti = matrix(nrow = n_esti) 
  
  pos_map = c( "algebp", "logarip") 
  if(mp_type %in% pos_map){
    x_i = seq(0, upper, length.out = n_esti)
  }else{
    x_i = seq(-upper, upper, length.out = n_esti) 
  }
  
  
  for(i in 1:n_esti){
    basis_xi = rbind(basis_xi, sqrt(mp_selection(mp_type, x_i[i], s)[2])*select_basis_timese(d, mp_selection(mp_type, x_i[i], s)[1], b_timese, df_basis))
  }
  basis_xi = basis_xi[-1, ]
  basis_xi = as.matrix(basis_xi, ncol = d_tri)
  
  
  if(b_time == "Cspli"){
    df_basis = Cspline_table(c, n_esti) # true is d - 2 + or 
    c = dim(df_basis)[2]
  }
  
  if(b_time %in% wavelet_basis){
    df_basis_1 = wavelet_kth_b(c, n_esti, b_time)
    c = dim(df_basis_1)[2]
  }
  
  if(b_time == "tri"){
    c = c*2 -1 
  }
  
  for(l1 in 1:c){
    if(b_time == "Cspli"){
      aux_bti = as.matrix(df_basis[, l1])
    }else if(b_time %in% wavelet_basis){
      aux_bti = as.matrix(df_basis_1[, l1])
    }else{
      aux_bti = bs.gene(b_time, l1, n_esti)
    }
    basis_ti = cbind(basis_ti, aux_bti)
  }
  
  colnames(basis_ti) = NULL
  basis_ti =  as.matrix(basis_ti[,-1], ncol = c) 
  
  
  # m_hat(i,j)
  esti_sep = list()
  for(k in 1:r){
    for(i in 1:n_esti){
      for(j in 1:n_esti){
        ti = matrix(as.numeric(basis_ti[i,]), ncol = 1)
        xj = matrix(as.numeric(basis_xi[j,]), ncol = 1)
        m_hat_ij[i,j] = t(matrix(beta_hat[((k-1)*c*d+1):(k*c*d_tri),], ncol = 1))%*%kronecker(ti, xj)
      }
    }
    
    index_x = unique(which(abs(x_i-f_x) == min(abs(x_i-f_x))))[1]
    
    f_t = (apply(m_hat_ij,1, sum)/n_esti)*2*upper
    g_x = (apply(m_hat_ij,2, sum)/n_esti)
    const = sum((apply(m_hat_ij,1, sum)/n_esti)*2*upper)/n_esti
    
    esti_sep[[k]] = (g_x[index_x]*(f_t)) / const
    
  }
  
  

  # basis_x = matrix(ncol = d) 
  basis_t = matrix(nrow = n_esti) 
  b = matrix(nrow = n_esti, ncol = 1)
  # x = seq(-10, 10, length.out = n_esti)
  
  # df_basis = 0
  # if(b_timese == "Cspli"){
  #  df_basis = Cspline_table(d) # true is d - 2 + or 
  #  d = dim(df_basis)[2]
  # }
  
  #for(i in 1:n_esti){
  #  basis_x = rbind(basis_x, sqrt(mp_selection(mp_type, x[i], s)[2])*select_basis_timese(d, mp_selection(mp_type, x[i], s)[1], b_timese, df_basis))
  #}
  #basis_x = basis_x[-1, ]
  
  
  #if(b_time == "Cspli"){
  #  df_basis = Cspline_table(c, n_esti) # true is d - 2 + or 
  #  c = dim(df_basis)[2]
  #}
  
  for(l1 in 1:c){
    if(b_time == "Cspli"){
      basis_t = as.matrix(df_basis[, l1])
    }else if(b_time %in% wavelet_basis){
      basis_t = as.matrix(df_basis_1[, l1])
    }else{
      basis_t = as.matrix(bs.gene(b_time, l1, n_esti))
    }
    for(l2 in 1:d_tri){
      b = cbind(b, basis_t*basis_xi[index_x, l2])
    }
  }
  colnames(b) = NULL
  b = as.matrix(b[,-1], ncol = c*d_tri) 
  
  # I_r = diag(r)
  # l_j = diag(c*d) 
  
  # esti = t(t(l_j%*%beta_hat)%*%(t(kronecker(b, I_r))))[,1]
  esti = list()
  for(k in 1:r){
    aux_beta_hat = matrix(beta_hat[((k-1)*c*d_tri+1):(k*c*d_tri), 1], ncol = 1)
    # esti[[k]] = t(t(l_j%*%aux_beta_hat)%*%(t(kronecker(b, I_r))))[,1]
    
    # aux_b = matrix(b[,((k-1)*c*d + 1):(k*c*d)], ncol = c*d)
    esti[[k]] = b%*%aux_beta_hat[,1]
  }
  
  f.df = list()
  for(k in 1:r){
    upper = esti[[k]] + width[[k]]
    lower = esti[[k]] - width[[k]]
    f.df[[k]] = data.frame(t = rep(seq(0, 1, length.out = n_esti), 3), y = c(esti_sep[[k]], upper, lower), order = as.factor(rep(c("esti_sep", "upper", "lower"), each = n_esti)))
    
  }
  
  return(f.df)
}


sepera_3d <- function(timese, c, d, m, b_time, b_timese, mp_type, r = 1, s = 1, n_point = 2000,  upper = 10){
  
  n_esti = n_point # length(timese)
  beta_hat = esti_beta(timese, c, d, b_time, b_timese, mp_type, r, s)[[1]]
  m_hat_ij = matrix(rep(0, n_esti*n_esti), nrow = n_esti, ncol = n_esti)
  aux_t_basis = seq(0, 1, length.out = n_point)
  aux_x_basis = seq(-10, 10, length.out = n_point)
  
  scr.df = list()
  for(j in 1:n_point){
    scr.df[[j]] = fix_x_scr(timese, aux_x_basis[j], c, d, m, b_time, b_timese, mp_type, r = r, s = s, n_point = n_point,  upper = upper)
    }
  
  
  df_basis = 0
  if(b_timese == "Cspli"){
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
  if(b_timese %in% wavelet_basis){
    df_basis = db_table(d, b_timese)
    d = dim(df_basis)[2]
  }
  
  if(b_timese == "tri"){
    d_tri = 2*d - 1
  } else{
    d_tri = d
  }
  
  m_sep_ij = c() # fix t
  basis_xi = matrix(ncol = d_tri) 
  basis_ti = matrix(nrow = n_esti) 
  
  
  pos_map = c( "algebp", "logarip") 
  if(mp_type %in% pos_map){
    x_i = seq(0, upper, length.out = n_esti)
  }else{
    x_i = seq(-upper, upper, length.out = n_esti) 
  }
  
  for(i in 1:n_esti){
    basis_xi = rbind(basis_xi, sqrt(mp_selection(mp_type, x_i[i], s)[2])*select_basis_timese(d, mp_selection(mp_type, x_i[i], s)[1], b_timese, df_basis))
  }
  basis_xi = basis_xi[-1, ]
  basis_xi = as.matrix(basis_xi, ncol = d_tri)
  
  
  if(b_time == "Cspli"){
    df_basis = Cspline_table(c, n_esti) # true is d - 2 + or 
    c = dim(df_basis)[2]
  }
  
  if(b_time %in% wavelet_basis){
    df_basis_1 = wavelet_kth_b(c, n_esti, b_time)
    c = dim(df_basis_1)[2]
  }
  
  if(b_time == "tri"){
    c = c*2 -1 
  }
  
  for(l1 in 1:c){
    if(b_time == "Cspli"){
      aux_bti = as.matrix(df_basis[, l1])
    }else if(b_time %in% wavelet_basis){
      aux_bti = as.matrix(df_basis_1[, l1])
    }else{
      aux_bti = bs.gene(b_time, l1, n_esti)
    }
    basis_ti = cbind(basis_ti, aux_bti)
  }
  colnames(basis_ti) = NULL
  basis_ti =  as.matrix(basis_ti[,-1], ncol = c) 
  
  
  sep_m_hat_ij = matrix(rep(0, n_esti*n_esti), nrow = n_esti, ncol = n_esti)
    
  # m_hat(i,j)
  esti_sep = list()
  for(k in 1:r){
    for(i in 1:n_esti){
      for(j in 1:n_esti){
        ti = matrix(as.numeric(basis_ti[i,]), ncol = 1)
        xj = matrix(as.numeric(basis_xi[j,]), ncol = 1)
        m_hat_ij[i,j] = t(matrix(beta_hat[((k-1)*c*d+1):(k*c*d),], ncol = 1))%*%kronecker(ti, xj)
      }
    }
    
    
    # func_t = seq(0, 1, length.out = n_esti)
    #index_t = unique(which(abs(func_t - t) == min(abs(func_t - t))))[1]
    
    f_t = (apply(m_hat_ij,1, sum)/n_esti)*2*upper 
    g_x = (apply(m_hat_ij,2, sum)/n_esti) 
    const = sum((apply(m_hat_ij,1, sum)/n_esti)*2*upper)/n_esti
    
    for(l in 1:n_esti){
      for(m in 1:n_esti){
        sep_m_hat_ij[l, m] = f_t[l]*g_x[m] / const
      }
    }
    esti_sep[[k]] = sep_m_hat_ij # there are 2 seperates
  }
  
  
  return(list(scr.df, esti_sep = esti_sep))
}


# Exact form Test # (input should be one data frame the exact function) 

### should create function basis on original code. 

exact_form_test <- function(ts, c, d, m, b_time, b_timese, mp_type, exact_func, r = 1, s = 1, upper = 10){ # exact_func is list which contain matrix contains m_k_0(i,j) in (i,j)-th elements, k = 1,2,...,r
  
  n = length(ts)
  n_esti =  dim(exact_func[[1]])[1]
  
  basis_ti = matrix(nrow = n_esti) 
  
  pos_map = c( "algebp", "logarip") 
  if(mp_type %in% pos_map){
    x_i = seq(0, upper, length.out = n_esti)
  }else{
    x_i = seq(-upper, upper, length.out = n_esti) 
  }
  
  esti_df = esti_beta(ts, c, d, b_time, b_timese, mp_type, r, s)
  beta_hat = esti_df[[1]]
  W =  esti_df[[2]]
  
  df_basis_d = 0
  df_basis = 0
  if(b_timese == "Cspli"){
    df_basis_d = Cspline_table(d) # true is d - 2 + or 
    d = dim(df_basis_d)[2]
  }
  
  # change for wavelet
  wavelet_basis = c("db1", "db2", "db3", "db4", "db5",
                    "db6", "db7", "db8", "db9", "db10",
                    "db11", "db12", "db13", "db14", "db15",
                    "db16", "db17", "db18", "db19", "db20",
                    "cf1", "cf2", "cf3", "cf4", "cf5"
  )  
  if(b_timese %in% wavelet_basis){
    df_basis_d = db_table(d, b_timese)
    d = dim(df_basis_d)[2]
  }
  
  if(b_timese == "tri"){
    d_tri = 2*d - 1
  } else{
    d_tri = d
  }
  
  basis_xi = matrix(ncol = d_tri) 
  
  for(i in 1:n_esti){
    basis_xi = rbind(basis_xi, sqrt(mp_selection(mp_type, x_i[i], s)[2])*select_basis_timese(d, mp_selection(mp_type, x_i[i], s)[1], b_timese, df_basis_d))
  }
  basis_xi = basis_xi[-1, ]
  basis_xi = as.matrix(basis_xi, ncol = d_tri)
  
  aux_c = c
  if(b_time == "Cspli"){
    df_basis = Cspline_table(aux_c, n_esti) # true is d - 2 + or 
    c = dim(df_basis)[2]
  }
  
  if(b_time %in% wavelet_basis){
    df_basis_1 = wavelet_kth_b(aux_c, n_esti, b_time)
    c = dim(df_basis_1)[2]
  }
  
  if(b_time == "tri"){
    c = c*2 -1 
  }
  
  for(l1 in 1:c){
    if(b_time == "Cspli"){
      aux_bti = as.matrix(df_basis[, l1])
    }else if(b_time %in% wavelet_basis){
      aux_bti = as.matrix(df_basis_1[, l1])
    }else{
      aux_bti = bs.gene(b_time, l1, n_esti)
    }
    basis_ti = cbind(basis_ti, aux_bti)
  }
  
  colnames(basis_ti) = NULL
  basis_ti =  as.matrix(basis_ti[,-1], ncol = c) 
  
  
  # 2 dimensinal integrate
  inte = 0
  t = matrix(seq(0, 1, length.out = n_esti))
  
  if(mp_type %in% pos_map){
    x = seq(0, upper, length.out = n_esti)
  }else{
    x = seq(-upper, upper, length.out = n_esti) 
  }
  
  #x = matrix(seq(-10, 10,length.out = n_esti))
  # m_hat_ij = matrix(rep(0, (n_esti)^2), nrow = n_esti, ncol = n_esti)
  aux_inte = list()
  aux_B_j = list()
  B_j = matrix(rep(0, (c*d_tri)^2), nrow = c*d_tri, ncol = c*d_tri)
  for(k in 1:r){
    aux_exact_func = exact_func[[k]]
    for(i in 1:n_esti){
      for(j in 1:n_esti){
        m_0 = aux_exact_func[i, j]
        ti = matrix(as.numeric(basis_ti[i,]), ncol = 1)
        xj = matrix(as.numeric(basis_xi[j,]), ncol = 1)
        B_j = B_j + kronecker(ti, xj)%*%t(kronecker(ti, xj))
        aux_m_hat_ij = t(matrix(beta_hat[((k-1)*c*d_tri+1):(k*c*d_tri),], ncol = 1))%*%kronecker(ti, xj)
        inte = inte + (aux_m_hat_ij  - m_0)^2
      }
    }
    aux_inte[[k]] = inte
    aux_B_j[[k]] = B_j
    inte = 0
    B_j = matrix(rep(0, (c*d_tri)^2), nrow = c*d_tri, ncol = c*d_tri)
  }
  
  nT2 = list()
  B_j = list()
  W_hat_j = list()
  # Sigma_hat = (1/n)*t(W)%*%W
  
  aux_Sigma_hat = list()
  
  for(k in 1:r){
    B_j[[k]] = 20*aux_B_j[[k]]/(n_esti^2)
    nT2[[k]] = 20*aux_inte[[k]]/(n_esti^2)
    nT2[[k]] = n*nT2[[k]]
    aux_W  = matrix(W[,((k-1)*c*d_tri + 1):(k*c*d_tri)], ncol = c*d_tri)
    aux_Sigma_hat[[k]] = 1/n*t(aux_W)%*%aux_W
  }
  
  for(k in 1:r){
    W_hat_j[[k]] = solve(aux_Sigma_hat[[k]],  tol = 1e-40)%*%B_j[[k]]%*%solve(aux_Sigma_hat[[k]],  tol = 1e-40)
  }
  
  
  #df_basis = 0
  #if(b_timese == "Cspli"){
  #  df_basis = Cspline_table(d) # true is d - 2 + or 
  #  d = dim(df_basis)[2]
  #}
  
  # basis_x = matrix(ncol = d)
  n_esti = n-r
  basis_x = matrix(nrow = n_esti, ncol = 1) 
  aux_phi_x = matrix(ncol = d_tri)
  
  for(j in 1:r){ # j for lag
    for(i in (r+1):n){
      aux_phi_x = rbind(aux_phi_x, sqrt(mp_selection(mp_type, ts[i-j], s)[2])*select_basis_timese(d, mp_selection(mp_type, ts[i-j], s)[1], b_timese, df_basis_d))
    }
    aux_phi_x = aux_phi_x[-1, ]
    aux_phi_x = as.matrix(aux_phi_x, ncol = d_tri)
    basis_x = cbind(basis_x, aux_phi_x)
    aux_phi_x = matrix(ncol = d_tri)
  }
  
  basis_x = matrix(basis_x[,-1], ncol = r*d_tri)
  
  #for(i in (r+1):n){
  #  basis_x = rbind(basis_x, sqrt(mp_selection(mp_type, ts[i-1], s)[2])*select_basis_timese(d, mp_selection(mp_type, ts[i-1], s)[1], b_timese, df_basis))
  #}
  
  #basis_x = basis_x[-1, ]
  #basis_x = as.matrix(basis_x, ncol = d)
  
  aux_bt  = matrix(nrow = n_esti)
  b = matrix(nrow = n_esti, ncol = 1)
  
  if(b_time == "Cspli"){
    df_basis = Cspline_table(aux_c, n) # true is d - 2 + or 
    c = dim(df_basis)[2]
  }
  
  if(b_time %in% wavelet_basis){
    df_basis_1 = wavelet_kth_b(aux_c, n, b_time)
    c = dim(df_basis_1)[2]
  }
  
  for(k in 1:r){
    for(l1 in 1:c){
      if(b_time == "Cspli"){
        basis_t = as.matrix(df_basis[(r+1):n, l1])
      }else if(b_time %in% wavelet_basis){
        basis_t = as.matrix(df_basis_1[(r+1):n, l1])
      }else{
        basis_t = bs.gene(b_time, l1, n)[(r+1):n, 1]
      }
      aux_bt = cbind(aux_bt, basis_t)
      for(l2 in 1:d_tri){
        b = cbind(b, matrix(basis_t*basis_x[, (k-1)*d_tri + l2], ncol = 1))
      }
    }
  }
  
  
  colnames(b) = NULL
  b = as.matrix(b[,-1], ncol = r*c*d_tri) 
  
  colnames(aux_bt) = NULL
  aux_bt = aux_bt[,-1]
  aux_bt = as.matrix(aux_bt[ , 1:c], ncol = c)
  
  # I_r = diag(r)
  # l_j = diag(c*d)  
  # aux.esti = t(t(l_j%*%beta_hat)%*%(t(kronecker(b, I_r))))[,1]
  esti = list()
  for(k in 1:r){
    aux_beta_hat = matrix(beta_hat[((k-1)*c*d_tri+1):(k*c*d_tri), 1], ncol = 1)
    # esti[[k]] = t(t(l_j%*%aux_beta_hat)%*%(t(kronecker(b, I_r))))[,1]
    
    aux_b = matrix(b[ ,((k-1)*c*d_tri + 1):(k*c*d_tri)], ncol = c*d_tri)
    esti[[k]] = (aux_b%*%aux_beta_hat)[,1]
  }
  
  aux.esti = rep(0, n_esti)
  for(i in 1:n_esti){
    for(k in 1:r){
      aux.esti[i] = aux.esti[i] + esti[[k]][i] 
    }
  }
  
  T_2k = list()
  
  
  
  # select m
  if(m == "MV"){
    m = min_vola(ts, aux.esti, basis_x, aux_bt, r)
  }
  
  # bootstrap 
  aux_pvalue = list()
  
  for(index_r in 1:r){
    for(k in 1:1000){
      Xi_i = matrix(rep(0, c*d_tri), ncol = 1)
      for(i in (r+1):(n-m)){
        R_i = rnorm(1, 0, 1)
        z_j = (ts[i] - aux.esti[i-r])*matrix(basis_x[i-r, ((index_r-1)*d_tri + 1):(index_r*d_tri)],ncol = 1)
        for(j in (i+1):(i+m)){
          z_j = z_j + (ts[j] - aux.esti[j-r])*matrix(basis_x[j-r, ((index_r-1)*d_tri + 1):(index_r*d_tri)], ncol = 1)
        }
        a_ti = matrix(as.numeric(aux_bt[i-r, ]), ncol =1)
        Xi_i = Xi_i + R_i*kronecker(z_j, a_ti)
      }
      Xi_i = (1/(sqrt((n-m-r)*(m))))*Xi_i
      T_2k[[k]] = t(Xi_i)%*%W_hat_j[[index_r]]%*%Xi_i
    }
    aux_pvalue[[index_r]] = 1- ifelse(length(which(sort(unlist(T_2k)) <= as.numeric(nT2[[index_r]]))) == 0, 0, max(which(sort(unlist(T_2k)) <= as.numeric(nT2[[index_r]]))))/1000
  }
  return(aux_pvalue)
} # you need to perform your own function for input. 


#  minimum volatility (MV) method
min_vola <- function(ts, aux.esti, basis_x, aux_bt, r){
  h.0 = 3
  m.li = c(1:25)
  #library(Matrix)
  # Design matrix
  # Y =  beta_l(timese, c, b)[[2]]
  n = length(ts)
  # li.res = list()
  # m = 6
  
  # Error, i = b* + 1... n
  #error.s = c()
  #es.alpha = # alpha.legen(timese, c, b, n)
  #aux.len = length(es.alpha)
  #for(i in (b+1):n){
  #  val.aux = es.alpha[[1]][i]
  #  for(j in 2:aux.len){
  #    val.aux = val.aux + es.alpha[[j]][i]*timese[i-j+1]
  #  }
  #  error.s[i-b] = timese[i] - val.aux
  # }
  

  Phi.li = list()
  for (m in m.li){
    aux_Phi = 0
    Phi = 0
    for(i in (r+1):(n-m)){
      z_j = (ts[i]-aux.esti[i-r])*matrix(basis_x[i-r, ],ncol = 1)
      for(j in i:(i+m)){
        z_j = z_j + (ts[j]-aux.esti[j-r])*matrix(basis_x[j-r, ], ncol = 1)
      }
      a_ti = matrix(aux_bt[i-r,], ncol =1)
      Phi = Phi + kronecker(z_j, a_ti)
      aux_Phi = aux_Phi + Phi%*%t(Phi)
    }
    Phi.li[[m]] = 1/((n-m-r+1)*m)*aux_Phi
  }
  
  se.li = list()
  for(mj in (min(m.li)+h.0):(max(m.li)-h.0)){
    av.Phi = 0
    se = 0
    for (k in -3:3){
      av.Phi = av.Phi + Phi.li[[mj + k]]
    }
    av.Phi = av.Phi/7
    
    for(k in -3:3){
      se = se + norm(av.Phi - Phi.li[[mj + k]], "2")^2
    }
    se.li[[mj-3]] = sqrt(se/6)
  }
  return(m.op = which(unlist(se.li) == min(unlist(se.li))) + 3)
}



