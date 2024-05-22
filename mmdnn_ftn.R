

obs_Overlap = function(u, v, Masking){
  # Takes two rows u, v with masking matrix and returns column indices 
  # where entries for rows u, v are both observed
  overlap = Masking[u, ] == Masking[v, ]
  return(overlap)
}

hatMMD_square = function(dat1, dat2){
  # Takes two data set and returns U-statistics estimator for MMD squared distance 
  val = 0
  N = nrow(dat1)
  square_kernel = function(x, y) (sum(x*y) + 1)^2
  for(i in seq(1:(N - 1))){
    for(j in seq(from = i + 1, to = N)){
      val = val + square_kernel(dat1[i, ], dat1[j, ]) + square_kernel(dat2[i, ], dat2[j, ]) - 2*square_kernel(dat1[i, ], dat2[j, ])    
    }
  }
  return(2*val/(N*(N - 1)))
}

dissim_square = function(k, l, dim, Data, Masking, discard, u){
  # Takes rows k, l containing data of dimension dim
  # and returns dissimilarity measure between the two rows
  # using U-statistics MMD estimate - for Kernel-NN
  m = nrow(Masking)
  n = ncol(Masking)
  
  if(discard == TRUE){
    overlap = which(obs_Overlap(k, l, Masking))
    overlap = overlap[! overlap %in% c(u)] # discard reference column u
    if(length(overlap) == 0){
      val = 10^5
    }else if(length(overlap) > 0){
      urow_ind = (k - 1)*n + overlap
      vrow_ind = (l - 1)*n + overlap
      
      val = 0
      for(i in c(1 : length(overlap))){
        mat1 = matrix(Data[urow_ind[i], ], N, dim, byrow = TRUE)
        mat2 = matrix(Data[vrow_ind[i], ], N, dim, byrow = TRUE)
        val = val + hatMMD_square(mat1, mat2)
      } 
      val = val/length(overlap)
    }
  }
  else if(discard == FALSE){
    overlap = which(obs_Overlap(k, l, Mask_train))
    if(length(overlap) == 0){
      val = 10^5
    }else if(length(overlap) > 0){
      urow_ind = (k - 1)*n + overlap
      vrow_ind = (l - 1)*n + overlap
      
      val = 0
      for(i in c(1 : length(overlap))){
        mat1 = matrix(Data[urow_ind[i], ], N, dim, byrow = TRUE)
        mat2 = matrix(Data[vrow_ind[i], ], N, dim, byrow = TRUE)
        val = val + hatMMD_square(mat1, mat2)
      } 
      val = val/length(overlap)
    }
  }
  return(val) 
  
}

ave_dim = function(vec, N, dim) apply(matrix(vec, N, dim, byrow = TRUE), 2, mean)

dissim_mean = function(k, l, dim, Data, Masking, discard, u){
  # Takes rows k, l containing data of dimension dim 
  # and returns dissimilarity measure between rows k, l
  # using simple squared distance - for Vanilla NN
  N = ncol(Data)/dim
  m = nrow(Masking)
  n = ncol(Masking)
  
  if(dim == 1){
    Data_ave = matrix(apply(Data, 1, function(x) ave_dim(x, N, dim)), (m * n), 1)
  }else if(dim > 1){
    Data_ave = t(apply(Data, 1, function(x) ave_dim(x, N, dim))) 
  }
  
  if(discard == TRUE){
    overlap = which(obs_Overlap(k, l, Masking))
    overlap = overlap[! overlap %in% c(u)] # discard reference column u
    if(length(overlap) == 0){
      val = 10^5
    }else if(length(overlap) > 0){
      urow_ind = (k - 1)*n + overlap
      vrow_ind = (l - 1)*n + overlap
      
      mat1 = Data_ave[urow_ind, ]
      mat2 = Data_ave[vrow_ind, ]
      
      if(length(overlap) == 1){
        mat1 = matrix(mat1, 1, dim)
        mat2 = matrix(mat2, 1, dim)
      }
      if(dim == 1){
        mat1 = matrix(mat1, length(urow_ind), 1)
        mat2 = matrix(mat2, length(urow_ind), 1)
      }
      
      val = sum(apply(mat1 - mat2, 1, function(x) sum(x^2)))/length(overlap) 
    }
  }
  else if(discard == FALSE){
    overlap = which(obs_Overlap(k, l, Masking))
    if(length(overlap) == 0){
      val = 10^5
    }else if(length(overlap) > 0){
      urow_ind = (k - 1)*n + overlap
      vrow_ind = (l - 1)*n + overlap
      
      mat1 = Data_ave[urow_ind, ]
      mat2 = Data_ave[vrow_ind, ]
      
      if(length(overlap) == 1){
        mat1 = matrix(mat1, 1, dim)
        mat2 = matrix(mat2, 1, dim)
      }
      
      if(dim == 1){
        mat1 = matrix(mat1, length(urow_ind), 1)
        mat2 = matrix(mat2, length(urow_ind), 1)
      }
      
      val = sum(apply(mat1 - mat2, 1, function(x) sum(x^2)))/length(overlap) 
    }
  }
  return(val) 
}

dissim_mat = function(Data_train, Mask_train, m, dim){
  # Returns dissimilarity between all the rows for both Kernel-NN and Vanilla NN
  dist_mat = matrix(0, m, m)
  mean_mat = matrix(0, m, m)
  for(k in c(1:(m-1))){
    for(l in c((k+1):m)){
      dist_mat[k, l] = dissim_square(k, l, dim, Data_train, Mask_train, discard = FALSE, u = 1)    
      mean_mat[k, l] = dissim_mean(k, l, dim, Data_train, Mask_train, discard = FALSE, u = 1)
    }
  }
  dist_mat = t(dist_mat) + dist_mat
  mean_mat = t(mean_mat) + mean_mat
  
  out = list(dist_mat, mean_mat)
  return(out)
}

mixture_gen = function(Data, N, dim){
  # Returns samples of mixture of empirical distributions
  nbhd_size = nrow(Data)
  NN = ncol(Data)/dim
  mixture = matrix(0, dim, N)
  sel_comp = sample(c(1:NN), N, replace = TRUE)
  for(i in seq(1:N)){
    sel_mixing_i = sample((1:nbhd_size), 1)
    sel_comp_i = sel_comp[i]
    mixture[, i] = Data[sel_mixing_i, c(((sel_comp_i - 1)*dim + 1) : (sel_comp_i*dim))]
  }  
  return(mixture)
}

spread = function(dat){
  if(nrow(dat) > 1){
    spreaded = c()
    for(i in c(1:nrow(dat))){
      spreaded = c(spreaded, dat[i, ])
    }
  }else if(nrow(dat) == 1){
    spreaded = dat
  }
  return(spreaded)
}

# Kernels
square_kernel = function(x, y) (sum(x*y) + 1)^2
exp_kernel = function(x, y) exp(-(0.5)*(sum((x - y)*(x - y))))


ustat_exp = function(dat1, dat2, dim){
  # Returns U-statistics of two d-dimensional datasets of exponential kernel
  N = ncol(dat1)/dim
  nbhd_size = nrow(dat2)
  Ne = N*nbhd_size
  flat_dat2 = spread(dat2)
  kern = exp_kernel
  
  mat1 = matrix(0, N, N)
  mat2 = matrix(0, Ne, Ne)
  for(i in c(1:N)){
    for(j in c(1:N)){
      mat1[i, j] = kern(dat1[c((((i-1)*dim + 1) : (i*dim)))], dat1[c((((j-1)*dim + 1) : (j*dim)))])
    }
  }
  for(i in c(1:Ne)){
    for(j in c(1:Ne)){
      mat2[i, j] =  kern(flat_dat2[c((((i-1)*dim + 1) : (i*dim)))], flat_dat2[c((((j-1)*dim + 1) : (j*dim)))])
    }
  }
  crosssum = 0
  for(i in c(1:N)){
    for(j in c(1:Ne)){
      crosssum = crosssum + kern(dat1[c((((i-1)*dim + 1) : (i*dim)))], flat_dat2[c((((j-1)*dim + 1) : (j*dim)))])
    }
  }
  
  val = (sum(mat1) - sum(diag(mat1)))/(N*(N - 1)) + (sum(mat2) - sum(diag(mat2)))/(Ne*(Ne - 1)) - 2*crosssum/(N*Ne)
  return(val)
}

ustat_square = function(dat1, dat2, dim){
  # Returns U-statistics of two d-dimensional datasets of squared poly kernel
  N = ncol(dat1)/dim
  nbhd_size = nrow(dat2)
  Ne = N*nbhd_size
  flat_dat2 = spread(dat2)
  kern = square_kernel
  
  mat1 = matrix(0, N, N)
  mat2 = matrix(0, Ne, Ne)
  for(i in c(1:N)){
    for(j in c(1:N)){
      mat1[i, j] = kern(dat1[c((((i-1)*dim + 1) : (i*dim)))], dat1[c((((j-1)*dim + 1) : (j*dim)))])
    }
  }
  for(i in c(1:Ne)){
    for(j in c(1:Ne)){
      mat2[i, j] =  kern(flat_dat2[c((((i-1)*dim + 1) : (i*dim)))], flat_dat2[c((((j-1)*dim + 1) : (j*dim)))])
    }
  }
  crosssum = 0
  for(i in c(1:N)){
    for(j in c(1:Ne)){
      crosssum = crosssum + kern(dat1[c((((i-1)*dim + 1) : (i*dim)))], flat_dat2[c((((j-1)*dim + 1) : (j*dim)))])
    }
  }
  
  val = (sum(mat1) - sum(diag(mat1)))/(N*(N - 1)) + (sum(mat2) - sum(diag(mat2)))/(Ne*(Ne - 1)) - 2*crosssum/(N*Ne)
  return(val)
}

mmdnn = function(Data, n, nbhd, k, l){
  # Returns the neighbors for Kernel-NN
  if(length(nbhd) == 0){
    mix = matrix(0, 1, ncol(Data))
  }else if(length(nbhd) == 1){
    mix = matrix(Data[(n*(nbhd - 1) + l), ], 1, ncol(Data))
  }else if(length(nbhd) > 1){
    mix = Data[(n*(nbhd - 1) + l), ]
  }
  return(mix)
}

vanillann = function(Data, n, dim, nbhd, k, l){
  # Returns the neighbors for Vanilla NN
  N = ncol(Data)/dim
  if(length(nbhd) == 0){
    res = rep(0, dim)
  }else if(length(nbhd) == 1){
    res = ave_dim(Data[(n*(nbhd - 1) + l), ], N, dim)
  }else if(length(nbhd) > 1){
    res = t(apply(Data[(n*(nbhd - 1) + l), ], 1, function(x) ave_dim(x, N, dim)))
    res = apply(res, 2, mean)
  }
  return(res)
}


NNcv = function(dissim_dist, dissim_mean, Data_test, Mask_test, dim, eta){
  # Takes in dissimilarity matrix constructed from training data
  # returns the performance of Kernel-NN for test data for each radius eta
  N = ncol(Data_test)/dim
  m = nrow(Mask_test)
  n_t = ncol(Mask_test)
  
  ave_size = 0
  sigmaE_d = sigmaS_d = sigma_m = 0
  
  for(k in c(1:m)){
    for(l in c(1:n_t)){
      if(Mask_test[k, l] == 1){
        nbhd_dist_kl = intersect(setdiff(which(dissim_dist[, k] < eta), c(k)), setdiff(which(Mask_test[, l] == 1), c(k)))
        nbhd_mean_kl = intersect(setdiff(which(dissim_mean[, k] < eta), c(k)), setdiff(which(Mask_test[, l] == 1), c(k)))
        
        Data_dist_kl = matrix(Data_test[(n_t*(k - 1) + l), ], 1, N*dim) # Data corresponding to (k, l)
        Data_mean_kl = c(apply(Data_dist_kl, 1, function(x) ave_dim(x, N, dim))) 
        
        MMD_mix = mmdnn(Data_test, n_t, nbhd_dist_kl, k, l)
        VAN_mean = vanillann(Data_test, n_t, dim, nbhd_mean_kl, k, l)
        # Neighbor_dist_kl = Data_test[n_t*(nbhd_dist_kl - 1) + l, ]
        # Neighbor_mean_kl = Data_test[n_t*(nbhd_mean_kl - 1) + l, ]
        # Neighbor_mean_kl = apply(Neighbor_mean_kl, 1, function(x) ave_dim(x, N, dim))
        
        sigmaE_d_kl = ustat_exp(Data_dist_kl, MMD_mix, dim)
        sigmaS_d_kl = ustat_square(Data_dist_kl, MMD_mix, dim)
        sigma_m_kl = sum((VAN_mean - Data_mean_kl)^2)
        
        sigmaE_d = sigmaE_d + sigmaE_d_kl
        sigmaS_d = sigmaS_d + sigmaS_d_kl
        sigma_m = sigma_m + sigma_m_kl
        
        ave_size = ave_size + 1
      }
    }
  }
  sigmaE_d = sigmaE_d/ave_size
  sigmaS_d = sigmaS_d/ave_size
  sigma_m = sigma_m/ave_size
  
  return(c(sigmaE_d, sigmaS_d, sigma_m))
}


# Below are all for calculating MMD squared distance between mixture distribution
# and Gaussian distribution with diagonal covariance matrix

inn_ET = function(final_Matrix, true_Mean, true_Var){
  # final_Matrix : nbhd_size / (dim*N) matrix, true_Mean : d dimensional vector
  # true_Var : d dimensional diagonal variance(specialized only for diag Gaussian case!)
  nbhd_size = nrow(final_Matrix)
  dim = length(true_Mean)
  N = ncol(final_Matrix)/dim
  dist = 0
  true_Mean_mat = matrix(rep(true_Mean, N), dim, N, byrow = FALSE)
  
  for(i in c(1:nbhd_size)){
    samples_mat = matrix(final_Matrix[i, ], dim, N, byrow = FALSE)
    if(dim == 1){
      dist = dist + mean(exp((-0.5)*((true_Mean_mat - samples_mat)^2)/(1 + true_Var)))
    }else if(dim > 1){
      scaled_diff = t(true_Mean_mat - samples_mat)%*%diag(c(1/(1 + true_Var)))%*%(true_Mean_mat - samples_mat)
      dist = dist + sum(exp((-0.5) * diag(scaled_diff)))/N 
    }
  }
  
  dist = (dist/nrow(final_Matrix))/sqrt(prod(1 + true_Var))
  
  return(dist)
}

exp_kern = function(x, y) exp(-0.5*(sum((x - y)*(x - y))))

inn_EE = function(final_Matrix){
  
  nbhd_size = nrow(final_Matrix)
  dim = length(true_Mean)
  N = ncol(final_Matrix)/dim
  dist = 0
  
  for(i in c(1:nbhd_size)){
    for(j in c(1:nbhd_size)){
      samples_mat1 = matrix(final_Matrix[i, ], dim, N, byrow = FALSE)
      samples_mat2 = matrix(final_Matrix[j, ], dim, N, byrow = FALSE)
      for(k in c(1:N)){
        for(l in c(1:N)){
          dist = dist + exp_kern(samples_mat1[, k], samples_mat2[, l])
        }
      }
    }
  }
  
  dist = dist/((N^2) * (nbhd_size^2))
  
  return(dist)
}

Gkernel_metric = function(final_Matrix, true_Mean, true_Var){
  # MMD squared distance where onedistribution is empirical mixture measure
  # the other is Gaussian with diagonal Covariance
  inn_tt = 1/sqrt(prod(1 + 2*true_Var))
  inn_et = inn_ET(final_Matrix, true_Mean, true_Var)
  inn_ee = inn_EE(final_Matrix)
  
  dist = inn_tt + inn_ee - 2*inn_et
  return(dist)
}
