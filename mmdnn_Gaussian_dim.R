### Analyze Gaussian Location Scale Family

# Can analyze by changing rows(m) / dimension of data(d) / missingness type(gendata)

# Partition data by half w.r.t. column of matrix - so input matrix column should be EVEN!

rm(list = ls())
library(mvtnorm)
library(ggplot2)
library(dplyr)
library(xtable)
library(MASS)

setwd("/Users/kyuseongchoi/Desktop/DNN")

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

source('MMD-DNN-DGP.R') # Source file `MMD-DNN-DGP.R` generates data
source('mmdnn_ftn.R') # Source file `mmdnn_ftn.R` contains necessary functions to run kernel-NN

m = 64 # Number of rows
n = 20 # Number of columns
N = 30 # Number of data points in each matrix entries
nsim = 25 # Number of simulations ran
eta_pool = c(exp(c(-2.8, -2.6, -2.4, -2.2)), exp(seq(-2, 0, by = 0.3))) # Candidates for optimal radius eta


for(dim in 2^c(0, 3)){
  sigmaS_traj = sigmaE_traj = sigmam_traj = matrix(0, nsim, length(eta_pool))
  etaE_star_pool = etaS_star_pool = etam_star_pool = c()
  mmdperf_E = mmdperf_S = vanperf = dist_vanperf = c()
  nbhd_Size = matrix(0, nsim, 3)
  
  par(mfrow = c(5, 5))
  par(mar=c(1,1,1,1))
  
  true_Mean_mat = matrix(0, nsim, dim)
  true_Cov_mat = matrix(0, nsim, dim)
  
  hatA_E_collective = list()
  hatA_S_collective = list()
  
  for(sim in c(1:nsim)){
    init = gendata3(m, n, N, dim, seed = sim)
    Data = init$Data
    Masking = init$Masking
    true_Mean = init$target_mean
    true_Cov = init$target_cov
    true_Mean_mat[sim, ] = true_Mean
    true_Cov_mat[sim, ] = true_Cov
    
    ind1 = c()
    for(i in c(1:m)){
      ind1 = c(ind1, c(((i-1)*n + 1) : (((i-1)*n) + n/2)))
    }
    ind2 = setdiff(c(1:(m*n)), ind1)
    
    Data_train = Data[ind1, ]
    Mask_train = Masking[, c(1 : (n/2))]
    Data_test = Data[ind2, ]
    Mask_test = Masking[, c(((n/2) + 1) : n)]
    
    dissim_out = dissim_mat(Data_train, Mask_train, m, dim)
    dist_dist = dissim_out[[1]]
    dist_mean = dissim_out[[2]]
    
    sigmaS_d_pool = sigmaE_d_pool = sigma_m_pool = c()
    for(eta in eta_pool){
      a = Sys.time()
      
      sigma = NNcv(dist_dist, dist_mean, Data_test, Mask_test, dim, eta)
      sigmaE_d_pool = c(sigmaE_d_pool, sigma[1])
      sigmaS_d_pool = c(sigmaS_d_pool, sigma[2])
      sigma_m_pool = c(sigma_m_pool, sigma[3])
      
      b = Sys.time()
      cat((eta), "iteration:", b-a, "seconds", "\n")
    }
    
    sigmaS_traj[sim, ] = sigmaS_d_pool
    sigmaE_traj[sim, ] = sigmaE_d_pool
    sigmam_traj[sim, ] = sigma_m_pool
    
    etaE_star = eta_pool[which.min(sigmaE_d_pool)]
    etaS_star = eta_pool[which.min(sigmaS_d_pool)]
    etam_star = eta_pool[which.min(sigma_m_pool)]
    
    etaE_star_pool = c(etaE_star_pool, etaE_star)
    etaS_star_pool = c(etaS_star_pool, etaS_star)
    etam_star_pool = c(etam_star_pool, etam_star)
    
    rho = c()
    hor = c()
    cand_rows = c(2:m)
    for(i in c(2:m)){
      rho_1i = dissim_square(1, i, dim, Data, Masking, discard = TRUE, u = 1)
      hor_1i = dissim_mean(1, i, dim, Data, Masking, discard = TRUE, u = 1)
      rho = c(rho, rho_1i)
      hor = c(hor, hor_1i)
    }
    
    setdiff(which(Masking[, 1] == 1), c(1))
    
    final_nbhdE = intersect(cand_rows[which(rho < etaE_star)], setdiff(which(Masking[, 1] == 1), c(1))) 
    final_nbhdS = intersect(cand_rows[which(rho < etaS_star)], setdiff(which(Masking[, 1] == 1), c(1)))
    final_nbhd_m =intersect(cand_rows[which(hor < etam_star)], setdiff(which(Masking[, 1] == 1), c(1)))
    
    nbhd_Size[sim, ] = c(length(final_nbhdE), length(final_nbhdS), length(final_nbhd_m))
    
    hatA_E = mmdnn(Data, n, final_nbhdE, 1, 1)
    hatA_S = mmdnn(Data, n, final_nbhdS, 1, 1)
    hatA_m = vanillann(Data, n, dim, final_nbhd_m, 1, 1)
    
    hatA_E_collective[[sim]] = hatA_E
    hatA_S_collective[[sim]] = hatA_S
    
    # Performance of hatA(1, 1) for learning Gaussian loc/scale --- let's change dim!!
    perf_dist_E = Gkernel_metric(hatA_E, true_Mean, true_Cov)
    perf_dist_S = Gkernel_metric(hatA_S, true_Mean, true_Cov)
    
    mmdperf_E = c(mmdperf_E, perf_dist_E)
    mmdperf_S = c(mmdperf_S, perf_dist_S)
    
    #Compare the mean learning performance
    perf_mean = sum((hatA_m - true_Mean)^2)
    perf_meanviadist = sum((vanillann(Data, n, dim, final_nbhdE, 1, 1) - true_Mean)^2)
    
    vanperf = c(vanperf, perf_mean)
    dist_vanperf = c(dist_vanperf, perf_meanviadist)
    
    # Denoising aspect
    if(dim == 1){
     arg = seq(-3, 3, by = 0.1) 
     prob = pnorm(arg, mean = true_Mean, sd = sqrt(true_Cov))
     plot(arg, prob, type = "l", col = "black", lwd = 2)
     empirical_dist = matrix(Data[1, ], N/2, dim)
     hatA_E_samples = t(mixture_gen(hatA_E, N/2, dim))
     lines(ecdf(empirical_dist), col = "red")
     lines(ecdf(hatA_E_samples), col = "blue")
    }
    
    if(dim == 2){
      hatA_E_samples = data.frame(t(mixture_gen(hatA_E, N/2, dim)))
      empirical_dist = data.frame(t(matrix(Data[1, ], dim, N/2, byrow = FALSE)))
      true_density = function(x_1, x_2){
        m1 = true_Mean[1]
        m2 = true_Mean[2]
        v1 = true_Cov[1]
        v2 = true_Cov[2]
        ((2*pi)^{-1})*((prod(true_Cov))^{-1/2})*exp(-(0.5)*(((x_1 - m1)^2)/v1 + ((x_2 - m2)^2)/v2))
      }
      x_1 = x_2 = seq(-4, 4, by = 0.2)
      
      zz = outer(x_1, x_2, true_density)
      contour(x_1, x_2, zz)
      points(empirical_dist[, 1], empirical_dist[, 2], col = "red")
      points(hatA_E_samples[, 1], hatA_E_samples[, 2], col = "blue") 
    }
    # Compare the following three distributions
    # mixture_gen(Data[(n*(final_nbhdE - 1) + l), ], N, dim) # Samples from hatA(1, 1)
    # Data[1, ] # Observed empirical dist
    # rmvnorm(N, true_Mean, diag(true_Cov))
  }
  
  sigma_traj = cbind(sigmaS_traj, sigmaE_traj, sigmam_traj) # Save this
  eta_star_pool = rbind(etaE_star_pool, etaS_star_pool, etam_star_pool) # Save this
  mmdperf = cbind(mmdperf_E, mmdperf_S)
  meanperf = cbind(vanperf, dist_vanperf)
  true_info = cbind(true_Mean_mat, true_Cov_mat)
  
  Collective_Info = list(sigma_traj, eta_star_pool, nbhd_Size, mmdperf, meanperf, true_info)
  Denoising_Info = list(hatA_E_collective, hatA_S_collective, true_info)
  file_Name = paste("m =", m, ",n =", n, ",N =", N, ",dim =", dim, "nsim =", nsim, ".Rdata")
  file_Name2 = paste("m =", m, ",n =", n, ",N =", N, ",dim =", dim, "nsim =", nsim, ",denoising.Rdata")
  save(Collective_Info, file = file_Name)
  save(Denoising_Info, file = file_Name2)
}

## Note : subscript S, E, m for each item corresponds to DNN via squared poly kernel / DNN via exponential kernel / Vanilla NN
## Output :
# sigma_traj - tracks CV statistics sigma for kernel-NN and vanilla-NN
# eta_star_pool - records optimal radius chosen
# mmdperf - MMD performance of kernel-NN when using exponential kernel or squared poly kernel
# meanperf - MSE performance of vanilla-NN
# true_info - true mean and variance of Gaussian distributions for entry (1, 1)




