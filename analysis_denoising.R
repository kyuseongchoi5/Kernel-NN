# load Data that is automatically saved from `mmdnn_Gaussian_rows` under name 'Denoising_Info'
# analyze the denoising property of Kernel-NN

rm(list = ls())


library(ggplot2)
library(dplyr)  
library(mvtnorm)

setwd('/Users/choiqs/Desktop/DNN') # Set working dir
source('/Users/choiqs/Desktop/DNN/mmdnn_ftn.R') # Source file `mmdnn_ftn.R`

load("m = 64 ,n = 20 ,N = 30 ,dim = 1 nsim = 25 ,denoising.Rdata")
hatA_E_d1 = Denoising_Info[[1]]
hatA_S_d1 = Denoising_Info[[2]]
true_d1 = Denoising_Info[[3]]

load("m = 64 ,n = 20 ,N = 30 ,dim = 2 nsim = 25 ,denoising.Rdata")
hatA_E_d2 = Denoising_Info[[1]]
hatA_S_d2 = Denoising_Info[[2]]
true_d2 = Denoising_Info[[3]]

load("m = 64 ,n = 20 ,N = 30 ,dim = 4 nsim = 25 ,denoising.Rdata")
hatA_E_d4 = Denoising_Info[[1]]
hatA_S_d4 = Denoising_Info[[2]]
true_d4 = Denoising_Info[[3]]

load("m = 64 ,n = 20 ,N = 30 ,dim = 8 nsim = 25 ,denoising.Rdata")
hatA_E_d8 = Denoising_Info[[1]]
hatA_S_d8 = Denoising_Info[[2]]
true_d8 = Denoising_Info[[3]]


## 1d-denoising

N = 30
nsim = 25


ks_hatA_d1 = c()
ks_emp_d1 = c()

#ks_hatA_mat_d1
#ks_emp_mat_d1

for(sim in c(1:nsim)){
  hatA_E = hatA_E_d1[[sim]]
  hatA_S = hatA_S_d1[[sim]]
  
  nbhd_size = nrow(hatA_E)
  
  true_mean = true_d1[sim, 1]; true_var = true_d1[sim, 2]
  
  hatA_E_long = matrix(0, 1, N*nbhd_size)
  
  arg = seq(-2, 2, by = 0.1) 
  prob = pnorm(arg, mean = true_mean, sd = sqrt(true_var))
  plot(arg, prob, type = "l", col = "black", lwd = 2)
  true_samples = rnorm(10^3, true_mean, sqrt(true_var))
  empirical_dist = rnorm(N, true_mean, sqrt(true_var))
  
  hatA_E_long = matrix(0, 1, N*nbhd_size)
  for(i in c(1:nbhd_size)){
    for(j in c(1:N)){
      hatA_E_long[, ((i - 1)*N + j)] = hatA_E[i, j] 
    }
  } 
  
  ks_emp = ks.test(empirical_dist, true_samples)$statistic
  ks_hatA = ks.test(hatA_E_long, true_samples)$statistic
  
  ks_emp_d1 = c(ks_emp_d1, ks_emp)
  ks_hatA_d1 = c(ks_hatA_d1, ks_hatA)
  
  lines(ecdf(empirical_dist), col = "red")
  lines(ecdf(c(hatA_E)), col = "blue")
}


hatA_ci_d1 = c(mean(ks_hatA_d1), mean(ks_hatA_d1) - 1.96*sd(ks_hatA_d1), mean(ks_hatA_d1) + 1.96*sd(ks_hatA_d1))
emp_ci_d1 = c(mean(ks_emp_d1), mean(ks_emp_d1) - 1.96*sd(ks_emp_d1), mean(ks_emp_d1) + 1.96*sd(ks_emp_d1))


plot(density(ks_hatA_d1))
lines(density(ks_emp_d1), col = "red")


## 2D-denoising

ks_hatA_mat_d2 = matrix(0, nsim, 2)
ks_emp_mat_d2 = matrix(0, nsim, 2)

for(sim in c(1:nsim)){
  hatA_E = hatA_E_d2[[sim]]
  hatA_S = hatA_S_d2[[sim]]
  true_mean = true_d2[sim, c(1, 2)]; true_var = true_d2[sim, c(3, 4)]
  
  # hatA_E_samples = t(mixture_gen(hatA_E, N, dim = 2))
  nbhd_size = nrow(hatA_E)
  hatA_E_long = matrix(0, 2, N*nbhd_size)
  for(i in c(1:nbhd_size)){
    for(j in c(1:N)){
      hatA_E_long[, ((i - 1)*N + j)] = hatA_E[i, c((2*j - 1), (2*j))] 
    }
  }  
  hatA_E_long = t(hatA_E_long)
  
  empirical_dist = rmvnorm(N, mean = true_mean, sigma = diag(true_var))
  true_density = function(x_1, x_2){
    m1 = true_mean[1]
    m2 = true_mean[2]
    v1 = true_var[1]
    v2 = true_var[2]
    ((2*pi)^{-1})*((prod(true_var))^{-1/2})*exp(-(0.5)*(((x_1 - m1)^2)/v1 + ((x_2 - m2)^2)/v2))
  }
  
  ks_hatA_mat_d2[sim, 1] = sum((apply(hatA_E_long, 2, mean) - true_mean)^2)
  ks_hatA_mat_d2[sim, 2] = norm(cov(hatA_E_long) - diag(true_var), "F")
  
  ks_emp_mat_d2[sim, 1] = sum((apply(empirical_dist, 2, mean) - true_mean)^2)
  ks_emp_mat_d2[sim, 2] = norm(cov(empirical_dist) - diag(true_var), "F")
  
  x_1 = x_2 = seq(-2, 2, by = 0.2)
  
  zz = outer(x_1, x_2, true_density)
  contour(x_1, x_2, zz)
  points(empirical_dist[, 1], empirical_dist[, 2], col = "red")
  points(hatA_E_long[, 1], hatA_E_long[, 2], col = "blue") 
}

print(ks_hatA_mat_d2)
print(ks_emp_mat_d2)

perf_hatA_d2 = (ks_hatA_mat_d2[, 1])^2 + ks_hatA_mat_d2[, 2]
perf_emp_d2 = (ks_emp_mat_d2[, 1])^2 + ks_emp_mat_d2[, 2]

print(c(mean(perf_hatA_d2), mean(perf_emp_d2)))

hatA_ci_d2 = c(mean(perf_hatA_d2), mean(perf_hatA_d2) - 1.96*sd(perf_hatA_d2), mean(perf_hatA_d2) + 1.96*sd(perf_hatA_d2))
emp_ci_d2 = c(mean(perf_emp_d2), mean(perf_emp_d2) - 1.96*sd(perf_emp_d2), mean(perf_emp_d2) + 1.96*sd(perf_emp_d2))

print(rbind(hatA_ci_d2, emp_ci_d2))

plot(density(perf_hatA_d2))
lines(density(perf_emp_d2), col = "red")


## 4d-denoising

nsim = 25
N = 30

ks_hatA_mat_d4 = matrix(0, nsim, 2)
ks_emp_mat_d4 = matrix(0, nsim, 2)

for(sim in c(1:nsim)){
  hatA_E = hatA_E_d4[[sim]]
  hatA_S = hatA_S_d4[[sim]]
  true_mean = true_d4[sim, c(1 : 4)]; true_var = true_d4[sim, c(5 : 8)]
  
  # hatA_E_samples = t(mixture_gen(hatA_E, N, dim = 2))
  nbhd_size = nrow(hatA_E)
  hatA_E_long = matrix(0, 4, N*nbhd_size)
  for(i in c(1:nbhd_size)){
    for(j in c(1:N)){
      hatA_E_long[, ((i - 1)*N + j)] = hatA_E[i, c((4*j - 3) : (4*j))] 
    }
  }  
  hatA_E_long = t(hatA_E_long)
  
  empirical_dist = rmvnorm(N, mean = true_mean, sigma = diag(true_var))
  
  ks_hatA_mat_d4[sim, 1] = sum((apply(hatA_E_long, 2, mean) - true_mean)^2)
  ks_hatA_mat_d4[sim, 2] = norm(cov(hatA_E_long) - diag(true_var), "F")
  
  ks_emp_mat_d4[sim, 1] = sum((apply(empirical_dist, 2, mean) - true_mean)^2)
  ks_emp_mat_d4[sim, 2] = norm(cov(empirical_dist) - diag(true_var), "F")
}

print(ks_hatA_mat_d4)
print(ks_emp_mat_d4)

perf_hatA_d4 = (ks_hatA_mat_d4[, 1])^2 + ks_hatA_mat_d4[, 2]
perf_emp_d4 = (ks_emp_mat_d4[, 1])^2 + ks_emp_mat_d4[, 2]

print(c(mean(perf_hatA_d4), mean(perf_emp_d4)))

hatA_ci_d4 = c(mean(perf_hatA_d4), mean(perf_hatA_d4) - 1.96*sd(perf_hatA_d4), mean(perf_hatA_d4) + 1.96*sd(perf_hatA_d4))
emp_ci_d4 = c(mean(perf_emp_d4), mean(perf_emp_d4) - 1.96*sd(perf_emp_d4), mean(perf_emp_d4) + 1.96*sd(perf_emp_d4))

print(rbind(hatA_ci_d4, emp_ci_d4))

plot(density(perf_hatA_d4))
lines(density(perf_emp_d4), col = "red")



ks_hatA_mat_d8 = matrix(0, nsim, 2)
ks_emp_mat_d8 = matrix(0, nsim, 2)

for(sim in c(1:nsim)){
  hatA_E = hatA_E_d8[[sim]]
  hatA_S = hatA_S_d8[[sim]]
  true_mean = true_d8[sim, c(1 : 8)]; true_var = true_d8[sim, c(9 : 16)]
  
  # hatA_E_samples = t(mixture_gen(hatA_E, N, dim = 2))
  nbhd_size = nrow(hatA_E)
  hatA_E_long = matrix(0, 8, N*nbhd_size)
  for(i in c(1:nbhd_size)){
    for(j in c(1:N)){
      hatA_E_long[, ((i - 1)*N + j)] = hatA_E[i, c((8*j - 7) : (8*j))] 
    }
  }  
  hatA_E_long = t(hatA_E_long)
  
  empirical_dist = rmvnorm(N, mean = true_mean, sigma = diag(true_var))
  
  ks_hatA_mat_d8[sim, 1] = sum((apply(hatA_E_long, 2, mean) - true_mean)^2)
  ks_hatA_mat_d8[sim, 2] = norm(cov(hatA_E_long) - diag(true_var), "F")
  
  ks_emp_mat_d8[sim, 1] = sum((apply(empirical_dist, 2, mean) - true_mean)^2)
  ks_emp_mat_d8[sim, 2] = norm(cov(empirical_dist) - diag(true_var), "F")
}

print(ks_hatA_mat_d8)
print(ks_emp_mat_d8)

perf_hatA_d8 = (ks_hatA_mat_d8[, 1])^2 + ks_hatA_mat_d8[, 2]
perf_emp_d8 = (ks_emp_mat_d8[, 1])^2 + ks_emp_mat_d8[, 2]


print(c(mean(perf_hatA_d8), mean(perf_emp_d8)))
print(c(sd(perf_hatA_d8), sd(perf_emp_d8)))

hatA_ci_d8 = c(mean(perf_hatA_d8), mean(perf_hatA_d8) - 1.96*sd(perf_hatA_d8), mean(perf_hatA_d8) + 1.96*sd(perf_hatA_d8))
emp_ci_d8 = c(mean(perf_emp_d8), mean(perf_emp_d8) - 1.96*sd(perf_emp_d8), mean(perf_emp_d8) + 1.96*sd(perf_emp_d8))


plot(density(perf_hatA_d8), main = "Closeness to true distribution - dimension 8")
lines(density(perf_emp_d8), col = "red")


hatA_center_ci = rbind(hatA_ci_d1, hatA_ci_d2, hatA_ci_d4, hatA_ci_d8)
hatA_center_ci = data.frame(rep("NN", 4), c(0:3), hatA_center_ci)
emp_center_ci = rbind(emp_ci_d1, emp_ci_d2, emp_ci_d4, emp_ci_d8)
emp_center_ci = data.frame(rep("Emp", 4), c(0:3), emp_center_ci)
colnames(hatA_center_ci) = c("Group", "Log dim", "MSE", "lower 95", "upper 95")
colnames(emp_center_ci) = c("Group", "Log dim", "MSE", "lower 95", "upper 95")

center_ci = rbind(hatA_center_ci, emp_center_ci)

pd = position_dodge(0.1)


center_ci %>%
  ggplot(aes(`Log dim`, MSE, colour = Group, group = Group)) + 
  geom_errorbar(aes(ymin = `lower 95`, ymax = `upper 95`), width=.1, position = pd) +
  geom_line(position = pd) +
  geom_point(aes(shape = Group), position = pd, size = 2) + theme_light(base_size = 17)

