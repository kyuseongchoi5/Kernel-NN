# load Data that is automatically saved from `mmdnn_Gaussian_rows` under name 'Collective_Info'
# analyze the performance of Kernel-NN

rm(list = ls())

library(ggplot2)
library(dplyr)


## MCAR

load("/Users/kyuseongchoi/Desktop/DNN/addsim/m =  4 , n = 10, N = 70, dim = 2 .Rdata")

sigma_traj_m4 = Collective_Info[[1]]
eta_star_pool_m4 = Collective_Info[[2]]
nbhd_Size_m4 = Collective_Info[[3]]
mmdperf_m4 = Collective_Info[[4]]
meanperf_m4 = Collective_Info[[5]]
trueinfo_m4 = Collective_Info[[6]]

load("/Users/kyuseongchoi/Desktop/DNN/addsim/m = 8, n = 10, N = 70, dim = 2 .Rdata")

sigma_traj_m8 = Collective_Info[[1]]
eta_star_pool_m8 = Collective_Info[[2]]
nbhd_Size_m8 = Collective_Info[[3]]
mmdperf_m8 = Collective_Info[[4]]
meanperf_m8 = Collective_Info[[5]]
trueinfo_m8 = Collective_Info[[6]]

load("/Users/kyuseongchoi/Desktop/DNN/addsim/m =  16 , n = 10, N = 70, dim = 2 .Rdata")

sigma_traj_m16 = Collective_Info[[1]]
eta_star_pool_m16 = Collective_Info[[2]]
nbhd_Size_m16 = Collective_Info[[3]]
mmdperf_m16 = Collective_Info[[4]]
meanperf_m16 = Collective_Info[[5]]
trueinfo_m16 = Collective_Info[[6]]

log_mmdperf_m4 = log(mmdperf_m4)
log_mmdperf_m8 = log(mmdperf_m8)
log_mmdperf_m16 = log(mmdperf_m16)

center_E_MCAR = c(mean(log_mmdperf_m4[, 1]), mean(log_mmdperf_m8[, 1]), mean(log_mmdperf_m16[, 1]))
Upper_MCAR = c(mean(log_mmdperf_m4[, 1]) +  1.96*sd(log_mmdperf_m4[, 1]), mean(log_mmdperf_m8[, 1]) +  1.96*sd(log_mmdperf_m8[, 1]), mean(log_mmdperf_m16[, 1]) +  1.96*sd(log_mmdperf_m16[, 1]))
Lower_MCAR = c(mean(log_mmdperf_m4[, 1]) -  1.96*sd(log_mmdperf_m4[, 1]), mean(log_mmdperf_m8[, 1]) - 1.96*sd(log_mmdperf_m8[, 1]), mean(log_mmdperf_m16[, 1]) -  1.96*sd(log_mmdperf_m16[, 1]))


## MNAR 
par(mfrow = c(1, 1))

load("/Users/kyuseongchoi/Desktop/DNN/addsim/m = 4 ,n = 10 ,N = 70 ,dim = 2 ,e_T = 0.5 ,e_N = 0.5 .Rdata")

sigma_traj_m4 = Collective_Info[[1]]
eta_star_pool_m4 = Collective_Info[[2]]
nbhd_Size_m4 = Collective_Info[[3]]
mmdperf_m4 = Collective_Info[[4]]
meanperf_m4 = Collective_Info[[5]]
trueinfo_m4 = Collective_Info[[6]]

load("/Users/kyuseongchoi/Desktop/DNN/addsim/m = 8 ,n = 10 ,N = 70 ,dim = 2 ,e_T = 0.5 ,e_N = 0.5 .Rdata")

sigma_traj_m8 = Collective_Info[[1]]
eta_star_pool_m8 = Collective_Info[[2]]
nbhd_Size_m8 = Collective_Info[[3]]
mmdperf_m8 = Collective_Info[[4]]
meanperf_m8 = Collective_Info[[5]]
trueinfo_m8 = Collective_Info[[6]]

load("/Users/kyuseongchoi/Desktop/DNN/addsim/m = 16 ,n = 10 ,N = 70 ,dim = 2 ,e_T = 0.5 ,e_N = 0.5 .Rdata")

sigma_traj_m16 = Collective_Info[[1]]
eta_star_pool_m16 = Collective_Info[[2]]
nbhd_Size_m16 = Collective_Info[[3]]
mmdperf_m16 = Collective_Info[[4]]
meanperf_m16 = Collective_Info[[5]]
trueinfo_m16 = Collective_Info[[6]]

log_mmdperf_m4 = log(mmdperf_m4)
log_mmdperf_m8 = log(mmdperf_m8)
log_mmdperf_m16 = log(mmdperf_m16)

center_E_MNAR = c(mean(log_mmdperf_m4[, 1]), mean(log_mmdperf_m8[, 1]), mean(log_mmdperf_m16[, 1]))
Upper_MNAR = c(mean(log_mmdperf_m4[, 1]) +  1.96*sd(log_mmdperf_m4[, 1]), mean(log_mmdperf_m8[, 1]) +  1.96*sd(log_mmdperf_m8[, 1]), mean(log_mmdperf_m16[, 1]) +  1.96*sd(log_mmdperf_m16[, 1]))
Lower_MNAR = c(mean(log_mmdperf_m4[, 1]) -  1.96*sd(log_mmdperf_m4[, 1]), mean(log_mmdperf_m8[, 1]) - 1.96*sd(log_mmdperf_m8[, 1]), mean(log_mmdperf_m16[, 1]) -  1.96*sd(log_mmdperf_m16[, 1]))

## MNAR + confounding

load("/Users/kyuseongchoi/Desktop/DNN/addsim/m = 4 ,n = 10 ,N = 70 ,dim = 2 ,e_T = 0.5 ,e_N = 0.5 confounding.Rdata")

sigma_traj_m4 = Collective_Info[[1]]
eta_star_pool_m4 = Collective_Info[[2]]
nbhd_Size_m4 = Collective_Info[[3]]
mmdperf_m4 = Collective_Info[[4]]
meanperf_m4 = Collective_Info[[5]]
trueinfo_m4 = Collective_Info[[6]]

load("/Users/kyuseongchoi/Desktop/DNN/addsim/m = 8 ,n = 10 ,N = 70 ,dim = 2 ,e_T = 0.5 ,e_N = 0.5 confounding.Rdata")

sigma_traj_m8 = Collective_Info[[1]]
eta_star_pool_m8 = Collective_Info[[2]]
nbhd_Size_m8 = Collective_Info[[3]]
mmdperf_m8 = Collective_Info[[4]]
meanperf_m8 = Collective_Info[[5]]
trueinfo_m8 = Collective_Info[[6]]

load("/Users/kyuseongchoi/Desktop/DNN/addsim/m = 16 ,n = 10 ,N = 70 ,dim = 2 ,e_T = 0.5 ,e_N = 0.5 confounding.Rdata")

sigma_traj_m16 = Collective_Info[[1]]
eta_star_pool_m16 = Collective_Info[[2]]
nbhd_Size_m16 = Collective_Info[[3]]
mmdperf_m16 = Collective_Info[[4]]
meanperf_m16 = Collective_Info[[5]]
trueinfo_m16 = Collective_Info[[6]]


log_mmdperf_m4 = log(mmdperf_m4)
log_mmdperf_m8 = log(mmdperf_m8)
log_mmdperf_m16 = log(mmdperf_m16)


center_E_MNAR_conf = c(mean(log_mmdperf_m4[, 1]), mean(log_mmdperf_m8[, 1]), mean(log_mmdperf_m16[, 1]))
# center_E_MNAR_conf = c(mean(log_mmdperf_m4[, 2]), mean(log_mmdperf_m8[, 2]), mean(log_mmdperf_m16[, 2]), mean(log_mmdperf_m32[, 2]))
Upper_MNAR_conf = c(mean(log_mmdperf_m4[, 1]) +  1.96*sd(log_mmdperf_m4[, 1]), mean(log_mmdperf_m8[, 1]) +  1.96*sd(log_mmdperf_m8[, 1]), mean(log_mmdperf_m16[, 1]) +  1.96*sd(log_mmdperf_m16[, 1]))
Lower_MNAR_conf = c(mean(log_mmdperf_m4[, 1]) -  1.96*sd(log_mmdperf_m4[, 1]), mean(log_mmdperf_m8[, 1]) - 1.96*sd(log_mmdperf_m8[, 1]), mean(log_mmdperf_m16[, 1]) -  1.96*sd(log_mmdperf_m16[, 1]))


Group = c(rep("MCAR", 3), rep("L shape", 3), rep("L & Conf", 3))
Logdim = c(c(2:4), c(2:4), c(2:4))
Center = c(center_E_MCAR, center_E_MNAR, center_E_MNAR_conf)
Upper = c(Upper_MCAR, Upper_MNAR, Upper_MNAR_conf)
Lower = c(Lower_MCAR, Lower_MNAR, Lower_MNAR_conf)

center_ci = data.frame(Group, Logdim, Center, Upper, Lower)
colnames(center_ci) = c("Pattern", "Log row", "Log MMD", "Upper 95", "Lower 95")

pd= position_dodge(0.1)

center_ci %>%
  ggplot(aes(`Log row`, `Log MMD`, colour = Pattern, group = Pattern)) + 
  geom_errorbar(aes(ymin = `Lower 95`, ymax = `Upper 95`), width = .1, position = pd) +
  geom_line(position = pd) + 
  geom_point( aes(shape = Pattern), position = pd, size = 2.5) + theme_light(base_size = 17)



