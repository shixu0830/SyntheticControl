###########################################################
##### Setting: Linear interactive fixed effects model #####
##### Synthetic control methods                       #####
##### (1) Constrained OLS (Abadie)                    #####
##### (2) Unconstrained OLS                           #####
##### (3) Proximal inference                          #####
##### shixu@umich.edu                                 #####
##### 09/26/2021                                      #####
###########################################################
rm(list=ls())
.libPaths("/home/qijunli/R_libraries")

library(gmm, lib.loc = c(.libPaths(), "/home/qijunli/R_libraries"))
library(sandwich, lib.loc = c(.libPaths(), "/home/qijunli/R_libraries"))
library(parallel, lib.loc = c(.libPaths(), "/home/qijunli/R_libraries"))

source("myfunctions_LinearSetting.R")
param_grid <- expand.grid(n.rep = 200,
                          t0 = c(50, 100, 200),
                          n.units = c(1 + 2, 1 + 6, 1 + 10),
                          dist.epsilon = c("iid", "AR"),
                          addcov = c(TRUE, FALSE),
                          alpha = c(0.1, 0.2),
                          batch = 1:10) 


for (job_id in 1:nrow(param_grid)) { 
#for (job_id in 1:10) { 
  print(job_id)
  n.rep <- param_grid[job_id, "n.rep"]
  t0 <- param_grid[job_id, "t0"]
  n.units <- param_grid[job_id, "n.units"]
  dist.epsilon <- param_grid[job_id, "dist.epsilon"]
  addcov <- param_grid[job_id, "addcov"]
  batch <- param_grid[job_id, "batch"]
  alpha <- param_grid[job_id, "alpha"]
  
  theta <- 1
  true.beta <- 2
  mysd <- 1
  epsilonARMAcor <- 0.1
  
  t_star <- 1
  t <- t_star + t0
  
  ## save results
  myseeds <- (batch - 1) * n.rep + 1:n.rep
  
  rslt.all <- 
    parallel::mclapply(myseeds,
                       function(seed) {
                         run.one(seed = seed, grid = grid, n.units = n.units, t = t, t0 = t0,
                                 mysd = mysd, theta = theta, true.beta = true.beta,
                                 addcov = addcov, dist.epsilon = dist.epsilon,
                                 epsilonARMAcor = epsilonARMAcor,
                                 alpha = alpha)
                       }, mc.cores = parallel::detectCores() - 2)
  
  saveRDS(rslt.all, file = paste0("results/LM_",
                                  "epsilon_", dist.epsilon,
                                  "_cov_", addcov,
                                  "_t0_", t0,
                                  "_nunits_", n.units,
                                  "_nrep_", n.rep,
                                  "_alpha_", alpha,
                                  "_batch_", batch,
                                  ".rds"))
}
