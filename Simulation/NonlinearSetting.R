###########################################################
##### Setting: nonlinear models                       #####
##### Synthetic control methods                       #####
##### (1) Proximal inference                          #####
##### shixu@umich.edu                                 #####
##### 09/26/2021                                      #####
###########################################################
rm(list=ls())
library(Synth)
library(sandwich)

source("myfunctions_NonlinearSetting.R")

param_grid <- expand.grid(n.rep = 200,
                          t0 = c(50),
                          n.units = c(1 + 2, 1 + 10, 1 + 20),
                          dist.epsilon = c("iid"),
                          addcov = c(TRUE, FALSE),
                          batch = 1:10)

for (job_id in 1:nrow(param_grid)) {
  
  n.rep <- param_grid[job_id, "n.rep"]
  t0 <- param_grid[job_id, "t0"]
  n.units <- param_grid[job_id, "n.units"]
  dist.epsilon <- param_grid[job_id, "dist.epsilon"]
  addcov <- param_grid[job_id, "addcov"]
  batch <- param_grid[job_id, "batch"]
  
  
  theta <- 0.1
  true.beta <- -2
  
  
  t_star <- t0
  t <- t_star + t0
  
  ## save results
  myseeds <- (batch - 1) * n.rep + 1:n.rep
  
  rslt.all <- t(sapply(myseeds,
                       function(seed) {
                         run.one(seed = seed, n.units = n.units, t = t, t0 = t0,
                                 theta = theta, true.beta = true.beta,
                                 addcov = addcov, dist.epsilon = dist.epsilon)
                       }))
  
  
  
  saveRDS(rslt.all, file = paste0("results/NP",
                                  "_cov_",addcov, 
                                  "_t0_",t0,
                                  "_nunits_", n.units,
                                  "_nrep_", n.rep,
                                  "_batch_", batch,
                                  ".rds"))
}