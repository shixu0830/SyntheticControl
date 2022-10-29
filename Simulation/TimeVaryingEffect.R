###########################################################
##### Setting: Linear interactive fixed effects model #####
##### Time-varying effect: beta_t = 1+t/T             #####
##### Synthetic control methods                       #####
##### (1) Unconstrained OLS                           #####
##### (2) Proximal inference                          #####
##### shixu@umich.edu                                 #####
##### 09/26/2021                                      #####
###########################################################
rm(list=ls())
source("myfunctions_TimeVaryingEffect.R")

param_grid <- expand.grid(n.rep = 200,
                          t0 = c(60, 150, 300),
                          n.units = c(1 + 2, 1 + 6, 1 + 10),
                          dist.epsilon = c("iid", "AR"),
                          addcov = c(TRUE, FALSE),
                          batch = 1:10) 

for (job_id in 1:nrow(param_grid)) {
  n.rep <- param_grid[job_id, "n.rep"]
  t0 <- param_grid[job_id, "t0"]
  n.units <- param_grid[job_id, "n.units"]
  dist.epsilon <- param_grid[job_id, "dist.epsilon"]
  addcov <- param_grid[job_id, "addcov"]
  batch <- param_grid[job_id, "batch"]
  
  
  theta <- 1
  true.beta <- c(1, 1); n.b <- 2
  mysd <- 2
  epsilonARMAcor <- 0.1
  
  t_star <- t0
  t <- t_star + t0
  
  ## save results
  myseeds <- (batch - 1) * n.rep + 1:n.rep
  
  rslt.all <- t(sapply(myseeds,
                       function(seed) {
                         run.one(seed = seed, n.units = n.units, t = t, t0 = t0,
                                 mysd = mysd, theta = theta, true.beta = true.beta,
                                 n.b = n.b, addcov = addcov, dist.epsilon = dist.epsilon,
                                 epsilonARMAcor = epsilonARMAcor)
                       }))
  
  saveRDS(rslt.all, file = paste0("results/TV_",
                                  "epsilon_", dist.epsilon,
                                  "_cov_", addcov,
                                  "_t0_", t0,
                                  "_nunits_", n.units,
                                  "_nrep_", n.rep,
                                  "_batch", batch,
                                  ".rds"))
}