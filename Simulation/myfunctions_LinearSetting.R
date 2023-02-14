
## implement one run of simulation
run.one <- function(seed, n.units, t, t0, mysd, theta, true.beta, addcov = F,
                    dist.epsilon, epsilonARMAcor = 0.1){
  set.seed(seed)
  .env <- environment()
  
  n.ctrl <- n.units - 1 ## one treated unit
  n.Z <- floor(n.ctrl / 2) ## use half of the control units as proxies
  n.W <- n.ctrl - n.Z ## number of donors
  ind.trt <- 1 ## index of treated units
  ind.ctrl <- setdiff(1:n.units, ind.trt)
  ind.Z <- 1:n.Z ## index of proxies among the control units
  
  ## indicator of treatment
  X <<- c(rep(0, t0), rep(1, t - t0))
  
  ### gen lambda and U
  U.lambda <- generate.U(n.units, t, mysd)
  U <- U.lambda$U ## loading matrix
  U.trt.prob <- U.lambda$U.trt.prob
  lambda.t <- U.lambda$lambda.t
  
  ### generate covariates
  if(addcov == T){
    C <- replicate(n = n.units, expr = rnorm(t, sd = mysd))
  }
  
  ### generate outcome
  Y.allunits = NULL
  
  for(i in 1:n.units){
    ### first generate stationary weakly dependent error term epsilon
    if(dist.epsilon=="iid"){
      epsilon = rnorm(t, mean = 0, sd = mysd)
    }else if(dist.epsilon=="AR"){
      epsilon = as.numeric(arima.sim(n = t, list(ar = epsilonARMAcor), innov = rnorm(t)))
    }else{
      print('unrecognized dist.epsilon')
      break
    }
    ### then generate outcome for all units
    if(addcov==T){
      Yi = c(C[,i] * theta + U[i, ] %*% t(lambda.t) + epsilon)
    }else{
      Yi = c(U[i, ] %*% t(lambda.t) + epsilon)
    }
    Y.allunits = cbind(Y.allunits, Yi)
  }
  colnames(Y.allunits)=1:(n.ctrl+1)
  Y.allunits <- Y.allunits
  
  ##
  Y <<- Y.allunits[,ind.trt] + true.beta * X
  V <<- Y.allunits[,ind.ctrl]
  Z <- cbind(Y.allunits[,ind.ctrl[ind.Z]])
  W <- cbind(Y.allunits[,ind.ctrl[-ind.Z]])
  U.Y <- cbind(U[ind.trt,])
  U.Z <- t(U[ind.ctrl[ind.Z],])
  U.W <- t(U[ind.ctrl[-ind.Z],])
  if(addcov == T){
    C.Y <<- cbind(C[,ind.trt])
    C.Z <- cbind(C[,ind.ctrl[ind.Z]])
    C.W <<- cbind(C[,ind.ctrl[-ind.Z]])
  }
  #### define vcov assumption on epsilon
  if(dist.epsilon == "iid"){
    vcov.epsilon <- "iid"
  }else{
    vcov.epsilon <- "HAC"
  }
  
  
  #### SC methods
  if(addcov == T){
    fml_g <- as.formula(Y ~ -1 + V + X + C.Y + C.W, env = .env)
    fml_x <- as.formula(~ -1 + V + X + C.Y + C.W, env = .env)
    SC.model <- try(
      gmm::gmm(g = fml_g,
               x = fml_x,
               vcov = vcov.epsilon, prewhite = FALSE),
      silent=TRUE)
  } else {
    fml_g <- as.formula(Y ~ -1 + V + X, env = .env)
    fml_x <- as.formula(~ -1 + V + X, env = .env)
    SC.model <- try(
      gmm::gmm(g = fml_g,
               x = fml_x,
               vcov = vcov.epsilon, prewhite = FALSE),
      silent=TRUE)
  }
  

  
  if (class(SC.model) != "try-error"){
    SC.ate <- summary(SC.model)$coef["X","Estimate"]
    SC.se <- sqrt(vcovHAC(SC.model, prewhite = FALSE)[n.ctrl + 1, n.ctrl + 1])
    # SC.se  <- summary(SC.model)$coef["X","Std. Error"]
  } else {
    SC.ate <- SC.se <- NA
  }
  # if(addcov==T & addsynth==T){
  #   #####################################
  #   #### SC methods with constraints ####
  #   #####################################
  #   tmp=capture.output(
  #     {
  #       SC.rslt <- synth(
  #         data.prep.obj = NULL,
  #         X1=cbind(C.Y[1:t0]), X0=C[1:t0,ind.ctrl], 
  #         Z1=cbind(Y[1:t0]), Z0=Y.allunits[1:t0,ind.ctrl]
  #       )
  #     }
  #   )
  #   SC.ate.constrained = mean(cbind(Y[(t0+1):t])-Y.allunits[(t0+1):t,ind.ctrl]%*%SC.rslt$solution.w)
  # }else{
  #   SC.ate.constrained=NA
  # }
  #### NC methods
  if(addcov == T){
    ### adjust for covariates
    data <- list(X = X, Y = Y, W = W, Z = Z, C.Y = C.Y, C.W = C.W, C.Z = C.Z)
    NC2 <- NC_SC(data)
    NC.ate2 <- as.numeric(NC2["NC.ate2"])  ## 3.170
    NC.se2 <- as.numeric(NC2["NC.se2"])  ## 0.771
    NC.se2.HAC <- as.numeric(NC2["NC.se2.HAC"])  ## 0.551
    
    
    
    ### ignore covariates
    NC_ignoreC <- NC_SC_nocov(data)
    NC.ate <- as.numeric(NC_ignoreC["NC.ate2"])
    NC.se <- as.numeric(NC_ignoreC["NC.se2"])
    NC.se.HAC <- as.numeric(NC_ignoreC["NC.se2.HAC"])
  }else{
    data <- list(X = X, Y = Y, Z = Z, W = W)
    NC2 <- NC_SC_nocov(data)
    NC.ate2 <- as.numeric(NC2["NC.ate2"])
    NC.se2 <- as.numeric(NC2["NC.se2"])
    NC.se2.HAC <- as.numeric(NC2["NC.se2.HAC"])
    
    NC.ate = NC.se = NC.se.HAC = NA
  }#if(addcov==T){
  
  return(c(SC.ate = SC.ate, ##unconstrained OLS
           SC.ate2 = NA, ##constrained OLS
           NC.ate = NC.ate,##NC ignore C
           NC.ate2 = NC.ate2,##NC adjust C
           SC.se = SC.se, ##unconstrained OLS
           NC.se = NC.se, NC.se.HAC = NC.se.HAC, ##NC ignore C
           NC.se2 = NC.se2, NC.se2.HAC = NC.se2.HAC ##NC adjust C
           ))
}

gen.lambda <- function(t, mysd){
  mylambda <- log(1:t) + rnorm(t, sd = mysd)
  return(mylambda)
}

generate.U <- function(n.units, t, mysd){
  n.Z <- floor((n.units-1)/2) ## one treated unit, half of the control units are proxies
  n.W <- n.units - 1 - n.Z
  U_0 <- rep(1, n.W)
  U.mat <- matrix(0,ncol = n.W, nrow = n.W)
  diag(U.mat)<-1
  U <- rbind(U_0, U.mat, U.mat)
  U.w <- U[1 + 1:n.W,]; true.weights <- solve(t(U.w)) %*% U[1,]
  U.trt.prob <- U %*% rep(1, n.W)
  lambda.t <- cbind(replicate(n = n.W, gen.lambda(t, mysd))) ## the number of unmeasured confounders same as number of proxies
  U.trt.prob <- U.trt.prob / sum(U.trt.prob) ## inverse propability of treatment weights
  return(list(U = U, U.trt.prob = U.trt.prob,lambda.t = lambda.t))
}

# Negative control approach to synthetic control
NC_SC <- function(data, q = 10) {
  S1 <- with(data, cbind(X, W, C.Y, C.W))
  S2 <- with(data, cbind(X, Z, C.Y, C.W))
  Y <- data$Y
  spsz <- length(data$X)
  theta_est <- MASS::ginv(t(S1) %*% S2 %*% t(S2) %*% S1) %*% t(S1) %*% S2 %*% t(S2) %*% Y
  
  bg <- c(Y - S1 %*% theta_est) * S2
  bG <- - MASS::ginv(t(S2) %*% S1 / spsz)
  
  # "sandwich" variance estimate and Newey-West HAC variance
  hacOmega <- Omega <- t(bg) %*% bg / spsz
  for(i in 1:q){
    Omega_i <- t(bg[-(1:i),]) %*% bg[1:(spsz-i),] / spsz
    hacOmega <- hacOmega + (1 - i/(q+1))*(Omega_i + t(Omega_i))
  }
  Sigma <- bG %*% Omega %*% t(bG)
  hacSigma <- bG %*% hacOmega %*% t(bG)
  
  VAR <- Sigma/spsz; HACVAR <- hacSigma/spsz
  
  return(c(NC.ate2 = theta_est[1], NC.se2 = sqrt(VAR[1, 1]), 
           NC.se2.HAC = sqrt(HACVAR[1, 1])))
}

# Negative control approach to synthetic control, ignoring covariates
NC_SC_nocov <- function(data, q = 10) {
  S1 <- with(data, cbind(X, W))
  S2 <- with(data, cbind(X, Z))
  Y <- data$Y
  spsz <- length(data$X)
  theta_est <- MASS::ginv(t(S1) %*% S2 %*% t(S2) %*% S1) %*% t(S1) %*% S2 %*% t(S2) %*% Y
  
  bg <- c(Y - S1 %*% theta_est) * S2
  bG <- - MASS::ginv(t(S2) %*% S1 / spsz)
  
  # "sandwich" variance estimate and Newey-West HAC variance
  hacOmega <- Omega <- t(bg) %*% bg / spsz
  for(i in 1:q){
    Omega_i <- t(bg[-(1:i),]) %*% bg[1:(spsz-i),] / spsz
    hacOmega <- hacOmega + (1 - i/(q+1))*(Omega_i + t(Omega_i))
  }
  Sigma <- bG %*% Omega %*% t(bG)
  hacSigma <- bG %*% hacOmega %*% t(bG)
  
  VAR <- Sigma/spsz; HACVAR <- hacSigma/spsz
  
  return(c(NC.ate2 = theta_est[1], NC.se2 = sqrt(VAR[1, 1]), 
           NC.se2.HAC = sqrt(HACVAR[1, 1])))
}
