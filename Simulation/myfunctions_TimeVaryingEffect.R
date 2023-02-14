run.one <- function(seed, n.units, t, t0, theta, true.beta, n.b,
                    mysd, dist.epsilon, epsilonARMAcor, addcov){
  set.seed(seed)
  n.ctrl <- n.units - 1
  n.Z <- floor(n.ctrl / 2)
  n.W <- n.ctrl - n.Z
  ind.trt <- 1
  ind.ctrl <- setdiff(1:n.units, ind.trt)
  ind.Z <- 1:n.Z
  
  ## treatment indicator
  X <<- c(rep(0, t0), rep(1, t - t0))
  Xt <<- X * (1:t) / t
  
  ### gen lambda and U
  U.lambda <- generate.U(n.units, t, mysd)
  U <- U.lambda$U
  U.trt.prob <- U.lambda$U.trt.prob
  lambda.t <- U.lambda$lambda.t
  ### generate covariates
  if(addcov == T){
    C <- replicate(n = n.units, expr = rnorm(t, sd = mysd))
  }
  ### generate outcome
  Y.allunits <- NULL
  for(i in 1:n.units){
    ### first generate stationary weakly dependent error term epsilon
    if(dist.epsilon=="iid"){
      epsilon <- rnorm(t, mean = 0, sd = mysd)
    } else if (dist.epsilon =="AR"){
      epsilon <- as.numeric(arima.sim(n = t, list(ar = epsilonARMAcor), innov = rnorm(t)))
    } else {
      print('unrecognized dist.epsilon')
      break
    }
    ### then generate outcome for all units
    if(addcov == T){
      Yi <- c(C[,i] * theta + U[i,] %*% t(lambda.t) + epsilon)
    } else {
      Yi <- c(U[i,] %*% t(lambda.t) + epsilon)
    }
    Y.allunits <- cbind(Y.allunits, Yi)
  }
  colnames(Y.allunits) = 1:(n.ctrl + 1)
  Y.allunits <<- Y.allunits
  ##
  Y <<- Y.allunits[,ind.trt] + true.beta[1] * X + true.beta[2] * Xt
  V <<- Y.allunits[, ind.ctrl]
  Z <<- cbind(Y.allunits[, ind.ctrl[ind.Z]])
  W <<- cbind(Y.allunits[, ind.ctrl[-ind.Z]])
  U.Y <<- cbind(U[ind.trt, ])
  U.Z <<- t(U[ind.ctrl[ind.Z], ])
  U.W <<- t(U[ind.ctrl[-ind.Z], ])
  if(addcov == T){
    C.Y <<- cbind(C[, ind.trt])
    C.Z <<- cbind(C[, ind.ctrl[ind.Z]])
    C.W <<- cbind(C[, ind.ctrl[-ind.Z]])
  }
  #### define vcov assumption on epsilon
  if (dist.epsilon == "iid"){
    vcov.epsilon <- "iid"
  } else {
    vcov.epsilon <- "HAC"
  }
  #### SC methods
  if (addcov == T){
    SC.model <- try(
      gmm::gmm(g = Y ~ -1 + V + X + Xt + C.Y,
               x = ~ -1 + V + X + Xt + C.Y,
               vcov = vcov.epsilon),
      silent = TRUE)
  } else {
    SC.model <- try(
      gmm::gmm(g = Y~ -1 + V + X + Xt,
               x = ~ -1 + V + X + Xt,
               vcov = vcov.epsilon),
      silent = TRUE)
  }
  
  if(!inherits(SC.model, "try-error")){
    SC.ate <- summary(SC.model)$coef[grep("X", names(coef(SC.model))), "Estimate"]
    SC.se  <- summary(SC.model)$coef[grep("X",names(coef(SC.model))), "Std. Error"]
  } else {
    SC.ate <- SC.se <- rep(NA, n.b)
  }
  SC.ate.constrained <- rep(NA,n.b)
  #### NC methods
  if(addcov == T){
    ### adjusting covariates
    data <- list(X = X, Xt = Xt, Y = Y, W = W, Z = Z, C.Y = C.Y, C.W = C.W, C.Z = C.Z)
    NC2 <- NC_SC(data)
    NC.ate2 <- as.numeric(NC2[c("NC.ate21", "NC.ate22")])
    NC.se2 <- as.numeric(NC2[c("NC.se21", "NC.se22")])
    NC.se2.HAC <- as.numeric(NC2[c("NC.se2.HAC1", "NC.se2.HAC2")])
    
    ## ignoring covariates
    data <- list(X = X, Xt = Xt, Y = Y, W = W, Z = Z)
    NC2 <- NC_SC_nocov(data)
    NC.ate1 <- as.numeric(NC2[c("NC.ate21", "NC.ate22")])
    NC.se1 <- as.numeric(NC2[c("NC.se21", "NC.se22")])
    NC.se1.HAC <- as.numeric(NC2[c("NC.se2.HAC1", "NC.se2.HAC2")])
    
  } else {
    data <- list(X = X, Xt = Xt, Y = Y, W = W, Z = Z)
    NC2 <- NC_SC_nocov(data)
    NC.ate2 <- as.numeric(NC2[c("NC.ate21", "NC.ate22")])
    NC.se2 <- as.numeric(NC2[c("NC.se21", "NC.se22")])
    NC.se2.HAC <- as.numeric(NC2[c("NC.se2.HAC1", "NC.se2.HAC2")])
    
    NC.ate1 <- NC.se1 <- NC.se1.HAC <- rep(NA, n.b)
  }#if(addcov==T){
  
  return(c(SC.ate = SC.ate, ## unconstrained OLS
           NC.ate1 = NC.ate1,##NC ignore C
           NC.ate2 = NC.ate2,##NC adjust C
           SC.se = SC.se, ## unconstrained OLS
           NC.se1 = NC.se1, NC.se1.HAC = NC.se1.HAC, ##NC ignore C
           NC.se2 = NC.se2, NC.se2.HAC = NC.se2.HAC##NC adjust C
  ))
}#run.one

gen.lambda <- function(t, mysd){
  mylambda <- log(1:t) + rnorm(t, sd = mysd)
  return(mylambda)
}

generate.U <- function(n.units, t, mysd){
  n.Z <- floor((n.units - 1) / 2)
  n.W <- n.units - 1 - n.Z
  U_0 <- rep(1, n.W)
  U.mat <- matrix(0, ncol = n.W, nrow = n.W)
  diag(U.mat) <- 1
  U <- rbind(U_0, U.mat, U.mat)
  U.w <- U[1+1:n.W,]; true.weights <- solve(t(U.w)) %*% U[1,]
  U.trt.prob <- U %*% rep(1, n.W)
  lambda.t <- cbind(replicate(n = n.W, gen.lambda(t, mysd)))
  U.trt.prob <- U.trt.prob / sum(U.trt.prob)
  return(list(U = U, U.trt.prob = U.trt.prob, lambda.t = lambda.t))
}

# Negative control approach to synthetic control
NC_SC <- function(data, q = 10) {
  S1 <- with(data, cbind(X, Xt, W, C.Y, C.W))
  S2 <- with(data, cbind(X, Xt, Z, C.Y, C.W))
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
  
  return(c(NC.ate2 = theta_est[1:2], NC.se2 = sqrt(c(VAR[1, 1], VAR[2, 2])), 
           NC.se2.HAC = sqrt(c(HACVAR[1, 1], c(HACVAR[2, 2])))))
}

# Negative control approach to synthetic control, ignoring covariates
NC_SC_nocov <- function(data, q = 10) {
  S1 <- with(data, cbind(X, Xt, W))
  S2 <- with(data, cbind(X, Xt, Z))
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
  
  return(c(NC.ate2 = theta_est[1:2], NC.se2 = sqrt(c(VAR[1, 1], VAR[2, 2])), 
           NC.se2.HAC = sqrt(c(HACVAR[1, 1], c(HACVAR[2, 2])))))
}
