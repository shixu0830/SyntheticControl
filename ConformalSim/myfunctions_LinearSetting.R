
## implement one run of simulation
run.one <- function(seed, grid, n.units, t, t0, mysd, theta, true.beta, addcov = F,
                    dist.epsilon, epsilonARMAcor = 0.1, alpha = 0.10){
  set.seed(seed)
  .env <- environment()
  
  SC.pval <- NC.pval <- grid * NA  ## grid of p-values
  
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
    if(dist.epsilon == "iid"){
      epsilon = rnorm(t, mean = 0, sd = mysd)
    } else if (dist.epsilon == "AR"){
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
  colnames(Y.allunits) <- 1:(n.ctrl+1)
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
    C.Y <<- cbind(C[, ind.trt])
    C.V <- cbind(C[, ind.ctrl])
    C.Z <- cbind(C[, ind.ctrl[ind.Z]])
    C.W <- cbind(C[, ind.ctrl[-ind.Z]])
  }
  
  #### define vcov assumption on epsilon
  if(dist.epsilon == "iid"){
    vcov.epsilon <- "iid"
  }else{
    vcov.epsilon <- "HAC"
  }
  
  for (jj in seq_along(grid)) {
    Y1 <<- Y - c(rep(0, t0), grid[jj])
    
    #### SC methods
    if(addcov == T){
      fml_g <- as.formula(Y1 ~ -1 + V + C.Y, env = .env)
      fml_x <- as.formula(~ -1 + V + C.Y, env = .env)
      SC.model <- try(
        gmm::gmm(g = fml_g,
                 x = fml_x,
                 vcov = vcov.epsilon, prewhite = FALSE),
        silent=TRUE)
    } else {
      fml_g <- as.formula(Y1 ~ -1 + V, env = .env)
      fml_x <- as.formula(~ -1 + V, env = .env)
      SC.model <- try(
        gmm::gmm(g = fml_g,
                 x = fml_x,
                 vcov = vcov.epsilon, prewhite = FALSE),
        silent = TRUE)
    }
    
    
    if (class(SC.model) != "try-error") {
      if (addcov == TRUE) {
        SC.weights <- summary(SC.model)$coef[1:n.ctrl,"Estimate"]
        SC.uest <- summary(SC.model)$coef[n.ctrl + 1, "Estimate"]
        SC.resid <- c((Y1 - C.Y * SC.uest) - (V - C.V * SC.uest) %*% SC.weights)
        SC.stat <- abs(SC.resid[t])
        SC.pval[jj] <- 1 - mean(SC.stat > abs(SC.resid))
      } else {
        SC.weights <- summary(SC.model)$coef[,"Estimate"]
        SC.resid <- c(Y1 - V %*% SC.weights)
        SC.stat <- abs(SC.resid[t])
        SC.pval[jj] <- 1 - mean(SC.stat > abs(SC.resid))
      }
      # SC.se  <- summary(SC.model)$coef["X","Std. Error"]
    } else {
      SC.pval[jj] <- NA
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
      data <- list(X = X, Y = Y1, W = W, Z = Z, C.Y = C.Y, C.W = C.W, C.Z = C.Z)
      NC2 <- try(optim(par = rep(0, 1 + ncol(W)),
                       fn = GMMF, mrf = NC.U0,
                       data = data,
                       method = "BFGS", 
                       hessian = FALSE)$par,
                 silent=TRUE)
      if (!inherits(NC2, "try-error")) {
        NC.weights <- NC2[1:ncol(W)]
        NC.uest <- NC2[ncol(W) + 1]
        NC.resid <- c((Y1 - C.Y * NC.uest) - (W - C.W * NC.uest) %*% NC.weights)
        NC.stat <- abs(NC.resid[t])
        NC.pval[jj] <- 1 - mean(NC.stat > abs(NC.resid))
      } else {
        NC.pval[jj] <- NA
      }
    }else{
      data <- list(Y = Y1, Z = Z, W = W)
      NC2 <- try(optim(par = rep(0, ncol(W)),
                       fn = GMMF, mrf = NC.U0.nocov,
                       data = data,
                       method = "BFGS", 
                       hessian = FALSE)$par,
                 silent=TRUE)
      if (!inherits(NC2, "try-error")) {
        NC.weights <- NC2[1:ncol(W)]
        NC.resid <- c(Y1  - W %*% NC.weights)
        NC.stat <- abs(NC.resid[t])
        NC.pval[jj] <- 1 - mean(NC.stat > abs(NC.resid))
      } else {
        NC.pval[jj] <- NA
      }
    }#if(addcov==T){
  }
  
  
  return(c(SC_lb = grid[min(which(SC.pval > alpha))],
           SC_ub = grid[max(which(SC.pval > alpha))],
           NC_lb = grid[min(which(NC.pval > alpha))],
           NC_ub = grid[max(which(NC.pval > alpha))]))
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


## estimating equation to estimate the SC weights
NC.U0.nocov <- function(para, data){
  Y <- cbind(data$Y);
  Z <- cbind(data$Z)
  W <- cbind(data$W)
  omega.u <- para[1:ncol(W)]

  NC.obj <- c(Y - W %*% omega.u) * Z
  return(NC.obj)
}

NC.U0 <- function(para, data){
  Y <- cbind(data$Y); C.Y=cbind(data$C.Y)
  Z <- cbind(data$Z); C.Z=cbind(data$C.Z)
  W <- cbind(data$W); C.W=cbind(data$C.W)
  theta.u <- para[ncol(W) + 1]
  omega.u <- para[1:ncol(W)]
  Ztheta <- Z - C.Z * theta.u
  Comega <- C.Y - C.W %*% omega.u
  instrument <- cbind(
    Ztheta,
    Comega
  )
  bridge <- W %*% omega.u + Comega * theta.u
  
  NC.obj <- c(Y - bridge) * instrument
  return(NC.obj)
}


### generic functions below
library(numDeriv)

# GMM function
GMMF <- function(mrf, para, data){
  g0 <- mrf(para = para, data = data)
  g <- apply(g0, 2, mean)
  gmmf <- sum(g ^ 2)
  
  return(gmmf)
}

# Derivative of score equations
G1 <- function(bfun,para,data){
  G1 <- apply(bfun(para,data), 2, mean)
  return(G1)
}

G <- function(bfun, para, data){
  G <- jacobian(func = G1, bfun = bfun, x = para, data = data)
  return(G)
}
