
## implement one run of simulation
run.one <- function(seed, grid, n.units, t, t0, mysd, theta, true.beta, addcov = F,
                    dist.epsilon, epsilonARMAcor = 0.1, alpha = 0.10,
                    output = "both") {
  set.seed(seed)
  .env <- environment()
  
  # SC.pval <- NC.pval <- grid * NA  ## grid of p-values
  n.ctrl <- n.units - 1 ## one treated unit
  n.Z <- floor(n.ctrl / 2) ## use half of the control units as proxies
  n.W <- n.ctrl - n.Z ## number of donors
  ind.trt <- 1 ## index of treated units
  ind.ctrl <- setdiff(1:n.units, ind.trt)
  ind.Z <- 1:n.Z ## index of proxies among the control units
  
  ## indicator of treatment
  X <<- c(rep(0, t0), 1)
  
  ### gen lambda and U
  U.lambda <- generate.U(n.units, t0 + 1, mysd)
  U <- U.lambda$U ## loading matrix
  U.trt.prob <- U.lambda$U.trt.prob
  lambda.t <- U.lambda$lambda.t
  ### generate covariates
  if(addcov == T){
    C <- replicate(n = n.units, expr = rnorm(t0 + 1, sd = mysd))
  }
  ### generate outcome
  Y.allunits = NULL
  for(i in 1:n.units){
    ### first generate stationary weakly dependent error term epsilon
    if(dist.epsilon == "iid"){
      epsilon = rnorm(t0 + 1, mean = 0, sd = mysd)
    } else if (dist.epsilon == "AR"){
      epsilon = as.numeric(arima.sim(n = t0 + 1, list(ar = epsilonARMAcor), innov = rnorm(t0 + 1)))
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
  
  if(addcov == T){
    ### adjust for covariates
    data <- list(Y = Y, W = W, Z = Z, C.Y = C.Y, C.W = C.W, C.Z = C.Z)
    
  } else {
    data <- list(Y = Y, Z = Z, W = W)
  }#if(addcov==T){
  
  ### NC method
  ## crude search
  crude_CI <- ConformalPointwiseCI(data = data, t0 = t0, grid = seq(-18, 22, 0.1),
                                   addcov = addcov, output = "both", 
                                   method = "NC", alpha = alpha)
  ## find search
  NC_CI_lb <- 
    ConformalPointwiseCI(data = data, t0 = t0, 
                         grid = seq(crude_CI[1] - 0.09, crude_CI[1], by = 0.01),
                         addcov = addcov, output = "lb", 
                         method = "NC", alpha = alpha)
  
  NC_CI_ub <- 
    ConformalPointwiseCI(data = data, t0 = t0, 
                         grid = seq(crude_CI[2], crude_CI[2] + 0.09, by = 0.01),
                         addcov = addcov, output = "ub", 
                         method = "NC", alpha = alpha)
  
  
  ### SC method
  ## crude search
  crude_CI <- ConformalPointwiseCI(data = data, t0 = t0, grid = seq(-18, 22, 0.1),
                                   addcov = addcov, output = "both", 
                                   method = "SC", alpha = alpha)
  ## find search
  SC_CI_lb <- 
    ConformalPointwiseCI(data = data, t0 = t0, 
                         grid = seq(crude_CI[1] - 0.09, crude_CI[1], by = 0.01),
                         addcov = addcov, output = "lb", 
                         method = "SC", alpha = alpha)
  SC_CI_ub <- 
    ConformalPointwiseCI(data = data, t0 = t0, 
                         grid = seq(crude_CI[2], crude_CI[2] + 0.09, by = 0.01),
                         addcov = addcov, output = "ub", 
                         method = "SC", alpha = alpha)
  
  return(c(NC_lb = NC_CI_lb, NC_ub = NC_CI_ub,
           SC_lb = SC_CI_lb, SC_ub = SC_CI_ub))
}

## Conformal pointwise confidence intervals
ConformalPointwiseCI <- 
  function(data, t0, grid, addcov, output = "both", alpha = 0.1, method = "NC") {
    X <- cbind(data$X); Y <- cbind(data$Y)
    Z <- cbind(data$Z); W <<- cbind(data$W)
    
    if (addcov == T) {
      C.Y <<- cbind(data$C.Y); C.Z <- cbind(data$C.Z)
      C.W <- cbind(data$C.W)
    }
    X1 <- c(rep(0, t0), 1)
    
    pval <- NA * grid
    for (bb in seq_along(grid)) {
      Y1 <- cbind(c(Y[1:t0], Y[t0 + 1] - grid[bb]))
      
      if (addcov == TRUE) {
        if (method == "NC") {
          dat.h0 <- list(X = X, Y = Y1, W = W, Z = Z,
                         C.Y = C.Y, C.W = C.W,C.Z = C.Z)
          NC2 <- NC.U0(dat.h0)
          if (!inherits(NC2, "try-error")) {
            NC.weights <- NC2[1:ncol(W)]
            NC.uY.est <- NC2[ncol(W) + 1]
            NC.uW.est <- NC2[ncol(W) + 1 + 1:ncol(W)]
            NC.resid <- c(Y1 - W %*% NC.weights - C.Y %*% NC.uY.est - C.W %*% NC.uW.est)
            NC.stat <- abs(NC.resid[t0 + 1])
            pval[bb] <- 1 - mean(NC.stat > abs(NC.resid))
          } else {
            pval[bb] <- NA
          }
        } else if (method == "SC") {
          dat.h0 <- list(X = X, Y = Y1, V = cbind(W, Z),
                         C.Y = C.Y, C.V = cbind(C.W, C.Z))
          SC2 <- SC.U0(dat.h0)
          if (!inherits(SC2, "try-error")) {
            SC.weights <- SC2[1:ncol(dat.h0$V)]
            SC.uY.est <- SC2[ncol(dat.h0$V) + 1]
            SC.uV.est <- SC2[ncol(dat.h0$V) + 1 + 1:ncol(dat.h0$V)]
            SC.resid <- with(dat.h0,
                             c(Y - V %*% SC.weights - C.Y %*% SC.uY.est - 
                                 C.V %*% SC.uV.est))
            SC.stat <- abs(SC.resid[t0 + 1])
            pval[bb] <- 1 - mean(SC.stat > abs(SC.resid))
          } else {
            pval[bb] <- NA
          }
        }
      } else {
        if (method == "NC") {
          dat.h0 <- list(Y = Y1, Z = Z, W = W)
          NC2 <- NC.U0.nocov(dat.h0)
          if (!inherits(NC2, "try-error")) {
            NC.weights <- NC2[1:ncol(W)]
            NC.resid <- c(Y1  - W %*% NC.weights)
            NC.stat <- abs(NC.resid[t0 + 1])
            pval[bb] <- 1 - mean(NC.stat > abs(NC.resid))
          } else {
            pval[bb] <- NA
          }
        } else if (method == "SC") {
          dat.h0 <- list(Y = Y1, V = cbind(W, Z))
          SC2 <- SC.U0.nocov(dat.h0)
          if (!inherits(SC2, "try-error")) {
            SC.weights <- SC2[1:ncol(dat.h0$V)]
            SC.resid <- with(dat.h0, c(Y  - V %*% SC.weights))
            SC.stat <- abs(SC.resid[t0 + 1])
            pval[bb] <- 1 - mean(SC.stat > abs(SC.resid))
          } else {
            pval[bb] <- NA
          }
        }
      }
    }
    # find upper bound
    if (output == "both") {
      #  SC_lb <- grid[min(which(SC.pval > alpha))]
      #  SC_ub <- grid[max(which(SC.pval > alpha))]
      lb_id <- max(which(pval[1:which.max(pval)] < alpha), na.rm = T) + 1
      ub_id <- which.max(pval) + min(which(pval[which.max(pval):length(pval)] < alpha), na.rm = T) - 2
      
      lb <- grid[lb_id]
      ub <- grid[ub_id]
      
      return(c(lb, ub))
      
    } else if (output == "lb") {
      lb_id <- min(which(pval >= alpha), na.rm = T) 
      lb <- grid[lb_id]
      return(lb)
    } else if (output == "ub") {
      ub_id <- max(which(pval >= alpha), na.rm = T)
      ub <- grid[ub_id]
      return(ub)
    }
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
NC.U0.nocov <- function(data){
  Y <- cbind(data$Y);
  Z <- cbind(data$Z)
  W <- cbind(data$W)
  
  NC_est <- c(MASS::ginv(t(W) %*% Z %*% t(Z) %*% W) %*% t(W) %*% Z %*% t(Z) %*% Y)
  return(NC_est)
}

NC.U0 <- function(data){
  Y <- cbind(data$Y)
  S1 <- with(data, cbind(W, C.Y, C.W))
  S2 <- with(data, cbind(W, C.Y, C.W))           
  
  NC_est <- c(MASS::ginv(t(S1) %*% S2 %*% t(S2) %*% S1) %*% t(S1) %*% S2 %*% t(S2) %*% Y)
  return(NC_est)
}

## estimating equation to estimate the SC weights
SC.U0.nocov <- function(data){
  Y <- cbind(data$Y);
  V <- cbind(data$V)
  
  SC_est <- c(MASS::ginv(t(V) %*% V) %*% t(V) %*% Y)
  return(SC_est)
}

SC.U0 <- function(data){
  Y <- cbind(data$Y)
  S <- with(data, cbind(V, C.Y, C.V))
  
  
  SC_est <- c(MASS::ginv(t(S) %*% S) %*% t(S) %*% Y)
  return(SC_est)
}
