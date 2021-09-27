library(Synth)
run.one = function(rep,n.units,t,t0){
  set.seed(rep)
  n.ctrl <<- n.units-1
  n.Z <<- floor(n.ctrl/2)
  n.W <<- n.ctrl-n.Z
  ind.trt <<- 1
  ind.ctrl <<- setdiff(1:n.units,ind.trt)
  ind.Z <<- 1:n.Z
  ### gen lambda and U
  U.lambda = generate.U(n.units,t)
  U = U.lambda$U
  U.trt.prob = U.lambda$U.trt.prob
  lambda.t = U.lambda$lambda.t
  ### generate covariates
  if(addcov==T){
    C=replicate(n=n.units,expr=runif(t,min=0,max=1))
  }
  ### generate outcome
  Y.allunits = NULL
  for(i in 1:n.units){
    ### generate outcome for all units
    if(i==ind.trt){### this is Y derived using moment generating function
      if(addcov==T){
        mean.Y = exp(as.numeric(C[,i]*theta + U[i,]%*%t(lambda.t) + true.beta*X))
      }else{
        mean.Y = exp(as.numeric(U[i,]%*%t(lambda.t) + true.beta*X))
      }
      Yi = rpois( n=t, lambda=mean.Y )
    }else{ ### these are W
      if(addcov==T){
        mean.Y = as.numeric(C[,i]*theta + U[i,]%*%t(lambda.t))/(exp(1)-1)
      }else{
        mean.Y = as.numeric(U[i,]%*%t(lambda.t))/(exp(1)-1)
      }
      Yi = rpois( n=t, lambda=mean.Y )
    }
    Y.allunits = cbind(Y.allunits,Yi)
  }#for(i in 1:n.units){
  colnames(Y.allunits)=1:(n.ctrl+1)
  Y.allunits <<- Y.allunits
  ##
  Y <<- Y.allunits[,ind.trt]
  Z <<- cbind(Y.allunits[,ind.ctrl[ind.Z]])
  W <<- cbind(Y.allunits[,ind.ctrl[-ind.Z]])
  U.Y <<- cbind(U[ind.trt,])
  U.Z <<- t(U[ind.ctrl[ind.Z],])
  U.W <<- t(U[ind.ctrl[-ind.Z],])
  if(addcov==T){
    C.Y <<- cbind(C[,ind.trt])
    C.Z <<- cbind(C[,ind.ctrl[ind.Z]])
    C.W <<- cbind(C[,ind.ctrl[-ind.Z]])
  }
  #### define vcov assumption on epsilon
  if(dist.epsilon=="iid"){
    vcov.epsilon="iid"
  }else{
    vcov.epsilon="HAC"
  }
  # #### SC methods
  # if(addcov==T){
  #   SC.model = try(
  #     gmm::gmm(g=Y~-1+Y.allunits[,ind.ctrl]+X+C.Y,
  #              x=~-1+Y.allunits[,ind.ctrl]+X+C.Y,
  #              vcov=vcov.epsilon),
  #     silent=TRUE)
  # }else{
  #   SC.model = try(
  #     gmm::gmm(g=Y~-1+Y.allunits[,ind.ctrl]+X,
  #              x=~-1+Y.allunits[,ind.ctrl]+X,
  #              vcov=vcov.epsilon),
  #     silent=TRUE)
  # }
  # if(!inherits(SC.model, "try-error")){
  #   SC.ate = summary(SC.model)$coef["X","Estimate"]
  #   SC.se  = summary(SC.model)$coef["X","Std. Error"]
  # }else{
  #   SC.ate = SC.se  = NA
  # }
  #### NC methods
  if(addcov==T){
    ### adjust for covariates
    data = list(X=X,Y=Y,W=W,Z=Z,C.Y=C.Y,C.W=C.W,C.Z=C.Z)
    NC2 = mygmm(NC.U.func=NC.U0,data=data,loc.X=1+ncol(W)+1)
    NC.ate2 = as.numeric(NC2["NC.ate2"])
    NC.se2 = as.numeric(NC2["NC.se2"])
    NC.se2.HAC = as.numeric(NC2["NC.se2.HAC"])
    
    ### ignore covariates
    NC_ignoreC = mygmm(NC.U.func=NC.U0.nocov,data=data,loc.X=1+ncol(W))
    NC.ate = as.numeric(NC_ignoreC["NC.ate2"])
    NC.se = as.numeric(NC_ignoreC["NC.se2"])
    NC.se.HAC = as.numeric(NC_ignoreC["NC.se2.HAC"])
  }else{
    data = list(X=X,Y=Y,Z=Z,W=W)
    NC2 = mygmm(NC.U.func=NC.U0.nocov,data=data,loc.X=1+ncol(W))
    NC.ate2 = as.numeric(NC2["NC.ate2"])
    NC.se2 = as.numeric(NC2["NC.se2"])
    NC.se2.HAC = as.numeric(NC2["NC.se2.HAC"])
    
    NC.ate = NC.se = NC.se.HAC = NA
  }#if(addcov==T){
  if(rep%%20==0){
    print(paste0(n.units,"_",t0,"_",rep,"_bias_",
                 paste(round((apply(ate.all,2,mean)+c(2,2,2,2,0,0,0,0,0))*1000),collapse="_")
    ))
  }
  return(c(SC.ate=NA, ##unconstrained OLS
           SE.ate2=NA, ##constrained OLS
           NC.ate=NC.ate,##NC ignore C
           NC.ate2=NC.ate2,##NC adjust C
           SC.se =NA, ##unconstrained OLS
           NC.se =NC.se,NC.se.HAC=NC.se.HAC, ##NC ignore C
           NC.se2=NC.se2,NC.se2.HAC=NC.se2.HAC##NC adjust C
  ))
}

gen.lambda = function(dist.lambda){
  mylambda = runif(t,min=0,max=1)
  return(mylambda)
}

generate.U = function(n.units,t){
  n.Z = floor((n.units-1)/2)
  n.W = n.units-1-n.Z
  U_0 = rep(1,n.W)
  U.mat = matrix(0,ncol=n.W,nrow=n.W)
  diag(U.mat)=1
  U = rbind(U_0,U.mat,U.mat)
  U.w = U[1+1:n.W,]; true.weights=solve(t(U.w))%*%U[1,]
  U.trt.prob=U%*%rep(1,n.W)
  lambda.t = cbind(replicate(n=n.W,gen.lambda(dist.lambda)))
  U.trt.prob = U.trt.prob/sum(U.trt.prob)
  return(list(U=U,U.trt.prob=U.trt.prob,lambda.t=lambda.t))
}

mygmm = function(NC.U.func,data,loc.X){
  NC.weights = try(
    optim(par = rep(0,loc.X),
          fn = GMMF, mrf = NC.U.func,
          data=data,
          method = "BFGS", hessian = FALSE)$par,
    silent=TRUE)
  var_est = try(HAC_VAREST(bfun=NC.U.func,para=NC.weights,data=data,q=10),silent=TRUE)
  if((!inherits(NC.weights, "try-error"))&
     (!inherits(var_est, "try-error"))){
    NC.se2 <- diag(var_est$var)
    NC.se2.HAC <- diag(var_est$hacvar)
    NC.ate2 = NC.weights[loc.X]
    NC.se2 = sqrt(NC.se2[loc.X])
    NC.se2.HAC = sqrt(NC.se2.HAC[loc.X])
  }else{
    NC.ate2 = NC.se2 = NC.se2.HAC = NA
  }
  return(c(NC.ate2=NC.ate2, NC.se2=NC.se2, NC.se2.HAC=NC.se2.HAC))
}

NC.U0.nocov = function(para,data){
  X=cbind(data$X)
  Y=cbind(data$Y)
  Z=cbind(data$Z)
  W=cbind(data$W)
  omega.u = para[1:ncol(W)]
  beta = para[ncol(W)+1]
  instrument = cbind(
    Z, #c(1-X)*Z, #if pre and post two-stage
    X
  )
  bridge = exp(
    X*beta + W%*%omega.u
  )
  NC.obj=c(  Y-bridge   )*(instrument)
  return(NC.obj)
}

NC.U0 = function(para,data){
  X=cbind(data$X)
  Y=cbind(data$Y); C.Y=cbind(data$C.Y)
  Z=cbind(data$Z); C.Z=cbind(data$C.Z)
  W=cbind(data$W); C.W=cbind(data$C.W)
  theta.u = para[1]
  omega.u = para[1+1:ncol(W)]
  beta = para[1+ncol(W)+1]
  Comega = C.Y-C.W%*%omega.u
  instrument = cbind(
    Z,
    Comega,
    X
  )
  bridge = exp(
    X*beta + W%*%omega.u + Comega*theta.u
  )
  NC.obj=c(  Y-bridge  )*(instrument)
  return(NC.obj)
}

### generic functions below
library(numDeriv)

# GMM function
GMMF <- function(mrf,para,data){
  g0 <- mrf(para=para,data=data)
  g <- apply(g0,2,mean)
  gmmf <- sum(g^2)
  
  return(gmmf)
}

# Derivative of score equations
G1 <- function(bfun,para,data){
  G1 <- apply(bfun(para,data),2,mean)
  return(G1)
}

G <- function(bfun,para,data){
  G <- jacobian(func=G1,bfun=bfun,x=para,data=data)
  return(G)
}

# Variance estimation
VAREST <- function(bfun,para,data){
  bG <- solve(G(bfun,para,data))
  bg <- bfun(para,data)
  spsz <- dim(bg)[1]
  Omega <- t(bg)%*%bg/spsz
  Sigma <- bG%*%Omega%*%t(bG)
  
  return(Sigma/spsz)
}

# Newey-West 1987 variance estimator for serially correlated data
HAC_VAREST <- function(bfun,para,q,data){
  bG <- solve(G(bfun,para,data))
  bg <- bfun(para,data)
  spsz <- dim(bg)[1]
  hacOmega <- Omega <- t(bg)%*%bg/spsz
  for(i in 1:q){
    Omega_i <- t(bg[-(1:i),])%*%bg[1:(spsz-i),]/spsz
    hacOmega <- hacOmega + (1 - i/(q+1))*(Omega_i + t(Omega_i))
  }
  Sigma <- bG%*%Omega%*%t(bG)
  hacSigma <- bG%*%hacOmega%*%t(bG)
  
  return(list(var=Sigma/spsz, hacvar=hacSigma/spsz))
}
