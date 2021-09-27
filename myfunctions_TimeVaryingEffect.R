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
    C=replicate(n=n.units,expr=rnorm(t,sd=mysd))
  }
  ### generate outcome
  Y.allunits = NULL
  for(i in 1:n.units){
    ### first generate stationary weakly dependent error term epsilon
    if(dist.epsilon=="iid"){
      epsilon = rnorm(t,mean=0,sd=2)
    }else if(dist.epsilon=="AR"){
      epsilon = as.numeric(arima.sim(n = t, list(ar = epsilonARMAcor), innov=rnorm(t)))
    }else{
      print('unrecognized dist.epsilon')
      break
    }
    ### then generate outcome for all units
    if(addcov==T){
      Yi = c(C[,i]*theta + U[i,]%*%t(lambda.t) + epsilon)
    }else{
      Yi = c(U[i,]%*%t(lambda.t) + epsilon)
    }
    Y.allunits = cbind(Y.allunits,Yi)
  }
  colnames(Y.allunits)=1:(n.ctrl+1)
  Y.allunits <<- Y.allunits
  ##
  Y <<- Y.allunits[,ind.trt]+true.beta[1]*X + true.beta[2]*Xt
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
  #### SC methods
  if(addcov==T){
    SC.model = try(
      gmm::gmm(g=Y~-1+Y.allunits[,ind.ctrl]+Xt+X+C.Y,
               x=~-1+Y.allunits[,ind.ctrl]+Xt+X+C.Y,
               vcov=vcov.epsilon),
      silent=TRUE)
  }else{
    SC.model = try(
      gmm::gmm(g=Y~-1+Y.allunits[,ind.ctrl]+Xt+X,
               x=~-1+Y.allunits[,ind.ctrl]+Xt+X,
               vcov=vcov.epsilon),
      silent=TRUE)
  }
  if(!inherits(SC.model, "try-error")){
    SC.ate = summary(SC.model)$coef[grep("X",names(coef(SC.model))),"Estimate"]
    SC.se  = summary(SC.model)$coef[grep("X",names(coef(SC.model))),"Std. Error"]
  }else{
    SC.ate = SC.se  = rep(NA,n.b)
  }
  SC.ate.constrained=rep(NA,n.b)
  #### NC methods
  if(addcov==T){
    ### adjust for covariates
    NC.model = try(
      gmm::gmm(g=Y~-1+W+Xt+X+C.Y,
               x=~-1+Z+Xt+X+C.Y,
               wmatrix="ident",
               vcov="iid"),
      silent=TRUE)
    NC.ate2 = summary(NC.model)$coef[grep("X",names(coef(NC.model))),"Estimate"]
    NC.se2  = summary(NC.model)$coef[grep("X",names(coef(NC.model))),"Std. Error"]
    NC.model = try(
      gmm::gmm(g=Y~-1+W+Xt+X+C.Y,
               x=~-1+Z+Xt+X+C.Y,
               wmatrix="ident",
               vcov="HAC"),
      silent=TRUE)
    NC.se2.HAC  = summary(NC.model)$coef[grep("X",names(coef(NC.model))),"Std. Error"]
    
    ### ignore covariates
    NC.model = try(
      gmm::gmm(g=Y~-1+W+Xt+X,
               x=~-1+Z+Xt+X,
               wmatrix="ident",
               vcov="iid"),
      silent=TRUE)
    NC.ate = summary(NC.model)$coef[grep("X",names(coef(NC.model))),"Estimate"]
    NC.se  = summary(NC.model)$coef[grep("X",names(coef(NC.model))),"Std. Error"]
    NC.model = try(
      gmm::gmm(g=Y~-1+W+Xt+X,
               x=~-1+Z+Xt+X,
               wmatrix="ident",
               vcov="HAC"),
      silent=TRUE)
    NC.se.HAC  = summary(NC.model)$coef[grep("X",names(coef(NC.model))),"Std. Error"]
    
    ### by hand GMM
    data = list(X=X,Xt=Xt,Y=Y,W=W,Z=Z,C.Y=C.Y,C.W=C.W,C.Z=C.Z)
    NC2 = mygmm(NC.U.func=NC.U0,data=data,loc.X=n.b+ncol(W)+1)
    NC.ate3 = as.numeric(NC2[grep("NC.ate2",names(NC2))])
    NC.se3 = as.numeric(NC2[grep("NC.se2",names(NC2))])
    NC.se3.HAC = as.numeric(NC2[grep("NC.se2.HAC",names(NC2))])
  }else{
    NC.model = try(
      gmm::gmm(g=Y~-1+W+Xt+X,
             x=~-1+Z+Xt+X,
             wmatrix="ident",
             vcov="iid"),
      silent=TRUE)
    NC.ate2 = summary(NC.model)$coef[grep("X",names(coef(NC.model))),"Estimate"]
    NC.se2  = summary(NC.model)$coef[grep("X",names(coef(NC.model))),"Std. Error"]
    NC.model = try(
      gmm::gmm(g=Y~-1+W+Xt+X,
               x=~-1+Z+Xt+X,
               wmatrix="ident",
               vcov="HAC"),
      silent=TRUE)
    NC.se2.HAC  = summary(NC.model)$coef[grep("X",names(coef(NC.model))),"Std. Error"]
    
    ### by hand GMM
    data = list(X=X,Xt=Xt,Y=Y,Z=Z,W=W)
    NC2 = mygmm(NC.U.func=NC.U0.nocov,data=data,loc.X=n.b+ncol(W))
    NC.ate3 = as.numeric(NC2[grep("NC.ate2",names(NC2))])
    NC.se3 = as.numeric(NC2[grep("NC.se2",names(NC2))])
    NC.se3.HAC = as.numeric(NC2[grep("NC.se2.HAC",names(NC2))])
    
    NC.ate = NC.se = NC.se.HAC = rep(NA,n.b)
  }#if(addcov==T){
  
  return(c(SC.ate=SC.ate, ## unconstrained OLS
           NC.ate=NC.ate,##NC ignore C
           NC.ate2=NC.ate2,##NC adjust C
           NC.ate3=NC.ate3, ##by hand
           SC.se =SC.se, ## unconstrained OLS
           NC.se =NC.se,NC.se.HAC=NC.se.HAC, ##NC ignore C
           NC.se2=NC.se2,NC.se2.HAC=NC.se2.HAC##NC adjust C
           ,NC.se3=NC.se3,NC.se3.HAC=NC.se3.HAC ##by hand
           ))
}#run.one

gen.lambda = function(dist.lambda){
  mylambda = log(1:t)+rnorm(t,sd=mysd)
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
    NC.ate2 = NC.weights[loc.X+0:(1-n.b)] #return slope followed by intercept
    NC.se2 = sqrt(NC.se2[loc.X+0:(1-n.b)])
    NC.se2.HAC = sqrt(NC.se2.HAC[loc.X+0:(1-n.b)])
  }else{
    NC.ate2 = NC.se2 = NC.se2.HAC = rep(NA,n.b)
  }
  return(c(NC.ate2=NC.ate2, NC.se2=NC.se2, NC.se2.HAC=NC.se2.HAC))
}

NC.U0.nocov = function(para,data){
  X=cbind(data$X)
  Xt=cbind(data$Xt)
  Y=cbind(data$Y);
  Z=cbind(data$Z)
  W=cbind(data$W)
  omega.u = para[1:ncol(W)]
  beta0 = para[ncol(W)+1]
  beta1 = para[ncol(W)+2]
  instrument = cbind(
    Z, #c(1-X)*Z, #if pre and post two-stage
    X,
    Xt
  )
  bridge = (
    X*beta0 + Xt*beta1 + W%*%omega.u
  )
  NC.obj=c(  Y-bridge   )*(instrument)
  return(NC.obj)
}

NC.U0 = function(para,data){
  X=cbind(data$X); Xt=cbind(data$Xt)
  Y=cbind(data$Y); C.Y=cbind(data$C.Y)
  Z=cbind(data$Z); C.Z=cbind(data$C.Z)
  W=cbind(data$W); C.W=cbind(data$C.W)
  theta.u = para[1]
  omega.u = para[1+1:ncol(W)]
  beta0 = para[1+ncol(W)+1]
  beta1 = para[1+ncol(W)+2]
  Ztheta = Z-C.Z*theta.u
  Comega = C.Y-C.W%*%omega.u
  instrument = cbind(
    Ztheta, #c(1-X)*Z, #if pre and post two-stage
    Comega,
    X,
    Xt
  )
  bridge = (
    X*beta0 + Xt*beta1 + W%*%omega.u + Comega*theta.u
  )
  NC.obj=c(  Y-bridge  )*(instrument)
  return(NC.obj)
}

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
