predY0_each = function(para,data=data.all){
  ### get data and estimates
  X=cbind(data$X)
  Y=cbind(data$Y);C.Y=cbind(data$C.Y)
  Z=cbind(data$Z);C.Z=cbind(data$C.Z)
  W=cbind(data$W);C.W=cbind(data$C.W)
  theta.est=para[1:ncol(C.Y)]
  omega.est=para[ncol(C.Y)+1:ncol(W)]
  ### predict Y0
  Womega = NULL
  for(i in predictor.name){
    Womega = cbind(
      Womega,
      C.W[,grep(i,colnames(C.W))]%*%omega.est
    )
  }
  # equivalent to below
  # Womega=NULL
  # for(t in 1:nrow(C.W)){
  #   ABC = matrix(C.W[t,],ncol=ncol(C.Y),byrow=T)
  #   Womega = rbind(
  #     Womega, t(omega.est)%*%ABC
  #   )
  # }
  Comega = C.Y-Womega
  if(ncol(C.Y)==1){
    Comega=cbind(apply(Comega,1,mean))
  }
  Y0 = as.matrix(W)%*%c(omega.est) + Comega%*%c(theta.est)
  return(Y0)
}

getWYZ = function(dat=d, Y.name, W.name,Z.name,
                  outcome.name="gdp", unit.name="country", time.name="year"){
  Synth.Y = reshape(dat[,c(unit.name,outcome.name,time.name)],timevar=unit.name,direction="wide",idvar=time.name)
  Synth.Y = Synth.Y[,-grep(time.name,names(Synth.Y))]
  ind.Y = which(names(Synth.Y)==c(sapply(outcome.name,FUN=function(x){paste0(x,".",Y.name)})))
  ind.W = which(names(Synth.Y)%in%c(sapply(outcome.name,FUN=function(x){paste0(x,".",W.name)})))
  ind.Z = which(names(Synth.Y)%in%c(sapply(outcome.name,FUN=function(x){paste0(x,".",Z.name)})))
  Y = cbind(Synth.Y[,ind.Y])
  W = cbind(Synth.Y[,ind.W])
  Z = cbind(Synth.Y[,ind.Z])
  return(list(W=W,Y=Y,Z=Z))
}


PI.U.tobacco_application = function(para,data){
  X=cbind(data$X)
  Y=cbind(data$Y);C.Y=cbind(data$C.Y)
  Z=cbind(data$Z);C.Z=cbind(data$C.Z)
  W=cbind(data$W);C.W=cbind(data$C.W)
  # print(dim(C.Z))
  theta = para[1:ncol(C.Y)]
  omega.u = para[ncol(C.Y)+1:ncol(W)]
  Z.state.num = as.character(Z.name)
  W.state.num = as.character(W.name)
  C.Z.num = unlist(lapply(strsplit(x=colnames(C.Z),split="[.]"),FUN=function(x){x[2]}))
  C.W.num = unlist(lapply(strsplit(x=colnames(C.W),split="[.]"),FUN=function(x){x[2]}))
  Ztheta = NULL
  for(i in 1:ncol(Z)){
    C.Z.i = C.Z[,C.Z.num==Z.state.num[i]]
    if(ncol(C.Y)==1){
      C.Z.i=cbind(apply(C.Z.i,1,mean))
    }
    Ztheta = cbind(
      Ztheta,
      Z[,i]-C.Z.i%*%theta
    )
  }
  Womega = NULL
  for(i in predictor.name){
    Womega = cbind(
      Womega,
      C.W[,grep(i,colnames(C.W))]%*%omega.u
    )
  }
  # equivalent to below
  # Womega=NULL
  # for(t in 1:nrow(C.W)){
  #   ABC = matrix(C.W[t,],ncol=ncol(C.Y),byrow=T)
  #   Womega = rbind(
  #     Womega, t(omega.est)%*%ABC
  #   )
  # }
  Comega = C.Y-Womega
  if(ncol(C.Y)==1){
    Comega=cbind(apply(Comega,1,mean))
  }
  instrument = cbind(
    Ztheta,
    Comega
  )
  PI.obj=c(  Y - W%*%omega.u - Comega%*%theta  )*(instrument)
  return(PI.obj)
}
PI.U.tobacco_application_addtrt=function(para,data=data.all){
  X=cbind(data$X)
  Y=cbind(data$Y);C.Y=cbind(data$C.Y)
  Z=cbind(data$Z);C.Z=cbind(data$C.Z)
  W=cbind(data$W);C.W=cbind(data$C.W)
  theta = para[1:ncol(C.Y)]
  omega.u = para[ncol(C.Y)+1:ncol(W)]
  beta=para[ncol(C.Y)+ncol(W)+1]
  Z.state.num = as.character(Z.name)
  W.state.num = as.character(W.name)
  C.Z.num = unlist(lapply(strsplit(x=colnames(C.Z),split="[.]"),FUN=function(x){x[2]}))
  C.W.num = unlist(lapply(strsplit(x=colnames(C.W),split="[.]"),FUN=function(x){x[2]}))
  Ztheta = NULL
  for(i in 1:ncol(Z)){
    C.Z.i = C.Z[,C.Z.num==Z.state.num[i]]
    if(ncol(C.Y)==1){
      C.Z.i=cbind(apply(C.Z.i,1,mean))
    }
    Ztheta = cbind(
      Ztheta,
      Z[,i]-C.Z.i%*%theta
    )
  }
  Womega = NULL
  for(i in predictor.name){
    Womega = cbind(
      Womega,
      C.W[,grep(i,colnames(C.W))]%*%omega.u
    )
  }
  Comega = C.Y-Womega
  if(ncol(C.Y)==1){
    Comega=cbind(apply(Comega,1,mean))
  }
  Ztheta = cbind(Ztheta)
  for(i in 1:ncol(Ztheta)){
    Ztheta[,i] = Ztheta[,i]*(1-X)
  }
  Comega = cbind(Comega)
  for(i in 1:ncol(Comega)){
    Comega[,i] = Comega[,i]*(1-X)
  }
  
  instrument = cbind(
    Ztheta,
    Comega,
    X
  )
  PI.obj=c(  Y - X*beta -W%*%omega.u - Comega%*%theta  )*(instrument)
  return(PI.obj)
}

GMMF <- function(mrf,para,data){
  g0 <- mrf(para=para,data=data)
  g <- apply(g0,2,mean)
  gmmf <- sum(g^2)
  
  return(gmmf)
}
G <- function(bfun,para,data){
  G <- numDeriv::jacobian(func=G1,bfun=bfun,x=para,data=data)
  return(G)
}
G1 <- function(bfun,para,data){
  G1 <- apply(bfun(para,data),2,mean,na.rm=T)
  return(G1)
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
  Sigma <- bG%*%Omega%*%t(bG)/spsz
  hacSigma <- bG%*%hacOmega%*%t(bG)/spsz
  
  return(list(var=Sigma, hacvar=hacSigma))
}
HAC_VAREST_prepost <- function(bfun,para,q,data){
  bG <- solve(G(bfun,para,data))
  bg <- bfun(para,data)
  spsz <- nrow(bg)
  bg[is.na(bg)]=0
  n.weight=diag(c(rep(1/t0,ncol(bg)-1), 1/t))
  hacOmega <- Omega <- t(bg)%*%bg%*%n.weight
  for(i in 1:q){
    Omega_i <- t(bg[-(1:i),])%*%bg[1:(spsz-i),]%*%n.weight
    hacOmega <- hacOmega + (1 - i/(q+1))*(Omega_i + t(Omega_i))
  }
  Sigma <- bG%*%Omega%*%t(bG)%*%n.weight
  hacSigma <- bG%*%hacOmega%*%t(bG)%*%n.weight
  # print((bg))
  # spsz <- nrow(bg) #2326.141
  # bg[is.na(bg)]=0
  # n.weight=diag(c(rep(1/t0,ncol(bg)-1), 1/t))
  # hacOmega <- Omega <- t(bg[,9])%*%bg[,9]%*%n.weight[9,9]
  # for(i in 1:q){
  #   Omega_i <- t(bg[-(1:i),][,9])%*%bg[1:(spsz-i),][,9]%*%n.weight[9,9]
  #   hacOmega <- hacOmega + (1 - i/(q+1))*(Omega_i + t(Omega_i))
  # }
  return(list(var=Sigma, hacvar=hacSigma))
}

