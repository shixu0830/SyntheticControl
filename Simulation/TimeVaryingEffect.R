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
n.rep=10000
t0.all=c(200,500,1000)
n.units.all=c(1+2,1+10,1+20)
arg1 = as.numeric(Sys.getenv("arg1")) ## arg1=1,2 dist.epsilon
arg2 = as.numeric(Sys.getenv("arg2")) ## arg2=0,1 addcov
arg3 = as.numeric(Sys.getenv("arg3")) ## arg3=1-10 myseeds
myfilepath = ""
source(paste0(myfilepath,"myfunctions_TimeVaryingEffect.R"))
dist.epsilon = c("iid","AR")[arg1]
addcov=(arg2==1)
myseeds = 1:n.rep+(arg3-1)*n.rep
theta=1
true.beta=c(1,1);n.b=length(true.beta)
mysd = 1
epsilonARMAcor=0.1

rslt.all=NULL
for(n.units in n.units.all){
  for(t0.i in 1:length(t0.all)){
    t0=t0.all[t0.i]
    addt=t0
    t=addt+t0
    ## pre/post treatment indicator
    X=c(rep(0,t0),rep(1,t-t0))
    Xt = X*(1:t)/t ##"(1:t)/t" is t/T
    ## save results
    ate.all = NULL
    for(rep in myseeds){
      rslt=run.one(rep,n.units,t,t0)
      ate.all=rbind(ate.all,rslt)
    }#for(rep in 1:n.rep){
    rslt.all = append(rslt.all,
                      list(ate.all=ate.all))
  }#for(t0.i in 1:3){
}#for(n.units in c(1+2,1+10,1+20))

save(rslt.all,file=paste0(myfilepath,"TV",
                          "epsilon_",dist.epsilon,
                          "cov_",addcov,
                          "t0_",paste(t0.all,collapse="_"),
                          "nunits_",paste(n.units.all,collapse="_"),
                          "nrep",n.rep,
                          "argseeds",arg3,
                          ".RData"))