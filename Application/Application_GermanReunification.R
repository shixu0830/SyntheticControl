rm(list=ls())
placebo.test = F
library(dplyr)
library(foreign)
library(tidyr)
# Load Code 
my.filepath=""
source(paste0(my.filepath,"Application_GermanReunification_functions.R"))
# Load Data 
d <- read.dta(file=paste0(my.filepath,"repgermany.dta"))
################################################
#####define variables
################################################
unit.name="country"
units.name=unique(d[,unit.name])
time.name="year"
W.name=c("Austria", "Japan", "Netherlands", "Switzerland", "USA")
Y.name=c("West Germany")
Z.name= setdiff(units.name,c(W.name,Y.name))
outcome.name = "gdp"
predictor.name = c("infrate","trade","industry")
## did not include schooling or invest because too many values missing
################################################
#####impute missing data
################################################
d = d %>% group_by(country) %>% fill(industry,.direction="updown")
d = d %>% group_by(country) %>% fill(trade,.direction="down")
d = d %>% group_by(country) %>% fill(infrate,.direction="updown")
d = as.data.frame(d)
################################################
#####subset data
################################################
first.time = 1960
if(placebo.test==F){
  last.ctrl.time = 1990
  last.trt.time = 2003
}else{
  last.ctrl.time = 1975
  last.trt.time = 1990
  # subset d to pre-trt time period
  nrow(d); d = d[d[,time.name]%in%c(first.time:last.trt.time),]; nrow(d)
}
# select columns of d
d = d[,c(outcome.name,predictor.name,unit.name,time.name)]
################################################
####define time
################################################
(pretime = first.time:last.ctrl.time)
(postime = (last.ctrl.time+1):last.trt.time)
(myt0 = length(pretime))
(myt1 = length(postime))
########################
#### Get data       ####
########################
### treatment = X
X = c(rep(0,myt0),rep(1,myt1))
### outcome = Y, PIO = W
outcomes = getWYZ(dat=d,Y.name,W.name,Z.name,outcome.name,unit.name,time.name)
Y <<- as.matrix(outcomes$Y)
W <<- as.matrix(outcomes$W)
Z <<- as.matrix(outcomes$Z)
### covariates
C.Y=C.W=C.Z=NULL
covariates = getWYZ(dat=d,Y.name,W.name,Z.name,outcome.name=predictor.name,unit.name,time.name)
C.Y <<- as.matrix(covariates$Y)
C.W <<- as.matrix(covariates$W)
C.Z <<- as.matrix(covariates$Z)
########################
#### do estimation ####
########################
data.pre = list(Y=Y[1:myt0],W=W[1:myt0,],Z=Z[1:myt0,],C.Y=C.Y[1:myt0,],C.W=C.W[1:myt0,],C.Z=C.Z[1:myt0,])
data.all = list(X=X,Y=Y,W=W,Z=Z,C.Y=C.Y,C.W=C.W,C.Z=C.Z)

#####################
#### estimation  ####
#####################
### estimate C.coef and W.coef
### using pre-treatment data
### method 1
PI.pre_trt = optim(par = rep(0,ncol(W)+ncol(C.Y)),
                   fn = GMMF, mrf = PI.U.tobacco_application,
                   data=data.pre
)
### method 2
PI.pre_trt.gmm=gmm::gmm(
  g=PI.U.tobacco_application,
  x=data.pre,t0=rep(0,ncol(W)+ncol(C.Y)),
  vcov="iid"
  ,wmatrix="ident" ##this implements "one-step GMM w/ fixed W"
)
PI.pre_trt.gmm.HAC=gmm::gmm(
  g=PI.U.tobacco_application,
  x=data.pre,t0=rep(0,ncol(W)+ncol(C.Y)),
  vcov="HAC"
  ,wmatrix="ident"
  ,kernel="Truncated" ##Kernel-based HAC Covariance Matrix Estimation 
)
## same point estimate
rbind(PI.pre_trt$par,coef(PI.pre_trt.gmm),coef(PI.pre_trt.gmm.HAC))
summary(PI.pre_trt.gmm)$coef
## get predicted Y0 from PI
my.coef.pre = coef(PI.pre_trt.gmm)
Y.PI.Sync = predY0_each(para=my.coef.pre,data=data.all)
# (ATT.pre = mean(Y[1:myt0]-Y.PI.Sync[1:myt0]))
(ATT.post = mean(Y[myt0+1:myt1]-Y.PI.Sync[myt0+1:myt1]))
est = c(coef(PI.pre_trt.gmm),ATT.post)

####################
#### inference  ####
####################
PI.prepost_trt.gmm=gmm::evalGmm(
  g=PI.U.tobacco_application_addtrt,
  x=data.all,
  t0=est,
  tetw=est,
  vcov="iid"
)
PI.prepost_trt.gmm.HAC=gmm::evalGmm(
  g=PI.U.tobacco_application_addtrt,
  x=data.all,
  t0=est,
  tetw=est,
  vcov="HAC"
  ,kernel="Bartlett"
)

### point estimates and se
summary(PI.prepost_trt.gmm)$coef[9,]
#HC Estimate    Std. Error       t value      Pr(>|t|) 
# -1.200199e+03  4.186189e+02 -2.867044e+00  4.143251e-03 
# Estimate  Std. Error     t value    Pr(>|t|) 
# 125.2275342 211.3094004   0.5926264   0.5534312 
summary(PI.prepost_trt.gmm.HAC)$coef[9,]
#HAC Estimate    Std. Error       t value      Pr(>|t|) 
# -1.200199e+03  5.962796e+02 -2.012812e+00  4.413437e-02 
# Estimate  Std. Error     t value    Pr(>|t|) 
# 125.2275342 207.6419081   0.6030937   0.5464463 

### 95% confidence interval
confint(PI.prepost_trt.gmm)$test[9,]
#HC 0.025      0.975 
# -2020.6769  -379.7209
# 0.025     0.975 
# -288.9313  539.3863
confint(PI.prepost_trt.gmm.HAC)$test[9,]
#HAC 0.025     0.975 
# -2368.88538   -31.51242 
# 0.025     0.975 
# -281.7431  532.1982 

########################
#### Prep plot      ####
########################
## (1) get avg of all 38 states
outcomes.tmp = getWYZ(dat=d,Y.name=NULL,W.name=setdiff(units.name,Y.name),Z.name=NULL,outcome.name,unit.name,time.name)
All_controls=outcomes.tmp$W; dim(All_controls); rm(outcomes.tmp)
Y0.avg = cbind(apply(All_controls,1,mean))
## (2) sync (Y.Sync_Abadie) and regression (Y.regress) using Abadie SC method
if(placebo.test==F){
  load(paste0(my.filepath,"SCrslt.RData"))
  Y.Sync_Abadie=cbind(synthY0[row.names(synthY0)%in%c(pretime,postime),])
  load(paste0(my.filepath,"SCrslt.RData"))
  Y.regress=cbind(as.matrix(All_controls)%*%regress.W)
}else{
  load(paste0(my.filepath,"SCrslt_placebo.RData"))
  Y.Sync_Abadie=cbind(synthY0[row.names(synthY0)%in%c(pretime,postime),])
  load(paste0(my.filepath,"SCrslt.RData"))
  Y.regress=cbind(as.matrix(All_controls)%*%regress.W)
}

## change row names to time
time = unique(d[,time.name])
row.names(Y.PI.Sync) = time
row.names(Y0.avg) = time
row.names(Y.Sync_Abadie) = time
row.names(Y.regress) = time
myrows = time
########################
#### End Prep plot  ####
########################
mygrey = "grey90"
lty.orig = 1; col.orig = 1;
lty.controls = 3; col.controls = mygrey; 
lty.average = 4; col.average = "grey50"
lty.PI.pre = 1; col.PI.pre = 2;
lty.SC.regress = 3; col.SC.regress = 4;
lty.SC.abadie = 5; col.SC.abadie = 4;  
if(placebo.test==F){
  pdf(file=paste0(my.filepath,"SC_NewApplication.pdf"),w=6,h=5)
  myrange=c(0,31000)#range(c(0,Y,Y.PI.Sync))
}else{
  pdf(file=paste0(my.filepath,"SC_NewApplication_placebo.pdf"),w=6,h=5)
  myrange=c(0,21000)#range(c(0,Y,Y.PI.Sync))
}

par(mfrow=c(1,1),mar=c(0,0,0,0),mai=c(0.7,0.7,0.1,0.1),omi=c(0,0,0,0))
plot(x=time,Y,ylim=myrange,type="l",ylab="",xlab="",col=col.orig,lty=lty.orig,lwd=2)
axis(side=2,line=1,tick=F,at=mean(myrange),labels="Per Capita GDP")
axis(side=1,line=1,tick=F,at=mean(time),labels="Year")
# for(i in 1:ncol(All_controls)){
#   lines(x=time,y=All_controls[,i],lty=lty.controls,col=col.controls)
# }
# lines(x=time,y=Y,col=col.orig,lty=lty.orig)
abline(v=last.ctrl.time,lty=3)
# axis(side=1,line=-2,tick=F,at=last.ctrl.time,labels=last.ctrl.time,cex.axis=0.8,font=4)
if(last.ctrl.time==1990){
  axis(side=1,line=-2,tick=F,at=1975,labels="Pre-treatment",cex.axis=0.8,font=4)
  axis(side=1,line=-2,tick=F,at=1999,labels="Post-treatment",cex.axis=0.8,font=4)
  arrows(x0=1993, y0=2500, x1 = 1990, y1=2500, length=0.1)
  text(x=1997,y=2500,labels="Reunification",cex=0.8)
}else{
  axis(side=1,line=-2,tick=F,at=1967,labels="Pre-treatment",cex.axis=0.8,font=4)
  axis(side=1,line=-2,tick=F,at=1985,labels="Post-treatment",cex.axis=0.8,font=4)
  arrows(x0=1977.2, y0=2500, x1 = 1975, y1=2500, length=0.1)
  text(x=1981.8,y=2500,labels="Placebo Reunification",cex=0.8)
}
lines(x=time,y=Y0.avg,lty=lty.average,col=col.average,lwd=2)
lines(x=time,y=Y.PI.Sync,lty=lty.PI.pre,col=col.PI.pre,lwd=2)
lines(x=time,y=Y.regress,col=col.SC.regress,lty=lty.SC.regress,lwd=2)
lines(x=time,y=Y.Sync_Abadie,lty=lty.SC.abadie,col=col.SC.abadie,lwd=2)
legend("topleft",
       legend = c("West Germany (treated)",
                  # "Control countries (N=16)",
                  "Average of 16 control countries",
                  "Syn West Germany (PI)",
                  "Syn West Germany (OLS)",
                  "Syn West Germany (Abadie et al. 2015)"),
       lty=c(lty.orig,#lty.controls,
             lty.average,lty.PI.pre,lty.SC.regress,lty.SC.abadie),
       col=c(col.orig,#col.controls,
             col.average,col.PI.pre,col.SC.regress,col.SC.abadie),
       bty="n",lwd=2,cex=1)
dev.off()


#### plot differences
myrange=range(c(Y-Y.PI.Sync,Y-Y.Sync_Abadie,Y-Y.regress,Y-Y0.avg))
plot(x=time,y=rep(0,length(time)),ylim=myrange,type="l",ylab="",xlab="",col=col.orig,lty=lty.orig,lwd=2)
axis(side=2,line=1,tick=F,at=mean(myrange),labels="Per Capita GDP")
axis(side=1,line=1,tick=F,at=mean(time),labels="Year")
# for(i in 1:ncol(All_controls)){
#   lines(x=time,y=All_controls[,i],lty=lty.controls,col=col.controls)
# }
# lines(x=time,y=Y,col=col.orig,lty=lty.orig)
abline(v=last.ctrl.time,lty=3)
# axis(side=1,line=-2,tick=F,at=last.ctrl.time,labels=last.ctrl.time,cex.axis=0.8,font=4)
if(placebo.test==F){
  axis(side=1,line=-2,tick=F,at=1975,labels="Pre-treatment",cex.axis=0.8,font=4)
  axis(side=1,line=-2,tick=F,at=1999,labels="Post-treatment",cex.axis=0.8,font=4)
  arrows(x0=1993, y0=-3500, x1 = 1990, y1=-3500, length=0.1)
  text(x=1997,y=-3500,labels="Reunification",cex=0.8)
}else{
  axis(side=1,line=-2,tick=F,at=1967,labels="Pre-treatment",cex.axis=0.8,font=4)
  axis(side=1,line=-2,tick=F,at=1985,labels="Post-treatment",cex.axis=0.8,font=4)
  arrows(x0=1977.2, y0=1000, x1 = 1975, y1=1000, length=0.1)
  text(x=1981.8,y=1000,labels="Placebo Reunification",cex=0.8)
}
lines(x=time,y=Y-Y0.avg,lty=lty.average,col=col.average,lwd=2)
lines(x=time,y=Y-Y.PI.Sync,lty=lty.PI.pre,col=col.PI.pre,lwd=2)
lines(x=time,y=Y-Y.regress,col=col.SC.regress,lty=lty.SC.regress,lwd=2)
lines(x=time,y=Y-Y.Sync_Abadie,lty=lty.SC.abadie,col=col.SC.abadie,lwd=2)
legend("topleft",
       legend = c("West Germany (treated)",
                  # "Control countries (N=16)",
                  "Average of 16 control countries",
                  "Syn West Germany (PI)",
                  "Syn West Germany (OLS)",
                  "Syn West Germany (Abadie et al. 2015)"),
       lty=c(lty.orig,#lty.controls,
             lty.average,lty.PI.pre,lty.SC.regress,lty.SC.abadie),
       col=c(col.orig,#col.controls,
             col.average,col.PI.pre,col.SC.regress,col.SC.abadie),
       bty="n",lwd=2,cex=1)




