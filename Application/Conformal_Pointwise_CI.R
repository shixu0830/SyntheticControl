rm(list = ls())

library(dplyr)
library(foreign)
library(tidyr)

## obtain conformal pointwise confidence intervals
ConformalPointwiseCI <- 
  function(data, t0, t1, grid_list) {
    X <- cbind(data$X)
    Y <- cbind(data$Y); C.Y <<- cbind(data$C.Y)
    Z <- cbind(data$Z); C.Z <- cbind(data$C.Z)
    W <<- cbind(data$W); C.W <- cbind(data$C.W)
    
    X1 <- c(rep(0, t0), 1)
    
    lb0.1 <- lb0.2 <- ub0.1 <- ub0.2 <- rep(NA, t1)
    
    
    for (tt in 1:t1) {
      C1.Y <- cbind(data$C.Y)[c(1:t0, t0 + tt),]
      Z1 <- Z[c(1:t0, t0 + tt),]
      C1.Z <- cbind(data$C.Z)[c(1:t0, t0 + tt),]
      W1 <- cbind(data$W)[c(1:t0, t0 + tt),]
      C1.W <- cbind(data$C.W)[c(1:t0, t0 + tt),]
      
      grid <- grid_list[[tt]]
      SC.pval <- NA * grid
      for (bb in seq_along(grid)) {
        Y1 <- cbind(c(data$Y[1:t0], data$Y[t0 + tt] - grid[bb]))
        dat.h0 <- list(X = X1, Y = Y1, W = W1, Z = Z1,
                       C.Y = C1.Y, C.W = C1.W,C.Z = C1.Z)
        
        PI.h0 <- optim(par = rep(0, ncol(W) + ncol(C.Y)),
                       fn = GMMF, mrf = PI.U.tobacco_application,
                       data = dat.h0)
        
        ## get predicted Y0 from PI
        est <- PI.h0$par
        Y.PI.SC <- predY0_each(para = est, data = dat.h0)
        # (ATT.pre = mean(Y[1:myt0]-Y.PI.Sync[1:myt0]))
        SC.resid <- Y1 - Y.PI.SC
        SC.stat <- abs(SC.resid[t0 + 1])
        SC.pval[bb] <- 1 - mean(SC.stat > abs(SC.resid))
      }
      # find upper bound
      lb0.1.id <- max(which(SC.pval[1:which.max(SC.pval)] < 0.1), na.rm = T) + 1
      lb0.2.id <- max(which(SC.pval[1:which.max(SC.pval)] < 0.2), na.rm = T) + 1
      
      ub0.1.id <- which.max(SC.pval) + min(which(SC.pval[which.max(SC.pval):length(SC.pval)] < 0.1), na.rm = T) - 2
      ub0.2.id <- which.max(SC.pval) + min(which(SC.pval[which.max(SC.pval):length(SC.pval)] < 0.2), na.rm = T) - 2
      
      lb0.1[tt] <- grid[lb0.1.id]; lb0.2[tt] <- grid[lb0.2.id]
      ub0.1[tt] <- grid[ub0.1.id]; ub0.2[tt] <- grid[ub0.2.id]
    } 
    return(list(lb0.1 = round(lb0.1), lb0.2 = round(lb0.2),
                ub0.1 = round(ub0.1), ub0.2 = round(ub0.2)))
  }


placebo.test = T

# Load Code 
my.filepath = ""
source(paste0(my.filepath,"Application_GermanReunification_functions.R"))
# Load Data 
d <- readRDS(file=paste0(my.filepath,"repgermany.rds"))
################################################
#####define variables
################################################
unit.name="country"
units.name=unique(d[,unit.name])
time.name="year"
W.name <<- c("Austria", "Japan", "Netherlands", "Switzerland", "USA")
Y.name <<- c("West Germany")
Z.name <<- setdiff(units.name,c(W.name,Y.name))
outcome.name <<- "gdp"
predictor.name <<- c("infrate","trade","industry")
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
covariates = getWYZ(dat=d,Y.name,W.name,Z.name,outcome.name=predictor.name,unit.name,time.name)
C.Y <<- as.matrix(covariates$Y)
C.W <<- as.matrix(covariates$W)
C.Z <<- as.matrix(covariates$Z)
########################
#### do estimation ####
########################
data.pre = list(Y=Y[1:myt0],W=W[1:myt0,],Z=Z[1:myt0,],C.Y=C.Y[1:myt0,],C.W=C.W[1:myt0,],C.Z=C.Z[1:myt0,])
data.all = list(X=X,Y=Y,W=W,Z=Z,C.Y=C.Y,C.W=C.W,C.Z=C.Z)

## crude search
crude_grid <- replicate(myt1, seq(-3500, 1000, 10), simplify = F)
CPCI_crude <- ConformalPointwiseCI(data = data.all, t0 = myt0, t1 = myt1,
                                   grid_list = crude_grid)
newgrid <- lapply(1:myt1, function(tt) {
  (CPCI_crude$lb0.1[tt] - 9):(CPCI_crude$ub0.1[tt] + 9)
})

CPCI <- ConformalPointwiseCI(data = data.all, t0 = myt0, t1 = myt1,
                             grid_list = newgrid)
if (placebo.test) {
  save(CPCI, file = "conformal_ci_rslt_placebo.rdata")
} else {
  save(CPCI, file = "conformal_ci_rslt.rdata")
}


## plot the pointwise confidence intervals of the treatment effect
