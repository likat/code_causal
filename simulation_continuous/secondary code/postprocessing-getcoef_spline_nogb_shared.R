library(nlme)
library(dplyr)
library(survey)
library(ggplot2)
library(sampling)
library(cmdstanr)
library("rootSolve")
library(LaplacesDemon)

type_y <- "cts"
type_outcome = "nogb"
type_overlap = "shared"
cluster = T
if(cluster){
  pathtocode = "./"
}else{
  pathtocode = "../../"
}
source(paste0(pathtocode,"CODE/helpers.R"))
source(paste0(pathtocode, "CODE/calculate_sate.R"))
source(paste0(pathtocode, "CODE/calculate_sate_x.R"))
source(paste0(pathtocode, "CODE_cts/calculate_pencomp_cts.R"))
# source(paste0("../causal-with-x2-lower-z/CODE_cts/calculate_pencomp_cts.R"))

source(paste0(pathtocode, "CODE_cts/calculate_wfpbb_pencomp_cts.R"))
# source(paste0("../causal-with-x2-lower-z/CODE_cts/calculate_wfpbb_pencomp_cts.R"))

# source(paste0(pathtocode, "CODE_cts/imputeF_stan_cts.R"))
source(paste0(pathtocode, "CODE_cts/imputeF_cts.R"))

# source(paste0(pathtocode, "CODE_cts/run_AIPTW_simulation_cts.R"))
# source(paste0(pathtocode, "CODE/run_AIPTW_new.R"))

source(paste0(pathtocode, "CODE/gcomputeFunc_sim.R"))
source(paste0(pathtocode, "CODE/calculate_pencomp_nospline.R"))
# source(paste0(pathtocode, "CODE/calculate_wtd.R"))
source(paste0(pathtocode, "CODE/calculate_wtd_2stage.R"))
source(paste0(pathtocode, "STEP0-poptingting.R"))

npop <- nrow(pop)
outcome_mod <- NA
propensity_mod <- NA
x_varnames = NA
y_varname = "Ysamp"
catecat_varname = "catecat"
# x_name <- "race4"
x_name <- c("X1")
# x_name_mod <- "factor(race4)"
x_name_mod <- c("X1")
z_name <- "Z"

### SIMULATION PARS
sim = 200
numboot=200
set.seed(42)
# numclus <- 400
# numclussamp <- 26
nknot = 15
# nstrat = 2
numclussamp_perstrat <- numclussamp / nstrat

## OUTCOME GENERATION
#***************************************************************
a0 <- 0
ax1 <- 1
ax2 <- 1
ax1x2 <- 0.5
az <- -0.75
b0 <- 0
bx1 <- 1 # X1
bx2 <- 1 # X2
bA <- 5
bAx1 <- 2
bAx2 <- 0
bz <- 1.5
bAz <- 2
bax1z <- 0
bax2z <- 0.8
#***************************************************************
if(type_outcome == "gb"){
  pop <- pop %>% mutate(
    mu1 = cluseff + bA + bx1*X1 + bAx1*X1 + bx2*X2 + bz*Z+ bAz*Z,
    mu0 = cluseff + bx1*X1 + bx2*X2 + bz*Z
  )
  print(type_outcome)
}else if(type_outcome == "nogb"){
  #  pop <- pop %>% mutate(
  #    mu1 = bA + bx1*X1  + bz*Z+ bAx1*X1,
  #    mu0 = bx1*X1 + bz*Z,
  #  )
  pop <- pop %>% mutate(
    mu1 = cluseff +bA + bx1*X1 + bAx1*X1 + bx2*X2 + bz*Z,
    mu0 = cluseff +bx1*X1 + bx2*X2 + bz*Z
  )
  print(type_outcome)
  
}
pop$Y1 <- rnorm(n=npop, pop$mu1)
pop$Y0 <- rnorm(n=npop, pop$mu0)

cluster_stratum_ref_pop <- pop %>% group_by(stratum, cluster) %>% summarize(cts=n()) 
stratum_ref_pop <- cluster_stratum_ref_pop %>% group_by(stratum) %>% summarize(cts=n()) 
pop$fpc1 <- numclus / nstrat
pop$fpc2 <- npop/numclus
## PROPENSITY GENERATION
if(type_overlap == "shared"){
  # propensity_mod <- as.formula("A ~ X1 + Z")
  # pop$ptrt <-  expit(ax1*pop$X1 + az*pop$Z)
  propensity_mod <- as.formula("A ~ X1 + X2 + Z + X1:X2")
  pop$ptrt <-  expit(ax1*pop$X1 + ax2*pop$X2 + ax1x2*pop$X1*pop$X2 +
                       az*pop$Z)
  
}else if(type_overlap == "notshared"){
  
  #  propensity_mod <- as.formula("A ~ X1 ")
  #  pop$ptrt <-  expit(ax1*pop$X1 )
  propensity_mod <- as.formula("A ~ X1 + X2 + X1:X2")
  pop$ptrt <-  expit(ax1*pop$X1 + ax2*pop$X2 + ax1x2*pop$X1*pop$X2 )
  
}

x_varnames <- c(x_name_mod, z_name, paste0(x_name_mod,":",z_name))

varnames_outcome <-  c("X2", "Z")
varnames_outcome_generalize <- c(x_name_mod, z_name, paste0(x_name_mod, ":", z_name))
outcome_mod <- as.formula("Ysamp ~ X2 + Z")
outcome_mod_nospline <- as.formula("Ysamp ~ X1 + X2 + Z")
# outcome_mod <- "Ysamp ~ X1 + X2 + Z"
# outcome_mod_nospline <- formulaF(varList=varnames_outcome_nospline, y.name=y_varname)
outcome_mod_generalize <- formulaF(varList=varnames_outcome_generalize, y.name=y_varname)

## PARAMETERS FOR AIPTW
covariateXnames <- c("X1","X2", "Z")
propensity_mod_aiptw <- propensity_mod
outcome_mod_aiptw <- outcome_mod_nospline
ncarlo_aiptw <- 2000


# true PATE
pate_true <- mean(pop$Y1 - pop$Y0)

### ESTIMATE CONTAINERS
sate <- rep(NA, sim)
resnames <- c("est","withinvar","btwnvar","mivar","ci_tlength","ci_emplength","cicov_t","cicov_emp") # must match names from rubin function
resnames_aiptw <- c("est", "var", "cilower","ciupper")
time_aiptw <- 0
time_pencompwfpbb <- 0
wtdest_pencomp <- vector(mode="list", length=8)
wtdest_pencompwfpbb <- vector(mode="list", length=8)
wtdest_nospline <- vector(mode="list", length=8)

aiptw <- vector(mode = "list", length = 4)
synthpop_pencompwfpbb <- vector(mode="list", length = 8)

twostage_pencomp <- vector(mode="list", length=8)
twostage_pencompwfpbb <- vector(mode="list", length=8)
twostage_nospline <- vector(mode="list", length=8)

ate_pencomp <- vector(mode="list", length=8)
ate_pencomp_wfpbb <- vector(mode="list", length=8)
ate_nospline <- vector(mode="list", length=8)

coef_pencomp <- vector(mode = "list", length = 4)
coef_pencompwfpbb <- vector(mode="list", length = 4)

names(wtdest_pencomp) <- resnames
names(wtdest_pencompwfpbb) <- resnames
names(wtdest_nospline) <- resnames
names(aiptw) <- resnames_aiptw
names(synthpop_pencompwfpbb) <- resnames
names(twostage_pencomp) <- resnames
names(twostage_pencompwfpbb) <- resnames
names(twostage_nospline) <- resnames
names(ate_pencomp) <- resnames
names(ate_pencomp_wfpbb) <- resnames
names(ate_nospline) <- resnames

bootflag <- 0
bootflagct <- 0
## To replicate a true experiment, recalculate the outcome every sample and treatment assignment
pop$popID <- 1:nrow(pop)

## To replicate a true experiment, recalculate the outcome every sample and treatment assignment
pop$popID <- 1:nrow(pop)

## To replicate a true experiment, recalculate the outcome every sample and treatment assignment
## HELPER FUNCTION: Make line segments for the splines #####

splinesegfun2 <- function(newdata,knotloc, fixedeff, randeff,sampdf,type_est){
  linearB=NULL
  xlocvec <- seq(min(newdata$pstar), max(newdata$pstar), by=0.01)
  linearB = outer(xlocvec, knotloc, "-")
  # linearB =outer(c(min(pstar),knotloc, max(pstar)), knotloc, "-")
  linearB =linearB * (linearB > 0)
  
  # designM = matrix(1,ncol=ncol(linearB),nrow=1)
  designM = cbind(1, xlocvec)
  designM=as.matrix(designM)
  fittedcoef <- fixedeff
  randomcoef <-randeff
  # predicted = designM %*% model$coefficients$fixed + as.matrix(linearB) %*% as.vector(unlist(model$random)) + rnorm(nrow(newdata), 0, model$sigma)
  predicted = designM %*%( fittedcoef %>% as.numeric())+ as.matrix(linearB) %*%as.numeric(randomcoef)
  plotdf <- data.frame(
    predicted = predicted,
    xlocs = xlocvec,
    type = type_est
  )
  res <- list(predicted = predicted,
              plotdf = plotdf)
  
  return(res)
}

logit <- function(x){
  return(log(x/(1-x)))
}

# draw the sample
set.seed(1)
pop$I <- 0
pop$Iclus <- 0
for(k in 1:nstrat){
  currentcluspool <- sort(unique(pop$cluster[pop$stratum ==k]))
  currentsampledclus <- sample(currentcluspool, size = numclussamp_perstrat)
  pop$Iclus[pop$cluster %in% currentsampledclus] <- 1
}
pop$I[pop$Iclus == 1] <- rbinom(n=sum(pop$Iclus == 1), size=1, prob = pop$pselectf2[pop$Iclus == 1])
# pop$I <- rbinom(n=npop, size=1, prob = nsamp/nrow(pop))

samp <- pop[pop$I==1,]
samp$ind <- 1:nrow(samp)
# samp$wts <- npop*normalize(1/samp$pselect)
samp$wts <- npop*normalize(1/samp$pselectf1f2)

nsamp_sim <- nrow(samp)
samp$Ysamp <- NA
# ### assign the treatment
samp$A <- NULL
samp$A <- rbinom(n=nsamp_sim, size=1, prob = samp$ptrt)
samp$Ysamp[samp$A==1] <- samp$Y1[samp$A==1]
samp$Ysamp[samp$A==0] <- samp$Y0[samp$A==0]
size_cnt <- sum(samp$A==0)
size_trt <- sum(samp$A==1)
ind_cnt <- which(samp$A==0)
ind_trt <- which(samp$A==1)
trtgrp0_popind <- samp$popID[samp$A==0]
trtgrp1_popind <- samp$popID[samp$A==1]

cluster_stratum_ref <- samp %>% group_by(stratum, cluster) %>% summarize(cts=n()) 
cluster_stratum_ref <- cluster_stratum_ref %>% group_by(stratum) %>% mutate(nclus = length(unique(cluster)),
                                                                            mstar = nclus / (nclus-1)) %>% select(-cts)

## Read in results
knots_pencompppsboot <- vector(mode = 'list', length = 2)
coef_pencompppsboot <- vector(mode = "list", length = 1)
knots_wfpbbpenppsboot <- vector(mode = 'list', length = 2)
coef_wfpbbpenppsboot <- vector(mode = "list", length = 1)    

for(boot in 1:200){ ## START bootstrap loop
  resboot <- readRDS(paste0("results-raw-2stage-coefonly/", type_outcome, type_overlap,"/results_tingting_noaiptw", type_outcome, type_overlap,boot, ".RDS"))

  #----- write to bootstrap containers-----
  knots_pencompppsboot<- Map(rbind, knots_pencompppsboot, resboot$knotloc_pencomp)
  knots_wfpbbpenppsboot<- Map(rbind, knots_wfpbbpenppsboot, resboot$knotloc_pencompwfpbb)
  coef_pencompppsboot<- Map(rbind, coef_pencompppsboot, resboot$coef_pencomp)
  coef_wfpbbpenppsboot<- Map(rbind, coef_wfpbbpenppsboot, resboot$coef_pencompwfpbb)

  print(boot)
}
knots_wfpbbpenpps <- knots_wfpbbpenppsboot %>% lapply(colMeans)
coef_wfpbbpenpps <- coef_wfpbbpenppsboot %>% lapply(colMeans, na.rm=T)
names(coef_wfpbbpenpps) <- c("fixed0", "fixed1", "rand0", "rand1")
knots_pencomppps <- knots_pencompppsboot %>% lapply(colMeans)
coef_pencomppps <- coef_pencompppsboot %>% lapply(colMeans)
names(coef_pencomppps) <-  c("fixed0", "fixed1", "rand0", "rand1")
  
  

## PLOTS #####
samp$pstar <- log(samp$ptrt/(1-samp$ptrt))
samp1_pps <- samp %>% filter(A==1)
samp0_pps <- samp %>% filter(A==0)

#--- For PPS plot ------
# use the pstar from pps samp instead
splines1pps_pstarpps <- splinesegfun2(newdata = samp1_pps,
                             knotloc = knots_pencomppps[[2]], 
                             fixedeff = coef_pencomppps$fixed1, 
                             randeff = coef_pencomppps$rand1,
                             type_est = "PENCOMP pps")
splines0pps_pstarpps <- splinesegfun2(newdata = samp0_pps,
                             knotloc = knots_pencomppps[[1]], 
                             fixedeff = coef_pencomppps$fixed0, 
                             randeff = coef_pencomppps$rand0,
                             type_est = "PENCOMP pps")
splines1wfpps_pstarpps <- splinesegfun2(newdata = samp1_pps,
                               knotloc = knots_wfpbbpenpps[[2]], 
                               fixedeff = coef_wfpbbpenpps$fixed1, 
                               randeff = coef_wfpbbpenpps$rand1,
                               type_est = "PENCOMP-WFPBB pps")
splines0wfpps_pstarpps <- splinesegfun2(newdata = samp0_pps,
                               knotloc = knots_wfpbbpenpps[[1]], 
                               fixedeff = coef_wfpbbpenpps$fixed0, 
                               randeff = coef_wfpbbpenpps$rand0,
                               type_est = "PENCOMP-WFPBB pps")

spline1df_pps <- bind_rows(splines1pps_pstarpps$plotdf, 
                           splines1wfpps_pstarpps$plotdf)
spline0df_pps <- bind_rows(splines0pps_pstarpps$plotdf, 
                           splines0wfpps_pstarpps$plotdf)

## GGPLOTS HERE ****************************
## SRS plot with all 4 ---------------------
library(latex2exp)
## ******************************************
## PPS plot with 2 methods ---------------------
spline1df_pps <- spline1df_pps %>% mutate(
  type = case_when(
    type == "PENCOMP pps" ~ "PENCOMP",
    type == "PENCOMP-WFPBB pps" ~ "WFPBB-PENCOMP"
  )
)
spline0df_pps <- spline0df_pps %>% mutate(
  type = case_when(
    type == "PENCOMP pps" ~ "PENCOMP",
    type == "PENCOMP-WFPBB pps" ~ "WFPBB-PENCOMP"
  )
)
pdf(file = paste0("results-clean/coef-trt1-",type_outcome, type_overlap, ".pdf"),width = 8*0.8, height = 6*0.8)
ggplot(data=samp1_pps,aes(y=Ysamp, x=pstar))+
  geom_point(alpha=0.2,color="gray", aes(size = wts)) + 
  # geom_point(data = spline1df, aes(y = predicted, x = xlocs, group = type, color = type)) + 
  geom_line(data = spline1df_pps, aes(y = predicted, x = xlocs, group = type, color = type)) + 
  # scale_colour_manual(values = c("PENCOMP srs"= "#4ea9c2",
  #                                "PENCOMP pps" = "#4ea9c2",
  #                                "PENCOMP-WFPBB srs" = "#e22f3d",
  #                                "PENCOMP-WFPBB pps" = "#e22f3d")) +
  labs(title = "Treatment Group 1- NO INT, Shared", subtitle = "Data points from 1 PPS sample",
       y= "Y", x = "True propensity score", color = NULL, size = "Selection weights") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_blank(),
        legend.key = element_blank())
dev.off()
pdf(file = paste0("results-clean/coef-trt0-",type_outcome, type_overlap, ".pdf"),width = 8*0.8, height = 6*0.8)

ggplot(data=samp0_pps,aes(y=Ysamp, x=pstar))+
  geom_point(alpha=0.2,color="gray", aes(size=wts)) + 
  # geom_point(data = spline0df, aes(y = predicted, x = xlocs, group = type, color = type)) + 
  geom_line(data = spline0df_pps, aes(y = predicted, x = xlocs, group = type, color = type)) + 
  labs(title = "Treatment Group 0- NO INT, Shared", subtitle = "Data points from 1 PPS sample",
       y= "Y", x = "True propensity score", color = NULL, size = "Selection weights") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background = element_blank(),
        legend.key = element_blank())
dev.off()

##  ****************************
## UNCOMMENT TO RE-SAVE RESULTS
# saveRDS(coef_pencompsrs, "nogbshared-coefres-raw/coef_pencompsrs.RDS")
# saveRDS(coef_wfpbbpensrs, "nogbshared-coefres-raw/coef_wfpbbpencompsrs.RDS")
# saveRDS(coef_pencomppps, "nogbshared-coefres-raw/coef_pencomppps.RDS")
# saveRDS(coef_wfpbbpenpps, "nogbshared-coefres-raw/coef_wfpbbpencomppps.RDS")

## Create tables ####
coef_pencompsrs <- readRDS("nogbshared-coefres-raw/coef_pencompsrs.RDS")
coef_wfpbbpensrs <- readRDS("nogbshared-coefres-raw/coef_wfpbbpencompsrs.RDS")
coef_pencomppps <- readRDS("nogbshared-coefres-raw/coef_pencomppps.RDS")
coef_wfpbbpenpps <- readRDS("nogbshared-coefres-raw/coef_wfpbbpencomppps.RDS")

resfixed0 <- data.frame(
  coef = c("int", "X2", "Z", "Phat"),
  srspencomp = coef_pencompsrs[[1]]$fixed,
  srswfpbbpencomp = coef_wfpbbpensrs$fixed0,
  ppspencomp = coef_pencomppps[[1]]$fixed,
  ppswfpbbpenpps = coef_wfpbbpenpps$fixed0
)
resfixed1 <- data.frame(
  coef = c("int", "X2", "Z", "Phat"),
  srspencomp = coef_pencompsrs[[2]]$fixed,
  srswfpbbpencomp = coef_wfpbbpensrs$fixed1,
  ppspencomp = coef_pencomppps[[2]]$fixed,
  ppswfpbbpenpps = coef_wfpbbpenpps$fixed1
)

resrand0 <- data.frame(
  coef = 1:10,
  srspencomp = coef_pencompsrs[[1]]$random$all %>% as.numeric,
  srswfpbbpencomp = coef_wfpbbpensrs$rand0 %>% as.numeric,
  ppspencomp = coef_pencomppps[[1]]$random$all %>% as.numeric,
  ppswfpbbpenpps = coef_wfpbbpenpps$rand0 %>% as.numeric
)

resrand1 <- data.frame(
  coef = 1:10,
  srspencomp = coef_pencompsrs[[2]]$random$all %>% as.numeric,
  srswfpbbpencomp = coef_wfpbbpensrs$rand1 %>% as.numeric,
  ppspencomp = coef_pencomppps[[2]]$random$all %>% as.numeric,
  ppswfpbbpenpps = coef_wfpbbpenpps$rand1 %>% as.numeric
)

