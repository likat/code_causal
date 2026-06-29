library(nlme)
library(dplyr)
library(survey)
library(ggplot2)
library(sampling)
library(cmdstanr)
library("rootSolve")
library(LaplacesDemon)

type_y <- "cts"
## CHANGE simulation type (gb / nogb)
type_outcome = "gb"
## CHANGE simulation type (shared / notshared)
type_overlap = "notshared"
pathtocode = "CODE_final/"
source(paste0(pathtocode,"helpers.R"))
source(paste0(pathtocode, "calculate_sate.R"))
source(paste0(pathtocode, "calculate_sate_x.R"))
source(paste0(pathtocode, "calculate_pencomp_cts.R"))
source(paste0(pathtocode, "calculate_wfpbb_pencomp_cts.R"))
source(paste0(pathtocode, "imputeF_cts.R"))
source(paste0(pathtocode, "gcomputeFunc_sim.R"))
source(paste0(pathtocode, "calculate_pencomp_nospline.R"))
source(paste0(pathtocode, "calculate_wtd_2stage.R"))
source(paste0(pathtocode, "STEP0-poptingting.R"))

npop <- nrow(pop)
outcome_mod <- NA
propensity_mod <- NA
x_varnames = NA
y_varname = "Ysamp"
x_name <- c("X1")
x_name_mod <- c("X1")
z_name <- "Z"

### SIMULATION PARS
sim = 200
numboot=200
set.seed(42)
nknot = 15
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
# for(i in 1:sim){

# sampsize <- rep(NA, sim)
# for(j in 1:sim){
#   set.seed(j)
#   pop$I <- NULL
#   pop$I <- rbinom(n=npop, size=1, prob = pop$pselect)
#   sampsize[j] <- sum(pop$I)
# }
# saveRDS(sampsize,paste0("results-clean/sampsize_", type_outcome, type_overlap, ".RDS"))

i <- Sys.getenv("SLURM_ARRAY_TASK_ID")
#i=1
i<- as.numeric(i)

# draw the sample
set.seed(i)
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

## DELETE LATER: check propensity score overlap
# ggplot(samp, aes(x=ptrt, fill=factor(A))) + 
#   geom_histogram(alpha=0.2, position="identity") + 
#   labs(x="True propensity score", fill="Trt group")+
#   theme_minimal()

### bootstrap containers
wtdest_pencomp_boot <- matrix(NA,nrow=numboot, ncol=2)
wtdest_pencompwfpbb_boot <- matrix(NA,nrow=numboot, ncol=2)
wtdest_nospline_boot <- matrix(NA,nrow=numboot, ncol=2)
ate_pencomp_boot <- matrix(NA,nrow=numboot, ncol=2)
ate_pencompwfpbb_boot <- matrix(NA,nrow=numboot, ncol=2)
ate_nospline_boot <- matrix(NA,nrow=numboot, ncol=2)
twostage_pencomp_boot <- matrix(NA,nrow=numboot, ncol=2)
twostage_pencompwfpbb_boot <- matrix(NA,nrow=numboot, ncol=2)
twostage_nospline_boot <- matrix(NA,nrow=numboot, ncol=2)
synthpop_pencompwfpbb_boot <- matrix(NA, nrow = numboot, ncol=2)
aiptw_boot <- rep(NA, numboot)
coef_pencomp_boot <- vector(mode="list", length = 4) ## extra 2 for cts for sigma_y
coef_pencompwfpbb_boot <- vector(mode="list", length = 4)

# BOOTSTRAP FOR FIRST STAGE IMPUTATION


boot <- 1
while(boot <= numboot){ ## START bootstrap loop
  
  # (1a) Bootstrap the sample (simple stratified bootstrap for non-WFPBB methods)
  bootind_cnt <- sample(ind_cnt, size=size_cnt, replace=T)
  bootind_trt <- sample(ind_trt, size=size_trt, replace=T)
  boot_all <- samp[c(bootind_cnt,bootind_trt),]
  boot_all$wts <- npop * normalize(boot_all$wts)
  

  # # step 1 imputation ===
  imputedlist_nospline <- calculate.pencomp.nospline(outc_model = outcome_mod)
  #   # Calculate ATE
  ate_nospline_boot[boot,]=c(mean(imputedlist_nospline[[2]]-imputedlist_nospline[[1]]), 
                             var(imputedlist_nospline[[2]]-imputedlist_nospline[[1]])/nrow(samp))
  wtdres_nospline <- calculate.wtd(imputedlist = imputedlist_nospline, bootdf=samp)
  wtdest_nospline_boot[boot,] <- c(wtdres_nospline[[1]],
                                   wtdres_nospline[[2]])
  
  print(paste0("Boot is ",boot))
  boot <- boot + 1
  bootflag <- 0
  # print(boot)
} ## END bootstrap

# step 3 rubin's combining rules and store results====
ate_nospline <- Map(c, ate_nospline, combine.rubin(ate_nospline_boot))
wtdest_nospline <- Map(c, wtdest_nospline, combine.rubin(wtdest_nospline_boot))

print(paste0("I is ",i))

# }
res <- list(patetrue = pate_true,
            ate_nospline = ate_nospline,
            wtdest_nospline = wtdest_nospline,
            bootflagct = bootflagct
)
# saveRDS(res, file=paste0("results/results_tingting_", type_outcome, type_overlap, ".RDS"))
saveRDS(res, file = paste0("results-raw-2stage/linimpute/", type_outcome, type_overlap,"/results_tingting_noaiptw", type_outcome, type_overlap,i, ".RDS"))
