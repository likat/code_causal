library(nlme)
library(dplyr)
library(survey)
library(ggplot2)
library(sampling)
library(cmdstanr)
library("rootSolve")

type_y <- "cts"
type_outcome = "nogb"
type_overlap = "shared"

## Change to directory containing helper code files
pathtocode = "./CODE-temp/"
source(paste0(pathtocode,"helpers.R"))
source(paste0(pathtocode, "run_AIPTW_simulation_cts.R"))
source(paste0(pathtocode, "gcomputeFunc_sim.R"))
source(paste0( "STEP0-poptingting.R"))

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
    mu1 = bA + bx1*X1 + bAx1*X1 + bx2*X2 + bz*Z+ bAz*Z,
    mu0 = bx1*X1 + bx2*X2 + bz*Z
    )
  print(type_outcome)
}else if(type_outcome == "nogb"){
#  pop <- pop %>% mutate(
#    mu1 = bA + bx1*X1  + bz*Z+ bAx1*X1,
#    mu0 = bx1*X1 + bz*Z,
#  )
  pop <- pop %>% mutate(
    mu1 = bA + bx1*X1 + bAx1*X1 + bx2*X2 + bz*Z,
    mu0 = bx1*X1 + bx2*X2 + bz*Z
    )
  print(type_outcome)
  
}
  pop$Y1 <- rnorm(n=npop, pop$mu1)
  pop$Y0 <- rnorm(n=npop, pop$mu0)
  
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
# outcome_mod <- "Ysamp ~ X1 + X2 + Z"
# outcome_mod_nospline <- formulaF(varList=varnames_outcome_nospline, y.name=y_varname)
outcome_mod_generalize <- formulaF(varList=varnames_outcome_generalize, y.name=y_varname)

## PARAMETERS FOR AIPTW
  covariateXnames <- c("X1","X2", "Z")
  propensity_mod_aiptw <- propensity_mod
  outcome_mod_aiptw <- outcome_mod
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

coef_pencomp <- vector(mode = "list", length = 6)
coef_pencompwfpbb <- vector(mode="list", length = 6)

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
i <- Sys.getenv("SLURM_ARRAY_TASK_ID")
#i=1
i<- as.numeric(i)

  # draw the sample
  set.seed(i)
  pop$I <- NULL
  pop$I <- rbinom(n=npop, size=1, prob = pop$pselect)
  # pop$I <- rbinom(n=npop, size=1, prob = nsamp/nrow(pop))
  
  samp <- pop[pop$I==1,]
  samp$wts <- npop*normalize(1/samp$pselect)
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
  

  ## DELETE LATER: check propensity score overlap
  # ggplot(samp, aes(x=ptrt, fill=factor(A))) + 
  #   geom_histogram(alpha=0.2, position="identity") + 
  #   labs(x="True propensity score", fill="Trt group")+
  #   theme_minimal()
  
  
  ### SATE
  ## CHECK: for the test outcome case, is the sate different from true pate?

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
  coef_pencomp_boot <- vector(mode="list", length = 6) ## extra 2 for cts for sigma_y
  coef_pencompwfpbb_boot <- vector(mode="list", length = 6)

  # BOOTSTRAP FOR FIRST STAGE IMPUTATION
  boot <- 1
  while(boot <= numboot){ ## START bootstrap loop
    # (1) Bootstrap the sample (stratified by treatment group)
    bootflag <- 0
    
    bootind_cnt <- sample(ind_cnt, size=size_cnt, replace=T)
    bootind_trt <- sample(ind_trt, size=size_trt, replace=T)
    boot_all <- samp[c(bootind_cnt,bootind_trt),]
    boot_all$wts <- npop * normalize(boot_all$wts)
    
    #   
    # # step 1 imputation ===
    # imputedlist_nospline <- calculate.pencomp.nospline(outc_model = outcome_mod_nospline)
    if(bootflag == 1){bootflagct = bootflagct + 1; next}
    aiptw_bootres <- AIPTW(dataInput=boot_all, 
                           covariateXnames= covariateXnames,
                           sampwts = boot_all$wts,
                           treat.varname= "A", 
                           outcome.varname= y_varname, 
                           propen.model2a = propensity_mod_aiptw,
                           modely1 = outcome_mod_aiptw,
                           modely0 = outcome_mod_aiptw,
                           numCarlo=ncarlo_aiptw)

    ## Store results 
    time_aiptw <- time_aiptw + aiptw_bootres$time/numboot

    aiptw_boot[boot] <- aiptw_bootres$result
    print(paste0("Boot is ",boot))
    boot <- boot + 1
    bootflag <- 0
    print(boot)
  } ## END bootstrap
  

  aiptw <- Map(c, aiptw, combine.aiptw(aiptw_boot))
  print(paste0("I is ",i))
  
# }
res <- list(patetrue = pate_true, 
            aiptw = aiptw
            )
saveRDS(res, file = paste0("results-raw/results_tingting_aiptw", type_outcome, type_overlap,i,".RDS"))
