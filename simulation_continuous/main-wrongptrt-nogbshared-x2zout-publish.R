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

pathtocode = "./CODE-temp/"
source(paste0(pathtocode,"helpers.R"))
source(paste0(pathtocode, "calculate_pencomp_cts.R"))
source(paste0(pathtocode, "calculate_wfpbb_pencomp_cts.R"))
source(paste0(pathtocode, "calculate_wfpbb_pencomp_selspl_cts.R"))
source(paste0(pathtocode, "calculate_ptrtspl.R"))
# source(paste0(pathtocode, "imputeF_stan_cts.R"))
source(paste0(pathtocode, "imputeF_cts.R"))

source(paste0(pathtocode, "run_AIPTW_simulation_cts.R"))
source(paste0(pathtocode, "gcomputeFunc_sim.R"))
source(paste0(pathtocode, "calculate_pencomp_nospline.R"))
source(paste0(pathtocode, "calculate_wtd.R"))
source(paste0("STEP0-poptingting.R"))
mod <- cmdstan_model(paste0(pathtocode, "outcome_mod_logistic.stan"))

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
numboot=50
set.seed(42)

## OUTCOME GENERATION
#***************************************************************
a0 <- 0
ax1 <- 1
ax2 <- 1
ax1x2 <- 0.5
az <- -0.75
# ax1z <- 0.3

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
    propensity_mod <- as.formula("A ~ X1 + X2 + X1:X2")
    propensity_mod_correct <- as.formula("A ~ X1 + X2 + Z + X1:X2")
  # propensity_mod <- as.formula("A ~ X1")
  # propensity_mod_correct <- as.formula("A ~ X1 + Z + Z:X1")
  
  pop$ptrt <-  expit(ax1*pop$X1 +ax2*pop$X2 + ax1x2*pop$X1*pop$X2+
                       az*pop$Z )

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
wtdest_noz <- vector(mode="list", length=8)
wtdest_correct <- vector(mode="list", length=8)
wtdest_pselspl <- vector(mode="list", length=8)

synthpop_noz <- vector(mode="list", length = 8)
synthpop_correct<- vector(mode="list", length = 8)
synthpop_pselspl<- vector(mode="list", length = 8)

twostage_noz <- vector(mode="list", length=8)
twostage_correct<- vector(mode="list", length=8)
twostage_pselspl<- vector(mode="list", length=8)

names(wtdest_noz) <- resnames
names(wtdest_correct) <- resnames
names(wtdest_pselspl) <- resnames

names(synthpop_noz) <- resnames
names(synthpop_correct) <- resnames
names(synthpop_pselspl) <- resnames

names(twostage_noz) <- resnames
names(twostage_correct) <- resnames
names(twostage_pselspl) <- resnames

bootflag <- 0
bootflagct <- 0
## To replicate a true experiment, recalculate the outcome every sample and treatment assignment
pop$popID <- 1:nrow(pop)
# for(i in 1:sim){
 i <- Sys.getenv("SLURM_ARRAY_TASK_ID")
# i=1
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
  
  ### SATE
  ## CHECK: for the test outcome case, is the sate different from true pate?
  ggplot(samp, aes(x=ptrt, fill=factor(A))) + 
      geom_histogram(alpha=0.2, position="identity") +
      labs(x="True propensity score", fill="Trt group")+
      theme_minimal()
  ### bootstrap containers
  wtdest_noz_boot <- matrix(NA,nrow=numboot, ncol=2)
  wtdest_correct_boot <- matrix(NA,nrow=numboot, ncol=2)
  wtdest_pselspl_boot <- matrix(NA,nrow=numboot, ncol=2)
  twostage_noz_boot <- matrix(NA,nrow=numboot, ncol=2)
  twostage_correct_boot <- matrix(NA,nrow=numboot, ncol=2)
  twostage_pselspl_boot <- matrix(NA,nrow=numboot, ncol=2)
  
  synthpop_noz_boot <- matrix(NA, nrow = numboot, ncol=2)
  synthpop_correct_boot <- matrix(NA, nrow = numboot, ncol=2)
  synthpop_pselspl_boot <- matrix(NA, nrow = numboot, ncol=2)
  
  aiptw_boot <- rep(NA, numboot)
  coef_pencomp_boot <- vector(mode="list", length = 4) ## extra 2 for cts for sigma_y
  coef_pencompwfpbb_boot <- vector(mode="list", length = 4)

  # BOOTSTRAP FOR FIRST STAGE IMPUTATION
  boot <- 1
  while(boot <= numboot){ ## START bootstrap loop
    # (1) Bootstrap the sample (stratified by treatment group)
    bootflag <- 0
    bootind_cnt <- sample(ind_cnt, size=size_cnt, replace=T)
    bootind_trt <- sample(ind_trt, size=size_trt, replace=T)
    boot_all <- samp[c(bootind_cnt,bootind_trt),]
    boot_all$wts <- npop * normalize(boot_all$wts)
    
    # WFPBB bootstrap
    ind_root <- c(ind_cnt, ind_trt) ## must be same order index as the bb weights 
    bbst0 <- length(ind_cnt) * normalize(BayesianBootstrap(ind_cnt, 1))
    bbst1 <- length(ind_trt) * normalize(BayesianBootstrap(ind_trt, 1))
    bbwt_norm <- nrow(samp) * normalize(c(bbst0, bbst1))
    bbdraw = rmultinom(n=1,size=nrow(samp), prob = bbwt_norm)
    
    ind_bb <- rep(ind_root, bbdraw)
    boot_allbb <- samp[ind_bb,]
    uniqueind_bb <- ind_root[bbdraw != 0]
    wts_bb <- npop * normalize(samp$wts[uniqueind_bb] * bbdraw[bbdraw != 0])
    # samp$bbwt[c(ind_cnt, ind_trt)] <- bbwt_norm

    # # step 1 imputation ===
    t0 <- Sys.time()
    imputedlist_noz <- calculate.wfpbb.pencomp(propen_model=propensity_mod,
                                                       outc_model = outcome_mod,
                                                       y.varname = y_varname,
                                                       uniqueind = uniqueind_bb,
                                                       polyawts = wts_bb,
                                                       trt.varname = "A",
                                                       num.knot = 10,
                                                       F_draw=10)
    t1 <- Sys.time()
    if(sum(is.na(imputedlist_noz)) > 0){bootflagct = bootflagct + 1; next}
    print("test1")
    imputedlist_correct <- calculate.wfpbb.pencomp(propen_model=propensity_mod_correct,
                                                       outc_model = outcome_mod,
                                                       y.varname = y_varname,
                                                       uniqueind = uniqueind_bb,
                                                       polyawts = wts_bb,
                                                       trt.varname = "A",
                                                       num.knot = 10,
                                                       F_draw=10)
    if(sum(is.na(imputedlist_correct)) > 0){bootflagct = bootflagct + 1; next}
    print("test2")
    t0_spl <- Sys.time()
    imputedlist_pselspl <- calculate.pwfpbb.pselspl(propen_model=propensity_mod,
                                                       outc_model = outcome_mod,
                                                       y.varname = y_varname,
                                                       uniqueind = uniqueind_bb,
                                                       polyawts = wts_bb,
                                                       trt.varname = "A",
                                                       num.knot = 10,
                                                       F_draw=10)
    t1_spl <- Sys.time()
    if(sum(is.na(imputedlist_pselspl)) > 0){bootflagct = bootflagct + 1; next}
    print("test3")
    
    imputedlist_noz$wtdest
    imputedlist_correct$wtdest
    imputedlist_pselspl$wtdest
    
    # aiptw_boot[boot] <- aiptw_bootres$result
    synthpop_noz_boot[boot,] <- c(imputedlist_noz$pate_gen_synthpop, imputedlist_noz$varpategen_synthpop)
    synthpop_correct_boot[boot,] <- c(imputedlist_correct$pate_gen_synthpop, imputedlist_correct$varpategen_synthpop)
    synthpop_pselspl_boot[boot,] <- c(imputedlist_pselspl$pate_gen_synthpop, imputedlist_pselspl$varpategen_synthpop)
    

    # print(ate_nospline_boot[boot,])
    # step 2 generalize! ====

    wtdest_noz_boot[boot,] <- c(imputedlist_noz$wtdest)
    wtdest_correct_boot[boot,] <- c(imputedlist_correct$wtdest)
    wtdest_pselspl_boot[boot,] <- c(imputedlist_pselspl$wtdest)
    
    ## how to generalize with weights? If the sample is not a probability sample from population, then we would have to estimate them here.
    twostage_noz_boot[boot,] <- c(imputedlist_noz$pate_gen, imputedlist_noz$varpategen)
    twostage_correct_boot[boot,] <- c(imputedlist_correct$pate_gen, imputedlist_correct$varpategen)
    twostage_pselspl_boot[boot,] <- c(imputedlist_pselspl$pate_gen, imputedlist_pselspl$varpategen)
    
    ## collect coefficient estimates
    print(paste0("Boot is ",boot))
    boot <- boot + 1
    bootflag <- 0
    # print(boot)
  } ## END bootstrap
  
  # step 3 rubin's combining rules and store results====
  wtdest_noz <- Map(c, wtdest_noz, combine.rubin(wtdest_noz_boot))
  wtdest_correct <- Map(c, wtdest_correct, combine.rubin(wtdest_correct_boot))
  wtdest_pselspl <- Map(c, wtdest_pselspl, combine.rubin(wtdest_pselspl_boot))
  
  twostage_noz <- Map(c, twostage_noz, combine.rubin(twostage_noz_boot))
  twostage_correct <- Map(c, twostage_correct, combine.rubin(twostage_correct_boot))
  twostage_pselspl <- Map(c, twostage_pselspl, combine.rubin(twostage_pselspl_boot))
  
  synthpop_noz <- Map(c, synthpop_noz, combine.rubin(synthpop_noz_boot))
  synthpop_correct <- Map(c, synthpop_correct, combine.rubin(synthpop_correct_boot))
  synthpop_pselspl <- Map(c, synthpop_pselspl, combine.rubin(synthpop_pselspl_boot))
  
  print(paste0("I is ",i))
  
# }
res <- list(patetrue = pate_true, sate = sate,
            wtdest_noz = wtdest_noz,
            wtdest_correct = wtdest_correct,
            wtdest_pselspl= wtdest_pselspl,
            twostage_noz = twostage_noz,
            twostage_correct = twostage_correct,
            twostage_pselspl = twostage_pselspl,
            synthpop_noz = synthpop_noz,
            synthpop_correct = synthpop_correct,
            synthpop_pselspl = synthpop_pselspl,
            # time_aiptw = time_aiptw,
            # aiptw = aiptw,
            bootflagct = bootflagct
            )
# saveRDS(res, file=paste0("results/results_tingting_", type_outcome, type_overlap, ".RDS"))
saveRDS(res, file = paste0("results-raw-wrongptrt/results_wrongptrt_x2z", type_outcome, type_overlap,i,".RDS"))
