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
terms_outcome = "x2z"
type_estimator = "pselspl"
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
source(paste0(pathtocode, "CODE_cts/calculate_wfpbb_pencomp_selspl_NEW.R"))
source(paste0(pathtocode, "CODE/calculate_ptrtspl.R"))

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
mod <- cmdstan_model(paste0(pathtocode, "CODE/outcome_mod_logistic.stan"))

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
numboot=100
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
wtdest_pencomp <- vector(mode="list", length=8)
wtdest_pencompwfpbb <- vector(mode="list", length=8)
wtdest_nospline <- vector(mode="list", length=8)

aiptw <- vector(mode = "list", length = 4)
synthpop_pencompwfpbb <-synthpop_pselspl <- synthpop_correct <-  vector(mode="list", length = 8)

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


### SATE
## CHECK: for the test outcome case, is the sate different from true pate?

# AIPW(
#   dataInput=samp, 
#   propen.mod = propensity_mod_aiptw,
#   modely1 = outcome_mod,
#   modely0 = outcome_mod,
#   outcomefamily = "gaussian")
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
synthpop_pencompwfpbb_boot <- synthpop_pselspl_boot <- synthpop_correct_boot <- matrix(NA, nrow = numboot, ncol=2)
aiptw_boot <- rep(NA, numboot)
coef_pencomp_boot <- vector(mode="list", length = 4) ## extra 2 for cts for sigma_y
coef_pencompwfpbb_boot <- vector(mode="list", length = 4)
time_pselspl <- rep(NA, numboot)
# BOOTSTRAP FOR FIRST STAGE IMPUTATION


boot <- 1
while(boot <= numboot){ ## START bootstrap loop
  
  # (1a) Bootstrap the sample (simple stratified bootstrap for non-WFPBB methods)
  bootind_cnt <- sample(ind_cnt, size=size_cnt, replace=T)
  bootind_trt <- sample(ind_trt, size=size_trt, replace=T)
  boot_all <- samp[c(bootind_cnt,bootind_trt),]
  boot_all$wts <- npop * normalize(boot_all$wts)
  
  # (1b) Bayesian bootstrap of the sample
  cluster_stratum_ref$bbresamp <- 0
  cluster_stratum_ref$bbclusmultiplier <- 0
  
  samp$bbclusmultiplier <- NULL
  samp$bbresamp <- NULL
  ## with clusters and strata, resample clusters
  for(k in 1:nstrat){
    cluspool <- cluster_stratum_ref$cluster[cluster_stratum_ref$stratum==k]
    tempnumclus <- length(cluspool)
    bootstrap <- BayesianBootstrap(cluspool, 1)
    resamp <- rmultinom(1, tempnumclus - 1, bootstrap)
    cluster_stratum_ref$bbresamp[cluster_stratum_ref$cluster %in% cluspool] <- resamp
  }
  cluster_stratum_ref$bbclusmultiplier <- cluster_stratum_ref$bbresamp*cluster_stratum_ref$mstar
  
  samp <- samp %>% left_join(cluster_stratum_ref %>% select(cluster, stratum,bbclusmultiplier, bbresamp))
  boot_allbb <- samp[samp$bbresamp != 0,]
  
  wts_bb <- npop * normalize(samp$wts*samp$bbclusmultiplier)
  wts_bb <- wts_bb[which(samp$bbresamp != 0)]
  uniqueind_bb <- which(samp$bbresamp != 0)
  bootflag <- 0
  # bootind_cnt <- sample(ind_cnt, size=size_cnt, replace=T)
  # bootind_trt <- sample(ind_trt, size=size_trt, replace=T)
  # boot_all <- samp[c(bootind_cnt,bootind_trt),]
  # boot_all$wts <- npop * normalize(boot_all$wts)
  
  # WFPBB bootstrap
  # ind_root <- c(ind_cnt, ind_trt) ## must be same order index as the bb weights 
  # bbst0 <- length(ind_cnt) * normalize(BayesianBootstrap(ind_cnt, 1))
  # bbst1 <- length(ind_trt) * normalize(BayesianBootstrap(ind_trt, 1))
  
  # bbwt_norm <-BayesianBootstrap(1:nrow(samp), 1)
  # bbdraw = rmultinom(n=1,size=nrow(samp), prob = bbwt_norm)
  # ind_bb <- rep(1:nrow(samp), bbdraw)
  # boot_allbb <- samp[ind_bb,]
  # uniqueind_bb <- ind_root[bbdraw != 0]
  # wts_bb <- npop * normalize(samp$wts[uniqueind_bb] * bbdraw[bbdraw != 0])
  # samp$bbwt[c(ind_cnt, ind_trt)] <- bbwt_norm
  
  # # step 1 imputation ===
  # imputedlist_nospline <- calculate.pencomp.nospline(outc_model = outcome_mod_nospline)
  # if(bootflag == 1){bootflagct = bootflagct + 1; next}
  t0_wfpbb <- Sys.time()
  imputedlist_pselspl <- calculate.pwfpbb.pselspl(propen_model=propensity_mod,
                                                  outc_model = outcome_mod,
                                                  y.varname = y_varname,
                                                  uniqueind = uniqueind_bb,
                                                  polyawts = wts_bb,
                                                  trt.varname = "A",
                                                  num.knot = nknot,
                                                  F_draw=20)
  # 
  t1_wfpbb <- Sys.time()
  time_pselspl[boot] <- imputedlist_pselspl$time_pencompwfpbb
  # if(sum(is.na(imputedlist_pselspl)) > 0){bootflagct = bootflagct + 1; next}
  # t0 <- Sys.time()
  # imputedlist_wfpbbspline <- calculate.wfpbb.pencomp(propen_model=propensity_mod,
  #                                                    outc_model = outcome_mod,
  #                                                    y.varname = y_varname,
  #                                                    uniqueind = uniqueind_bb,
  #                                                    polyawts = wts_bb,
  #                                                    trt.varname = "A",
  #                                                    num.knot = nknot,
  #                                                    F_draw=20)
  # t1 <- Sys.time()
  # if(sum(is.na(imputedlist_wfpbbspline)) > 0){bootflagct = bootflagct + 1; next}

  # imputedlist_correct <- calculate.wfpbb.pencomp(propen_model=propensity_mod_correct,
  #                                                outc_model = outcome_mod,
  #                                                y.varname = y_varname,
  #                                                uniqueind = uniqueind_bb,
  #                                                polyawts = wts_bb,
  #                                                trt.varname = "A",
  #                                                num.knot = nknot,
  #                                                F_draw=20)
  # if(sum(is.na(imputedlist_correct)) > 0){bootflagct = bootflagct + 1; next}
  # aiptw_bootres <- AIPTW(dataInput=samp, 
  #                        covariateXnames= covariateXnames,
  #                        sampwts = samp$wts,
  #                        treat.varname= "A", 
  #                        outcome.varname= y_varname, 
  #                        propen.model2a = propensity_mod_aiptw,
  #                        modely1 = outcome_mod_aiptw,
  #                        modely0 = outcome_mod_aiptw,
  #                        numCarlo=ncarlo_aiptw)
  # if(sum(c(is.na(imputedlist_wfpbbspline)))>0){bootflagct = bootflagct + 1; next}
  # if(sum(c(is.na(imputedlist_spline), is.na(imputedlist_wfpbbspline)))>0){print("error!")}
  #if(bootflag == 1){bootflagct = bootflagct + 1; next}
  
  ## Store results 
  # time_aiptw <- time_aiptw + aiptw_bootres$time/numboot

  # aiptw_boot[boot] <- aiptw_bootres$result
  # synthpop_pencompwfpbb_boot[boot,] <- c(imputedlist_wfpbbspline$pate_gen_synthpop, imputedlist_wfpbbspline$varpategen_synthpop)
  synthpop_pselspl_boot[boot,] <- c(imputedlist_pselspl$pate_gen_synthpop, imputedlist_pselspl$varpategen_synthpop)
  # synthpop_correct_boot[boot,] <- c(imputedlist_correct$pate_gen_synthpop, imputedlist_correct$varpategen_synthpop)
  
  print(paste0("Boot is ",boot))
  boot <- boot + 1
  bootflag <- 0
  # print(boot)
} ## END bootstrap

# step 3 rubin's combining rules and store results====
# synthpop_pencompwfpbb <- Map(c, synthpop_pencompwfpbb, combine.rubin(synthpop_pencompwfpbb_boot))
# synthpop_correct <- Map(c, synthpop_correct, combine.rubin(synthpop_correct_boot))
synthpop_pselspl <- Map(c, synthpop_pselspl, combine.rubin(synthpop_pselspl_boot))

# aiptw <- Map(c, aiptw, combine.aiptw(aiptw_boot))
print(paste0("I is ",i))

# }
res <- list(patetrue = pate_true, sate = sate,
            # synthpop_noz = synthpop_pencompwfpbb,
            # synthpop_correct = synthpop_correct,
            synthpop_pselspl = synthpop_pselspl,
            # aiptw = aiptw,
            time = mean(time_pselspl),
            bootflagct = bootflagct
)
saveRDS(res, file = paste0("results-raw-2stage/wrongptrt/results_wrongptrt_", terms_outcome,"_" ,type_outcome, type_overlap,"_", type_estimator,i,".RDS"))
# saveRDS(res, file = paste0("results-raw-2stage/", type_outcome, type_overlap,"/results_tingting_noaiptw", type_outcome, type_overlap,i, ".RDS"))
