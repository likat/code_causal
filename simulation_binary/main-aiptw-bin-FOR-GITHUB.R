library(nlme)
library(dplyr)
library(survey)
library(ggplot2)
library(sampling)
library(cmdstanr)
library("rootSolve")
library(LaplacesDemon)

type_y <- "bin"
type_outcome = "gb" ## CHANGE simulation type here (gb / nogb)
type_overlap = "notshared" ## CHANGE simulation type here (shared / notshared)
pathtocode = "./CODE_final/"
source(paste0(pathtocode,"helpers.R"))
source(paste0(pathtocode, "run_AIPTW_new.R"))
source(paste0(pathtocode, "STEP0-poptingting-bin.R"))

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

## OUTCOME GENERATION
#***************************************************************
a0 <- 0
ax1 <- 1
ax2 <- 1
ax1x2 <- 0.5
az <- -0.75
b0 <- 0
bx1 <- 0.2 # X1
bx2 <- 0.1 # X2
bA <- 0.5
bAx1 <- 0.8
bAx2 <- 0.5

bz <- -0.2
bAz <- 0.4
#***************************************************************
if(type_outcome == "gb"){
  pop <- pop %>% mutate(
    mu1 = cluseff + bA + bx1*X1 + bAx1*X1 + bx2*X2 + bz*Z+ bAz*Z,
    mu0 = cluseff + bx1*X1 + bx2*X2 + bz*Z
  )
  print(type_outcome)
}else if(type_outcome == "nogb"){
  pop <- pop %>% mutate(
    mu1 = cluseff +bA + bx1*X1 + bAx1*X1 + bx2*X2 + bz*Z,
    mu0 = cluseff +bx1*X1 + bx2*X2 + bz*Z
  )
  print(type_outcome)
  
}
pop$Y1 <- rbinom(n = nrow(pop), size = 1, prob = expit(pop$mu1))
pop$Y0 <- rbinom(n = nrow(pop), size = 1, prob = expit(pop$mu0))

## PROPENSITY GENERATION
if(type_overlap == "shared"){
  propensity_mod <- as.formula("A ~ X1 + X2 + Z + X1:X2")
  pop$ptrt <-  expit(ax1*pop$X1 + ax2*pop$X2 + ax1x2*pop$X1*pop$X2 +
                       az*pop$Z)
  
}else if(type_overlap == "notshared"){
  propensity_mod <- as.formula("A ~ X1 + X2 + X1:X2")
  pop$ptrt <-  expit(ax1*pop$X1 + ax2*pop$X2 + ax1x2*pop$X1*pop$X2 )
}

x_varnames <- c(x_name_mod, z_name, paste0(x_name_mod,":",z_name))

varnames_outcome <-  c("X2", "Z")
varnames_outcome_generalize <- c(x_name_mod, z_name, paste0(x_name_mod, ":", z_name))
outcome_mod <- as.formula("Ysamp ~ X2 + Z")
outcome_mod_correct <- as.formula("Ysamp ~ X1 + X2 + Z")
outcome_mod_generalize <- formulaF(varList=varnames_outcome_generalize, y.name=y_varname)

## PARAMETERS FOR AIPTW
covariateXnames <- c("X1","X2", "Z")
propensity_mod_aipw <- propensity_mod
propensity_mod_mis <- as.formula("A~X1 + X2 + X1:X2")
outcome_mod_aipw <- outcome_mod
outcome_mod_aipw_noz <- as.formula("Ysamp ~ X2 + X1")


# true PATE
pate_true <- mean(pop$Y1 - pop$Y0)

### ESTIMATE CONTAINERS
sate <- rep(NA, sim)
resnames <- c("est","withinvar","btwnvar","mivar","ci_tlength","ci_emplength","cicov_t","cicov_emp") # must match names from rubin function
resnames_aipw <- c("est", "var", "cilower","ciupper")
time_aipw <- 0
time_pencompwfpbb <- 0

aipw_unwtd_correcttrt<-aipw_wrong <-aipw_wtd_wtdtrt <-aipw_wtd_mistrt <-aipw_unwtd_wtdtrt <- vector(mode = "list", length = 4)
aipw_unwtd_correcttrt_noz <-aipw_wrong_noz<-aipw_wtd_wtdtrt_noz <-aipw_wtd_mistrt_noz <-aipw_unwtd_wtdtrt_noz <- vector(mode = "list", length = 4)
aipw_correct_bench  <-aipw_ptrtwrong_bench <-aipw_outwrong_bench  <- vector(mode = "list", length = 4)

names(aipw_unwtd_correcttrt) <-names(aipw_wtd_wtdtrt) <- names(aipw_wtd_mistrt) <- names(aipw_unwtd_wtdtrt) <-  resnames_aipw
names(aipw_unwtd_correcttrt_noz) <-names(aipw_wtd_wtdtrt_noz) <- names(aipw_wtd_mistrt_noz) <- names(aipw_unwtd_wtdtrt_noz) <-  resnames_aipw
names(aipw_correct_bench) <-names(aipw_wrong) <-names(aipw_wrong_noz)<-names(aipw_ptrtwrong_bench) <-names(aipw_outwrong_bench)  <- vector(mode = "list", length = 4)
trtcoef_aipw <- matrix(nrow=numboot, ncol = 4)
trtcoef_aipw_unwtd <- matrix(nrow=numboot, ncol = 4)

bootflag <- 0
bootflagct <- 0
pop$popID <- 1:nrow(pop)

### BEGIN SIM
## Note: this was made for the SLURM environment. Change to for loop here if running locally.
i <- Sys.getenv("SLURM_ARRAY_TASK_ID") 
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


## DELETE LATER: check propensity score overlap
# ggplot(samp, aes(x=ptrt, fill=factor(A))) + 
#   geom_histogram(alpha=0.2, position="identity") + 
#   labs(x="True propensity score", fill="Trt group")+
#   theme_minimal()
aipw_temp <- AIPW( dataInput=samp, 
                   propen.mod = propensity_mod_aipw,
                   propen.mod.mis = propensity_mod_mis,
                   modely1 = outcome_mod,
                   modely0 = outcome_mod,
                   outcomefamily = "binomial")
aipw_est_unwtd_correcttrt = aipw_temp$pate_unwtd_correcttrt ## both treatment and outcome models weighted
aipw_est_wtd_wtdtrt = aipw_temp$pate_wtd_wtdtrt ## only outcome model weighted
aipw_est_wtd_mistrt = aipw_temp$pate_wtd_mistrt ## both models unweighted
aipw_est_unwtd_wtdtrt = aipw_temp$pate_unwtd_wtdtrt ## misspecified treatment model unweighted treatment model, weighted outcome model
aipw_est_allwrong <- aipw_temp$pate_allwrong
aipw_temp_noz <- AIPW( dataInput=samp, 
                       propen.mod = propensity_mod_aipw,
                       propen.mod.mis = propensity_mod_mis,
                       modely1 = outcome_mod_aipw_noz,
                       modely0 = outcome_mod_aipw_noz,
                       outcomefamily = "binomial")
aipw_est_unwtd_correcttrt_noz = aipw_temp_noz$pate_unwtd_correcttrt ## both treatment and outcome models weighted
aipw_est_wtd_wtdtrt_noz = aipw_temp_noz$pate_wtd_wtdtrt ## only outcome model weighted
aipw_est_wtd_mistrt_noz = aipw_temp_noz$pate_wtd_mistrt ## both models unweighted
aipw_est_unwtd_wtdtrt_noz = aipw_temp_noz$pate_unwtd_wtdtrt ## misspecified treatment model unweighted treatment model, weighted outcome model
aipw_est_allwrong_noz <- aipw_temp_noz$pate_allwrong

aipw_temp_benchmark <- AIPW( dataInput=samp,
                             propen.mod = propensity_mod_aipw,
                             propen.mod.mis = propensity_mod_mis,
                             modely1 = outcome_mod,
                             modely0 = outcome_mod,
                             # modely1 = as.formula("Ysamp~X1+X2"),
                             # modely0 = as.formula("Ysamp~X1+X2"),
                             modely1_correct = outcome_mod_correct,
                             modely0_correct = outcome_mod_correct,
                             benchmark=T,
                             outcomefamily = "binomial")
aipw_est_correct_bench = aipw_temp_benchmark$pate_correct ## both treatment and outcome models weighted
aipw_est_ptrtwrong_bench =aipw_temp_benchmark$pate_ptrtwrong ## only outcome model weighted
aipw_est_outwrong_bench = aipw_temp_benchmark$pate_outwrong ## both models unweighted, ptrtmod correct. not sure why there is significant bias

### SATE
## CHECK: for the test outcome case, is the sate different from true pate?

aipw_bt_unwtd_correcttrt <- rep(NA, numboot)
aipw_bt_wtd_wtdtrt <- rep(NA,numboot)
aipw_bt_wtd_mistrt <- rep(NA, numboot)
aipw_bt_unwtd_wtdtrt <- rep(NA, numboot)
aipw_bt_allwrong <- rep(NA,numboot)

aipw_bt_unwtd_correcttrt_noz <- rep(NA, numboot)
aipw_bt_wtd_wtdtrt_noz <- rep(NA,numboot)
aipw_bt_wtd_mistrt_noz <- rep(NA, numboot)
aipw_bt_unwtd_wtdtrt_noz <- rep(NA, numboot)
aipw_bt_allwrong_noz <- rep(NA,numboot)

aipw_bt_correct <- rep(NA, numboot)
aipw_bt_ptrtwrong <- rep(NA,numboot)
aipw_bt_outwrong <- rep(NA, numboot)
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
  aipw_bt <- AIPW( dataInput=boot_all, 
                   propen.mod = propensity_mod_aipw,
                   propen.mod.mis = propensity_mod_mis,
                   modely1 = outcome_mod,
                   modely0 = outcome_mod,
                   outcomefamily = "binomial")
  aipw_bt_unwtd_correcttrt[boot] = aipw_bt$pate_unwtd_correcttrt ## unweighted misspecified outcome model, correct treatment model
  aipw_bt_wtd_wtdtrt[boot] = aipw_bt$pate_wtd_wtdtrt ## weighted misspecified outcome model, weighted misspecified treatmnet model
  aipw_bt_wtd_mistrt[boot] = aipw_bt$pate_wtd_mistrt ## weighted misspecified outcome model, unweighted misspecified treatment model
  aipw_bt_unwtd_wtdtrt[boot] = aipw_bt$pate_unwtd_wtdtrt ## unweighted misspecified outcome model, weighted misspecified treatment model
  aipw_bt_allwrong[boot] <- aipw_bt$pate_allwrong
  
  aipw_bt_temp_noz <- AIPW( dataInput=boot_all, 
                            propen.mod = propensity_mod_aipw,
                            propen.mod.mis = propensity_mod_mis,
                            modely1 = outcome_mod_aipw_noz,
                            modely0 = outcome_mod_aipw_noz,
                            outcomefamily = "binomial")
  aipw_bt_unwtd_correcttrt_noz[boot] = aipw_bt_temp_noz$pate_unwtd_correcttrt 
  aipw_bt_wtd_wtdtrt_noz[boot] = aipw_bt_temp_noz$pate_wtd_wtdtrt
  aipw_bt_wtd_mistrt_noz[boot] = aipw_bt_temp_noz$pate_wtd_mistrt 
  aipw_bt_unwtd_wtdtrt_noz[boot] = aipw_bt_temp_noz$pate_unwtd_wtdtrt 
  aipw_bt_allwrong_noz[boot] <- aipw_bt_temp_noz$pate_allwrong
  
  # aipw_correct_bench <-aipw_ptrtwrong_bench <-aipw_outwrong_bench  <- vector(mode = "list", length = 4)
  
  aipw_bt_temp_benchmark <- AIPW( dataInput=boot_all,
                                  propen.mod = propensity_mod_aipw,
                                  propen.mod.mis = propensity_mod_mis,
                                  modely1 = outcome_mod,
                                  modely0 = outcome_mod,
                                  modely1_correct = outcome_mod_correct,
                                  modely0_correct = outcome_mod_correct,
                                  benchmark=T,
                                  outcomefamily = "binomial")
  aipw_bt_correct[boot] = aipw_bt_temp_benchmark$pate_correct ## both treatment and outcome models correct
  aipw_bt_ptrtwrong[boot] =aipw_bt_temp_benchmark$pate_ptrtwrong ## only outcome model correct, ptrt wrong
  aipw_bt_outwrong[boot] = aipw_bt_temp_benchmark$pate_outwrong ## both models unweighted, ptrtmod correct. not sure why there is significant bias
  
  ## Store results 
  time_aipw <- time_aipw + aipw_bt$time/numboot
  # time_aipw_noz <- time_aipw + aipw_bt_temp_noz$time/numboot
  # time_aipw_benchmark <- time_aipw + aipw_bt_temp_benchmark$time/numboot
  
  print(paste0("Boot is ",boot))
  boot <- boot + 1
  bootflag <- 0
  print(boot)
} ## END bootstrap

aipw_unwtd_correcttrt <- combine.aipw(est = aipw_est_unwtd_correcttrt, boottbl = aipw_bt_unwtd_correcttrt,
                                      popsize = npop, sampsize = nrow(samp)) #est, boottbl, popsize,sampsize, num_boot = numboot
aipw_unwtd_noz <- combine.aipw(est = aipw_est_unwtd_correcttrt_noz, boottbl = aipw_bt_unwtd_correcttrt_noz,
                               popsize = npop, sampsize = nrow(samp))
aipw_wtd_wtdtrt <- combine.aipw(est = aipw_est_wtd_wtdtrt, boottbl = aipw_bt_wtd_wtdtrt,
                                popsize = npop, sampsize = nrow(samp))
aipw_wtd_wtdtrt_noz <- combine.aipw(est = aipw_est_wtd_wtdtrt_noz, boottbl = aipw_bt_wtd_wtdtrt_noz,
                                    popsize = npop, sampsize = nrow(samp))
aipw_wtd_mistrt <- combine.aipw(est = aipw_est_wtd_mistrt, boottbl = aipw_bt_wtd_mistrt,
                                popsize = npop, sampsize = nrow(samp))
aipw_wtd_mistrt_noz <- combine.aipw(est = aipw_est_wtd_mistrt_noz, boottbl = aipw_bt_wtd_mistrt_noz,
                                    popsize = npop, sampsize = nrow(samp))
aipw_unwtd_wtdtrt <- combine.aipw(est = aipw_est_unwtd_wtdtrt, boottbl = aipw_bt_unwtd_wtdtrt,
                                  popsize = npop, sampsize = nrow(samp))
aipw_unwtd_wtdtrt_noz <- combine.aipw(est = aipw_est_unwtd_wtdtrt_noz, boottbl = aipw_bt_unwtd_wtdtrt_noz,
                                      popsize = npop, sampsize = nrow(samp))
aipw_correct_bench <- combine.aipw(est = aipw_est_correct_bench, boottbl = aipw_bt_correct,
                                   popsize = npop, sampsize = nrow(samp))
aipw_ptrtwrong_bench <- combine.aipw(est = aipw_est_ptrtwrong_bench, boottbl = aipw_bt_ptrtwrong,
                                     popsize = npop, sampsize = nrow(samp))
aipw_outwrong_bench <- combine.aipw(est = aipw_est_outwrong_bench, boottbl = aipw_bt_outwrong,
                                    popsize = npop, sampsize = nrow(samp))

aipw_wrong <- combine.aipw(est=aipw_est_allwrong, boottbl=aipw_bt_allwrong,
                           popsize=npop, sampsize=nrow(samp))
aipw_wrong_noz <- combine.aipw(est = aipw_est_allwrong_noz, boottbl = aipw_bt_allwrong_noz,
                               popsize = npop,sampsize=nrow(samp))
print(paste0("I is ",i))

# }
res <- list(patetrue = pate_true, 
            aipw_unwtd_correcttrt =aipw_unwtd_correcttrt,
            aipw_unwtd_noz =aipw_unwtd_noz,
            aipw_wtd_wtdtrt =aipw_wtd_wtdtrt,
            aipw_wtd_wtdtrt_noz = aipw_wtd_wtdtrt_noz,
            aipw_wtd_mistrt = aipw_wtd_mistrt,
            aipw_wtd_mistrt_noz =aipw_wtd_mistrt_noz,
            aipw_unwtd_wtdtrt =aipw_unwtd_wtdtrt,
            aipw_unwtd_wtdtrt_noz= aipw_unwtd_wtdtrt_noz,
            aipw_allwrong = aipw_wrong,
            aipw_allwrong_noz = aipw_wrong_noz,
            aipw_correct_bench =aipw_correct_bench,
            aipw_ptrtwrong_bench = aipw_ptrtwrong_bench,
            aipw_outwrong_bench = aipw_outwrong_bench,
            time_aipw = time_aipw
            # time_aipw_noz = time_aipw_noz,
            # time_aipw_benchmark = time_aipw_benchmark
)
# saveRDS(res, file=paste0("results/results_tingting_", type_outcome, type_overlap, ".RDS"))
saveRDS(res, file = paste0("results-raw-2stage-bin/aipw/results_tingting_aipw", type_outcome, type_overlap,i,".RDS"))
