library(nlme)
library(dplyr)
library(survey)
library(ggplot2)
library(sampling)
library(cmdstanr)
library("rootSolve")
library(LaplacesDemon)

type_y <- "cts"
type_outcome = "gb"
type_overlap = "notshared"

pathtocode = "./CODE-publish/"
source(paste0(pathtocode,"helpers.R"))
source(paste0(pathtocode, "calculate_pencomp_bin_sim.R"))
source(paste0(pathtocode, "calculate_wfpbb_pencomp_bin_sim.R"))
source(paste0(pathtocode, "imputeF_bin_sim.R"))
source(paste0(pathtocode, "calculate_pencomp_nospline_bin_sim.R"))
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
numboot=20
set.seed(42)

## OUTCOME GENERATION
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
# bAx1 <- 0.8
# bAx2 <- 0.5
bAx1 <- 0.8
bAx2 <- 0.5

bz <- -0.2
bAz <- 0.4
# bax1z <- 0
# bax2z <- 0
#***************************************************************
if(type_outcome == "gb"){
  # pop <- pop %>% mutate(
  #   mu1 = bA + bx1*X1 + bAx1*X1 + bx2*X2 + bz*Z+ bAz*Z,
  #   mu0 = bx1*X1 + bx2*X2 + bz*Z
  #   )
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
  # pop$Y1 <- rnorm(n=npop, pop$mu1)
  # pop$Y0 <- rnorm(n=npop, pop$mu0)
pop$Y1 <- rbinom(n = nrow(pop), size = 1, prob = expit(pop$mu1))
pop$Y0 <- rbinom(n = nrow(pop), size = 1, prob = expit(pop$mu0))

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

coef_pencomp <- vector(mode = "list", length = 6)
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

sampsize <- rep(NA, sim)
for(j in 1:sim){
  set.seed(j)
  pop$I <- NULL
  pop$I <- rbinom(n=npop, size=1, prob = pop$pselect)
  sampsize[j] <- sum(pop$I)
}
saveRDS(sampsize,paste0("../results-clean/sampsize_", type_outcome, type_overlap, ".RDS"))
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
  
mean(samp$Ysamp)
mean(samp$Ysamp[samp$A==1])-mean(samp$Ysamp[samp$A==0])
pate_true
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
    imputedlist_nospline <- calculate.linimpute(outc_model = outcome_mod_nospline,
                                                trt.varname = "A", popdat = pop,
                                                y.varname = y_varname)
    # if(bootflag == 1){bootflagct = bootflagct + 1; next}
    imputedlist_spline   <- calculate.pencomp(propen_model=propensity_mod,
                                              outc_model = outcome_mod,
                                              y.varname = y_varname,
                                              trt.varname = "A",
                                              popdat = pop,
                                              num.knot=10)
    if(sum(is.na(imputedlist_spline)) > 0){bootflagct = bootflagct + 1; next}
    imputedlist_wfpbbspline <- calculate.wfpbb.pencomp(propen_model=propensity_mod,
                                                       outc_model = outcome_mod,
                                                       y.varname = y_varname,
                                                       bootdf = boot_allbb,
                                                       uniqueind = uniqueind_bb,
                                                       polyawts = wts_bb,
                                                       trt.varname = "A",
                                                       num.knot = 10,
                                                       F_draw=1)
    if(sum(c(is.na(imputedlist_wfpbbspline)))>0){bootflagct = bootflagct + 1; next}
    # aiptw_bootres <- AIPTW(dataInput=samp, 
    #                        covariateXnames= covariateXnames,
    #                        sampwts = samp$wts,
    #                        treat.varname= "A", 
    #                        outcome.varname= y_varname, 
    #                        propen.model2a = propensity_mod_aiptw,
    #                        modely1 = outcome_mod_aiptw,
    #                        modely0 = outcome_mod_aiptw,
    #                        numCarlo=ncarlo_aiptw)
    # if(sum(c(is.na(imputedlist_spline), is.na(imputedlist_wfpbbspline)))>0){print("error!")}
    #if(bootflag == 1){bootflagct = bootflagct + 1; next}
    
    ## Store results 
    # time_aiptw <- time_aiptw + aiptw_bootres$time/numboot
    time_pencompwfpbb <- time_pencompwfpbb + imputedlist_wfpbbspline$time_pencompwfpbb/numboot
    
    # aiptw_boot[boot] <- aiptw_bootres$result
    synthpop_pencompwfpbb_boot[boot,] <- c(imputedlist_wfpbbspline$pate_gen_synthpop, imputedlist_wfpbbspline$varpategen_synthpop)
    #   # Calculate ATE
    ate_pencomp_boot[boot,]=c(mean(imputedlist_spline[[2]]-imputedlist_spline[[1]]),
                              var(imputedlist_spline[[2]]-imputedlist_spline[[1]])/nrow(samp))
    ate_nospline_boot[boot,]=c(mean(imputedlist_nospline[[2]]-imputedlist_nospline[[1]]), 
                               var(imputedlist_nospline[[2]]-imputedlist_nospline[[1]])/nrow(samp))
    ate_pencompwfpbb_boot[boot,]=c(mean(imputedlist_wfpbbspline[[2]]-imputedlist_wfpbbspline[[1]]),
                              var(imputedlist_wfpbbspline[[2]]-imputedlist_wfpbbspline[[1]])/nrow(samp))
    
    # print(ate_nospline_boot[boot,])
    # step 2 generalize! ====
    wtdres <- calculate.wtd(imputedlist = imputedlist_spline, bootdf=samp)
    wtdres_nospline <- calculate.wtd(imputedlist = imputedlist_nospline, bootdf=samp)
    wtdres_wfpbb <- imputedlist_wfpbbspline$wtdest
    wtdest_pencomp_boot[boot,] <- c(wtdres[[1]],
                                    wtdres[[2]])
    wtdest_pencompwfpbb_boot[boot,] <- wtdres_wfpbb
    wtdest_nospline_boot[boot,] <- c(wtdres_nospline[[1]],
                                     wtdres_nospline[[2]])

    ## how to generalize with weights? If the sample is not a probability sample from population, then we would have to estimate them here.
    twostage_pencomp_boot[boot,] <- c(imputedlist_spline$pate_gen, imputedlist_spline$varpategen)
    twostage_pencompwfpbb_boot[boot,] <- c(imputedlist_wfpbbspline$pate_gen, imputedlist_wfpbbspline$varpategen)
    twostage_nospline_boot[boot,] <- c(imputedlist_nospline$pate_gen, imputedlist_nospline$varpategen)
    
    ## collect coefficient estimates
    coef_pencomp_boot <- Map(rbind, coef_pencomp_boot, imputedlist_spline$coeflistres)
    coef_pencompwfpbb_boot <- Map(rbind, coef_pencompwfpbb_boot, imputedlist_wfpbbspline$coeflistres)
    print(paste0("Boot is ",boot))
    boot <- boot + 1
    bootflag <- 0
     print(paste0("Boot is ", boot))
  } ## END bootstrap
  
  # step 3 rubin's combining rules and store results====
  ate_pencomp <- Map(c,ate_pencomp, combine.rubin(ate_pencomp_boot))
  ate_pencomp_wfpbb <- Map(c, ate_pencomp_wfpbb, combine.rubin(ate_pencompwfpbb_boot))
  ate_nospline <- Map(c, ate_nospline, combine.rubin(ate_nospline_boot))
  wtdest_pencomp <- Map(c, wtdest_pencomp, combine.rubin(wtdest_pencomp_boot))
  wtdest_pencompwfpbb <- Map(c, wtdest_pencompwfpbb, combine.rubin(wtdest_pencompwfpbb_boot))
  wtdest_nospline <- Map(c, wtdest_nospline, combine.rubin(wtdest_nospline_boot))
  twostage_pencomp <- Map(c, twostage_pencomp, combine.rubin(twostage_pencomp_boot))
  twostage_pencompwfpbb <- Map(c, twostage_pencompwfpbb, combine.rubin(twostage_pencompwfpbb_boot))
  twostage_nospline<- Map(c, twostage_nospline, combine.rubin(twostage_nospline_boot))
  
  # aiptw <- Map(c, aiptw, combine.aiptw(aiptw_boot))
  synthpop_pencompwfpbb <- Map(c, synthpop_pencompwfpbb, combine.rubin(synthpop_pencompwfpbb_boot))
  coef_pencomp <- Map(rbind, coef_pencomp,lapply(coef_pencomp_boot, colMeans))
  coef_pencompwfpbb <- Map(rbind, coef_pencompwfpbb, lapply(coef_pencompwfpbb_boot, colMeans, na.rm=T) )
  print(paste0("I is ",i))
  
# }
res <- list(patetrue = pate_true, sate = sate,
            ate_pencomp = ate_pencomp, 
            ate_nospline = ate_nospline, 
            ate_pencomp_wfpbb = ate_pencomp_wfpbb,
            wtdest_pencomp = wtdest_pencomp, 
            wtdest_pencompwfpbb = wtdest_pencompwfpbb,
            wtdest_nospline = wtdest_nospline,
            twostage_pencomp = twostage_pencomp,
            twostage_pencompwfpbb = twostage_pencompwfpbb,
            twostage_nospline = twostage_nospline,
            synthpop_pencompwfpbb = synthpop_pencompwfpbb,
            coef_pencomp = coef_pencomp,
            coef_pencompwfpbb = coef_pencompwfpbb,
            time_pencompwfpbb = time_pencompwfpbb,
            # time_aiptw = time_aiptw,
            # aiptw = aiptw,
            bootflagct = bootflagct
            )
# saveRDS(res, file=paste0("results/results_tingting_", type_outcome, type_overlap, ".RDS"))
saveRDS(res, file = paste0("../results-raw/results_tingting_noaiptw", type_outcome, type_overlap,i,".RDS"))
