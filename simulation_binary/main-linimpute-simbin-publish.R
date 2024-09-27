library(nlme)
library(dplyr)
library(survey)
library(ggplot2)
library(sampling)
library(cmdstanr)
library("rootSolve")
library(LaplacesDemon)

type_y <- "bin"
type_outcome = "gb"
type_overlap = "notshared"

pathtocode = "./CODE-publish/"
source(paste0(pathtocode,"helpers.R"))
source(paste0(pathtocode, "calculate_pencomp_nospline_bin_sim.R"))
source(paste0(pathtocode, "calculate_wtd.R"))
source(paste0("STEP0-poptingting.R"))
mod <- cmdstan_model(paste0(pathtocode, "outcome_mod_logistic.stan"))

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
numboot=20
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
  pop <- pop %>% mutate(
    mu1 = bA + bx1*X1 + bAx1*X1 + bx2*X2 + bz*Z+ bAz*Z,
    mu0 = bx1*X1 + bx2*X2 + bz*Z
  )
  
  print(type_outcome)
}else if(type_outcome == "nogb"){
  pop <- pop %>% mutate(
    mu1 = bA + bx1*X1 + bAx1*X1 + bx2*X2 + bz*Z,
    mu0 = bx1*X1 + bx2*X2 + bz*Z
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
outcome_mod_nospline <- as.formula("Ysamp ~ X1 + X2 + Z")
outcome_mod_generalize <- formulaF(varList=varnames_outcome_generalize, y.name=y_varname)

## PARAMETERS FOR AIPTW
covariateXnames <- c("X1","X2", "Z")
propensity_mod_aiptw <- propensity_mod
outcome_mod_aiptw <- outcome_mod_nospline
ncarlo_aiptw <- 2000


# true PATE
pate_true <- mean(pop$Y1 - pop$Y0)

### ESTIMATE CONTAINERS
resnames <- c("est","withinvar","btwnvar","mivar","ci_tlength","ci_emplength","cicov_t","cicov_emp") # must match names from rubin function
wtdest_nospline <- vector(mode="list", length=8)
twostage_nospline <- vector(mode="list", length=8)
ate_nospline <- vector(mode="list", length=8)
names(wtdest_nospline) <- resnames
names(twostage_nospline) <- resnames
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
wtdest_nospline_boot <- matrix(NA,nrow=numboot, ncol=2)
ate_nospline_boot <- matrix(NA,nrow=numboot, ncol=2)
twostage_nospline_boot <- matrix(NA,nrow=numboot, ncol=2)

# BOOTSTRAP FOR FIRST STAGE IMPUTATION
boot <- 1
  while(boot <= numboot){ ## START bootstrap loop
    # (1) Bootstrap the sample (stratified by treatment group)
    bootflag <- 0
    bootind_cnt <- sample(ind_cnt, size=size_cnt, replace=T)
    bootind_trt <- sample(ind_trt, size=size_trt, replace=T)
    boot_all <- samp[c(bootind_cnt,bootind_trt),]
    boot_all$wts <- npop * normalize(boot_all$wts)
    
    imputedlist_nospline <- calculate.linimpute(outc_model = outcome_mod,
                                                trt.varname = "A", popdat = pop,
                                                y.varname = y_varname)
    
    ## Store results 
    ate_nospline_boot[boot,]=c(mean(imputedlist_nospline[[2]]-imputedlist_nospline[[1]]), 
                               var(imputedlist_nospline[[2]]-imputedlist_nospline[[1]])/nrow(samp))
    wtdres_nospline <- calculate.wtd(imputedlist = imputedlist_nospline, bootdf=samp)
    wtdest_nospline_boot[boot,] <- c(wtdres_nospline[[1]],
                                     wtdres_nospline[[2]])

    ## how to generalize with weights? If the sample is not a probability sample from population, then we would have to estimate them here.
    twostage_nospline_boot[boot,] <- c(imputedlist_nospline$pate_gen, imputedlist_nospline$varpategen)
    
    ## collect coefficient estimates
    print(paste0("Boot is ",boot))
    boot <- boot + 1
    bootflag <- 0
    # print(boot)
  } ## END bootstrap
  
  # step 3 rubin's combining rules and store results====
  ate_nospline <- Map(c, ate_nospline, combine.rubin(ate_nospline_boot))
  wtdest_nospline <- Map(c, wtdest_nospline, combine.rubin(wtdest_nospline_boot))
  twostage_nospline<- Map(c, twostage_nospline, combine.rubin(twostage_nospline_boot))
  
  # aiptw <- Map(c, aiptw, combine.aiptw(aiptw_boot))
  print(paste0("I is ",i))
  
# }
res <- list(patetrue = pate_true, 
            ate_nospline = ate_nospline, 
            wtdest_nospline = wtdest_nospline,
            twostage_nospline = twostage_nospline,
            bootflagct = bootflagct
            )
saveRDS(res, file = paste0("../results-raw/results_tingting_linimpute", type_outcome, type_overlap,i,".RDS"))
