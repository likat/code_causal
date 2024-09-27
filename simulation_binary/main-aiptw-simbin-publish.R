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

source(paste0(pathtocode, "gcomputeFunc_bin.R"))
source(paste0(pathtocode, "run_AIPTW_bin.R"))
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
outcome_mod_aiptw <- outcome_mod
ncarlo_aiptw <- 2000


# true PATE
pate_true <- mean(pop$Y1 - pop$Y0)

### ESTIMATE CONTAINERS
sate <- rep(NA, sim)
resnames <- c("est","withinvar","btwnvar","mivar","ci_tlength","ci_emplength","cicov_t","cicov_emp") # must match names from rubin function
resnames_aiptw <- c("est", "var", "cilower","ciupper")
time_aiptw <- 0

aiptw <- vector(mode = "list", length = 4)
names(aiptw) <- resnames_aiptw

bootflag <- 0
bootflagct <- 0
## To replicate a true experiment, recalculate the outcome every sample and treatment assignment
pop$popID <- 1:nrow(pop)
numboot = 20



### BEGIN SIM
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


### bootstrap containers
aiptw_boot <- rep(NA, numboot)

# BOOTSTRAP FOR FIRST STAGE IMPUTATION
boot <- 1
while(boot <= numboot){ ## START bootstrap loop
  # (1) Bootstrap the sample (stratified by treatment group)
  bootflag <- 0
  
  bootind_cnt <- sample(ind_cnt, size=size_cnt, replace=T)
  bootind_trt <- sample(ind_trt, size=size_trt, replace=T)
  boot_all <- samp[c(bootind_cnt,bootind_trt),]
  boot_all$wts <- npop * normalize(boot_all$wts)
  
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
  # time_aiptw <- time_aiptw + aiptw_bootres$time/numboot

  aiptw_boot[boot] <- aiptw_bootres
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
saveRDS(res, file = paste0("../results-raw/results_tingting_aiptw", type_outcome, type_overlap,i,".RDS"))
