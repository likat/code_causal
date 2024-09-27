## load libraries
library(dplyr)
library(ggplot2)
library(survey)
library(cmdstanr)
library(rootSolve)
#################################################################################
## Read in code
dat = readRDS("data/dat_use_faps.RDS")
pathtocode = "./CODE-publish/"
## AIPTW
source(paste0(pathtocode, "gcomputeFunc_application.R"))
source(paste0(pathtocode, "run_AIPTW_application.R"))
source(paste0(pathtocode, "helpers.R"))
#################################################################################

## Read in data
ptrtmod_aiptw <- "fsbenefit_12mo ~ pctpov + age + wicelig + race + educ + as.vector(hhsize_scaled) + age:wicelig"
outcome_mod_aiptw <- "adltfs ~ pctpov + age + wicelig + educ +as.vector(hhsize_scaled) + race" %>% as.formula
## parameters
treatment_name <- "fsbenefit_12mo"
y_varname <- "adltfs"
numboot = 200
ncarlo_aiptw = 4000
covariateXnames <- c("pctpov" , "age" , "wicelig" , "race" , "educ" , "hhsize_scaled")
cluster_stratum_ref <- dat %>% group_by(tsstrata, tspsu) %>% summarize(cts=n()) 
cluster_stratum_ref <- cluster_stratum_ref %>% group_by(tsstrata) %>% mutate(nclus = length(unique(tspsu)),
                                                                             mstar = nclus / (nclus-1))
nstrata <- dat$tsstrata %>% unique() %>% length()

## Begin analysis
aiptw_boot <- NA

## Setup for sample
size_cnt <- sum(dat$fsbenefit_12mo==0)
size_trt <- sum(dat$fsbenefit_12mo==1)
ind_cnt <- which(dat$fsbenefit_12mo==0)
ind_trt <- which(dat$fsbenefit_12mo==1)

## Begin bootstrap
# boot <- 1
bootflagct = 0


# while(boot <= numboot){ ## START bootstrap loop
boot <- Sys.getenv("SLURM_ARRAY_TASK_ID")
boot <- as.numeric(boot)
time0 <- NA
time1 <- NA
set.seed(boot)
  # (1) Bootstrap the sample (stratified by treatment group)
  bootflag <- 0
  cluster_stratum_ref$resamp <- 0
  cluster_stratum_ref$clusmultiplier <- 0
  dat$clusmultiplier <- NULL
  dat$resamp <- NULL
  ## with clusters and strata, resample clusters
  for(h in 1:nstrata){
    cluspool <- cluster_stratum_ref$tspsu[cluster_stratum_ref$tsstrata==h]
    tempnumclus <- length(cluspool)
    resamp <- sort(sample(cluspool,tempnumclus-1,replace=T)) %>% unique()
    cluster_stratum_ref$resamp[resamp] <- table(resamp)
  }
  
  cluster_stratum_ref$clusmultiplier <- cluster_stratum_ref$resamp*cluster_stratum_ref$mstar
  dat <- dat %>% left_join(cluster_stratum_ref %>% dplyr::select(tspsu, tsstrata,clusmultiplier,resamp))
  boot_all <- dat[rep(1:nrow(dat), dat$resamp),]
  boot_all$wts <- sum(dat$wts) * normalize(boot_all$wts)
  time0 <- Sys.time()
  aiptw_boot  <- AIPTW(dataInput=boot_all,
                         covariateXnames= covariateXnames,
                         sampwts = boot_all$wts,
                         treat.varname= treatment_name,
                         outcome.varname= y_varname,
                         propen.model2a =ptrtmod_aiptw,
                         modely1 = outcome_mod_aiptw,
                         modely0 = outcome_mod_aiptw,
                         numCarlo = ncarlo_aiptw)
  time1 <- Sys.time()
  print(paste0("Boot is ",boot))
  bootflag <- 0
# } ## END bootstrap

### Combine
res <- list(est = aiptw_boot, time = time1-time0)
  saveRDS(res,paste0("results/step1_res_faps_aiptwbootraw", boot,".RDS"))
  
