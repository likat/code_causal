## Read in data
library(dplyr)
library(ggplot2)
library(survey)
library(cmdstanr)
require(LaplacesDemon)
require(polyapost)
require(parallel)
casenum <- Sys.getenv("SLURM_ARRAY_TASK_ID")
casenum <- as.numeric(casenum)
set.seed(1)
pathtocode <- "CODE/"
source(paste0(pathtocode, "helpers.R"))
dat <- readRDS("dat_use_faps.RDS")

## Obtain coefficients 
wtdes <- svydesign(ids=~tspsu,strata=~tsstrata,weights=~wts,data=dat)
outcome_coef <- svyglm(adltfs~ pctpov + age + wicelig + hhsize_scaled + educ + race, design=wtdes, family = "binomial")
trtment_coef <- svyglm(fsbenefit_12mo~pctpov + age + wicelig + race + educ + hhsize_scaled + 
         age:wicelig,family = "binomial", design=wtdes)
outcome_u <- c(-1, -0.6, -0.3)
trtment_u <- c(-1.5, -1, -0.6) # try 0.5, 1, 1.5
case_df <- data.frame(
  outcome = rep(outcome_u, 3),
  trtment = rep(trtment_u, each = 3)
)
case_df <- rbind(case_df, c(0,0))
case_df$name <- paste0(abs(case_df$outcome)*10, "_",abs(case_df$trtment)*10)


## Assign case number
case <- case_df[casenum,]

## Simulate U
dat$U <- rbinom(n=nrow(dat), size = 1,prob=0.5)

## DGP for Y
dat$mu <- expit(as.numeric(predict(outcome_coef,newdata=dat)) +case$outcome*dat$U)
dat$Y <- rbinom(nrow(dat), size=1,prob=dat$mu)

## DGP for A
dat$ptrt <- expit(predict(trtment_coef, newdata = dat) + case$trtment*dat$U)
dat$A <- rbinom(nrow(dat), size=1, prob=dat$ptrt)
wtdes <- svydesign(ids=~tspsu,strata=~tsstrata,weights=~wts,data=dat)

tempby <- svyby(~Y, by=~A, design = wtdes, FUN=svymean)

tempby$Y[2] - tempby$Y[1]


#### Start analysis #####
## load libraries
pathtocode <- "../causal-application-noaiptw/"
source(paste0(pathtocode, "helpers.R"))
source(paste0(pathtocode, "imputeF_application.R"))
# classical pencomp
source(paste0(pathtocode,"calculate_pencomp_application.R"))
# PENCOMP-WFPBB no spline in ptrt
source(paste0(pathtocode,"calculate_pwfpbb_wtdestonly_application_mclapply.R"))
# PENCOMP-WFPBB with spline in ptrt
source(paste0(pathtocode,"calculate_pwfpbb_ptrtsp_application.R"))
source(paste0(pathtocode,"calculate_ptrtspl_application.R"))
# linear impute
source(paste0(pathtocode, "calculate_pencomp_nospline_application.R"))
## calculating weighted estimates
source(paste0(pathtocode,"calculate_wtd.R"))
## compile stan file
mod <- cmdstan_model(paste0(pathtocode, "outcome_mod_logistic.stan"))

## AIPTW
pathtocode <- "../causal-bin-sim-new/CODE/"
source(paste0(pathtocode, "run_AIPTW_new.R"))

# ## Read in data
# dat <- readRDS("data/faps_cleaned.RDS") 
dat$tsstrata <- dplyr::dense_rank(dat$tsstrata)
dat$tspsu <- dplyr::dense_rank(dat$tspsu)
nstrata <- dat$tsstrata %>% unique() %>% length()
cluster_stratum_ref <- dat %>% group_by(tsstrata, tspsu) %>% summarize(cts=n()) 
cluster_stratum_ref <- cluster_stratum_ref %>% group_by(tsstrata) %>% mutate(nclus = length(unique(tspsu)),
                                                                             mstar = nclus / (nclus-1)) %>% 
  select(-cts)
ptrtvars <- readRDS("ptrtmodelvars_bic.RDS")
# ptrtmod <- "fsbenefit_12mo ~ pctpov + age + race + educ + hhsize_scaled" %>% as.formula
ptrtmod <- ptrtvars %>% as.formula
ptrtmod <- "A ~ pctpov + age + wicelig + race + educ + hhsize_scaled + 
    age:wicelig" %>% as.formula
ptrtmod_aiptw <- "fsbenefit_12mo ~ pctpov + age + wicelig + race + educ + as.vector(hhsize_scaled) + age:wicelig"

outcome_mod <- "Y ~ pctpov + age + wicelig + hhsize_scaled + educ + race" %>% as.formula
# outcome_mod <- readRDS("outmodeltxt_faps.RDS") %>% as.formula
outcome_mod_linimpute = outcome_mod
# outcome_mod <- "adltfs ~ pctpov + age + educ + hhsize_scaled + race" %>% as.formula
# outcome_mod_linimpute <- "adltfs ~ pctpov + age + educ + hhsize_scaled + race" %>% as.formula
## parameters
treatment_name <- "A"
y_varname <- "Y"
wfpbb_fdraw <- 20
Tfact <- 20
numboot = 200
## Begin analysis
#- Methods we are looking at: 
# classical PENCOMP
# WFPBB-PENCOMP with $Z$ variables only (wtdest)
# WFPBB-PENCOMP using the spline for sampling weights (wtdest)
resboot_aipw <- data.frame(
  est = rep(NA, numboot),
  est_wtd = rep(NA, numboot)
)
res_aipw <- data.frame(
  est = rep(NA, numboot),
  est_wtd = rep(NA, numboot),
  time = rep(NA, numboot)
)
resnames <- c("est","withinvar","btwnvar","mivar","ci_tlength","ci_emplength","cicov_t","cicov_emp") # must match names from rubin function
wtdest_pencomp<- vector(mode="list", length=8)
wtdest_linimpute <- vector(mode="list", length=8)
wtdest_ptrtzonly <- vector(mode="list", length=8)
wtdest_ptrtspl <- vector(mode="list", length=8)
iptwest_ptrtspl <- vector(mode="list", length=8)
ate_pencomp<- vector(mode="list", length=8)
ate_ptrtzonly <- vector(mode="list", length=8)
ate_ptrtspl <- vector(mode="list", length=8)
ate_linimpute <- vector(mode="list", length=8)

coef_pencomp_boot <- vector(mode="list", length = 4)
coef_ptrtspl_boot<- vector(mode="list", length = 4)
coef_ptrtzonly_boot <- vector(mode="list", length = 4)

ate_linimpute_boot <- matrix(NA,nrow=numboot, ncol=2)
ate_pencomp_boot <- matrix(NA,nrow=numboot, ncol=2)
ate_ptrtzonly_boot  <- matrix(NA,nrow=numboot, ncol=2)
ate_ptrtspl_boot  <- matrix(NA,nrow=numboot, ncol=2)
wtdest_linimpute_boot <- matrix(NA,nrow=numboot, ncol=2)
wtdest_pencomp_boot <- matrix(NA,nrow=numboot, ncol=2)
wtdest_ptrtzonly_boot  <- matrix(NA,nrow=numboot, ncol=2)
wtdest_ptrtspl_boot  <- matrix(NA,nrow=numboot, ncol=2)
iptwest_ptrtspl_boot <- matrix(NA,nrow=numboot, ncol=2)
pate_synthpop_ptrtzonly_boot <- rep(NA, nrow = numboot)
names(wtdest_linimpute) <- resnames
names(wtdest_pencomp) <- resnames
names(wtdest_ptrtzonly) <- resnames
names(wtdest_ptrtspl) <- resnames
names(iptwest_ptrtspl) <- resnames
names(ate_linimpute) <- resnames
names(ate_pencomp) <- resnames
names(ate_ptrtzonly) <- resnames
names(ate_ptrtspl) <- resnames
bootflag <- 0

## Setup for sample
size_cnt <- sum(dat$fsbenefit_12mo==0)
size_trt <- sum(dat$fsbenefit_12mo==1)
ind_cnt <- which(dat$fsbenefit_12mo==0)
ind_trt <- which(dat$fsbenefit_12mo==1)
npop <- sum(dat$wts)
phat_boot <- list(
  pencomp = matrix(nrow = numboot,ncol=nrow(dat)),
  ptrtzonly = matrix(nrow=numboot, ncol = nrow(dat)),
  ptrtspl = matrix(nrow=numboot, ncol=nrow(dat))
)

## Begin bootstrap
boot <- 1
bootflagct = 0
aipwest <- AIPW.application(dataInput = dat,
                             propen.mod = ptrtmod_aiptw,
                             out.name = y_varname,
                             trt.name = treatment_name,
                             modely1 = outcome_mod,
                             modely0 = outcome_mod,
                             outcomefamily = "binomial")
res_aipw$est <- aipwest$pate_unwtd
res_aipw$est_wtd <- aipwest$pate_wtd
res_aipw$time <- aipwest$time

while(boot <= numboot){ ## START bootstrap loop
  # (1) Bootstrap the sample (stratified by treatment group)
  set.seed(boot)
  
  bootflag <- 0
  cluster_stratum_ref$resamp <- 0
  cluster_stratum_ref$bbresamp <- 0
  cluster_stratum_ref$clusmultiplier <- 0
  cluster_stratum_ref$bbclusmultiplier <- 0
  
  dat$clusmultiplier <- NULL
  dat$resamp <- NULL
  dat$bbclusmultiplier <- NULL
  dat$bbresamp <- NULL
  ## with clusters and strata, resample clusters
  bootind_cnt <- sample(ind_cnt, size=size_cnt, replace=T)
  bootind_trt <- sample(ind_trt, size=size_trt, replace=T)
  boot_all <- dat[c(bootind_cnt,bootind_trt),]
  
  boot_all$wts <-  nrow(boot_all) *normalize(boot_all$wts)
  boot_all$wts_scaled <- nrow(boot_all) * normalize(boot_all$wts)
  
  for(k in 1:nstrata){
    cluspool <- cluster_stratum_ref$tspsu[cluster_stratum_ref$tsstrata==k]
    tempnumclus <- length(cluspool)
    bootstrap <- BayesianBootstrap(cluspool, 1)
    resamp <- rmultinom(1, tempnumclus - 1, bootstrap)
    cluster_stratum_ref$bbresamp[cluster_stratum_ref$tspsu %in% cluspool] <- resamp
  }
  cluster_stratum_ref$clusmultiplier <- cluster_stratum_ref$resamp*cluster_stratum_ref$mstar
  cluster_stratum_ref$bbclusmultiplier <- cluster_stratum_ref$bbresamp*cluster_stratum_ref$mstar
  
  dat <- dat %>% left_join(cluster_stratum_ref %>% select(tspsu, tsstrata,clusmultiplier,bbclusmultiplier,
                                                          resamp, bbresamp))
  boot_allbb <- dat[dat$bbresamp != 0,]
  
  # # step 1 imputation ===
  T_pop = 20
  bootdf = boot_allbb
  samp = dat
  propen_model=ptrtmod
  bootdf = boot_allbb
  outc_model = outcome_mod
  y.varname = y_varname
  F_draw = 20
  trt.varname = treatment_name
  num.knot=15
  imputedlist_ptrtzonly <- calculate.pencomp.wfpbb(propen_model=ptrtmod,
                                                   bootdf = boot_allbb,
                                                   outc_model = outcome_mod,
                                                   y.varname = y_varname,
                                                   F_draw = 20,
                                                   trt.varname = treatment_name,
                                                   num.knot=15)
  
  if(sum(is.na( imputedlist_ptrtzonly))>0){bootflagct = bootflagct + 1; next}
  if(bootflag == 1){bootflagct = bootflagct + 1; next}
  
  
  pate_synthpop_ptrtzonly_boot[boot] <- imputedlist_ptrtzonly$ate_pop

  temp <- AIPW.application(dataInput = boot_all,
                           propen.mod = ptrtmod_aiptw,
                           out.name = y_varname,
                           trt.name = treatment_name,
                           modely1 = outcome_mod,
                           modely0 = outcome_mod,
                           outcomefamily = "binomial")
  
  resboot_aipw$est[boot] <- temp$pate_unwtd
  resboot_aipw$est_wtd[boot] <- temp$pate_wtd
  # res_aipw$time[boot] <- temp$time
  
  print(paste0("Boot is ",boot))
  boot <- boot + 1
  bootflag <- 0
} ## END bootstrap

### Combine
# ate_linimpute <- combine.rubin.applied(ate_linimpute_boot)
# ate_pencomp <-combine.rubin.applied(ate_pencomp_boot)
# ate_ptrtzonly <-combine.rubin.applied(ate_ptrtzonly_boot)
# ate_ptrtspl <- combine.rubin.applied(ate_ptrtspl_boot)
# wtdest_linimpute <- combine.rubin.applied(wtdest_linimpute_boot)
# wtdest_pencomp <-combine.rubin.applied(wtdest_pencomp_boot)
# wtdest_ptrtzonly <-combine.rubin.applied(wtdest_ptrtzonly_boot)
# wtdest_ptrtspl <- combine.rubin.applied(wtdest_ptrtspl_boot)
# iptwest_ptrtspl <- combine.rubin.applied(iptwest_ptrtspl_boot)
# phat <- Map(colMeans, phat_boot)
res <- list(
  pate_synthpop_ptrtzonly_boot=pate_synthpop_ptrtzonly_boot,
  aipw_est = res_aipw,
  aipw_boot = resboot_aipw,
  case = case$name,
  outcomecoef = case$outcome,
  trtcoef = case$trtment,
  bootflagct = bootflagct
)
saveRDS(res,paste0("results/step3_simulation_", case$name,"_new.RDS"))
