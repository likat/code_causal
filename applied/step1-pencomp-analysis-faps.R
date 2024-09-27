## load libraries
library(dplyr)
library(ggplot2)
library(survey)
library(cmdstanr)
library(LaplacesDemon)

#################################################################################
## Read in code
pathtocode <- "./CODE-publish/"
source(paste0(pathtocode, "helpers.R"))
source(paste0(pathtocode, "imputeF_application.R"))
# classical pencomp
source(paste0(pathtocode,"calculate_pencomp_application.R"))
# PENCOMP-WFPBB no spline in ptrt
source(paste0(pathtocode,"calculate_pwfpbb_wtdestonly_application.R"))
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

#################################################################################

## Read in data
dat <- readRDS("data/faps_cleaned.RDS") 
dat$tsstrata <- dplyr::dense_rank(dat$tsstrata)
dat$tspsu <- dplyr::dense_rank(dat$tspsu)

## Arrange
nstrata <- dat$tsstrata %>% unique() %>% length()
cluster_stratum_ref <- dat %>% group_by(tsstrata, tspsu) %>% summarize(cts=n()) 
cluster_stratum_ref <- cluster_stratum_ref %>% group_by(tsstrata) %>% mutate(nclus = length(unique(tspsu)),
                                                                             mstar = nclus / (nclus-1)) 
ptrtvars <- readRDS("ptrtmodelvars_bic.RDS")
# ptrtmod <- "fsbenefit_12mo ~ pctpov + age + race + educ + hhsize_scaled" %>% as.formula
ptrtmod <- ptrtvars %>% as.formula
ptrtmod_aiptw <- "fsbenefit_12mo ~ pctpov + age + wicelig + race + educ + as.vector(hhsize_scaled) + age:wicelig"
# outcome_mod <- readRDS("outmodeltxt_faps.RDS") %>% as.formula
outcome_mod <- "adltfs ~ pctpov + age + wicelig + hhsize_scaled + educ + race" %>% as.formula
names_fixedcoef <- model.matrix(outcome_mod, dat) %>% colnames

outcome_mod_linimpute = outcome_mod
outcome_mod_aiptw <- "adltfs ~ pctpov + age + wicelig + race + educ + 
    age:wicelig" %>% as.formula
# outcome_mod_linimpute <- "adltfs ~ pctpov + age + educ + hhsize_scaled + race" %>% as.formula
## parameters
treatment_name <- "fsbenefit_12mo"
y_varname <- "adltfs"
wfpbb_fdraw <- 10
Tfact <- 20
numboot = 200

## Begin analysis
#- Methods we are looking at: 
# classical PENCOMP
# WFPBB-PENCOMP with $Z$ variables only (wtdest)
# WFPBB-PENCOMP using the spline for sampling weights (wtdest)
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
pate_synthpop_ptrtzonly <- vector(mode="list", length = 8)
pate_synthpop_ptrtspl <- vector(mode= "list", length= 8)

coef_pencomp_boot <- vector(mode="list", length = 6)
coef_ptrtspl_boot<- vector(mode="list", length = 5)
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
pate_synthpop_ptrtzonly_boot <- rep(NA, numboot)
pate_synthpop_ptrtspl_boot <- rep(NA, numboot)

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


## The cutoffs are based on the model for treatment propensity w log(1/wts) spline
boot_all <- dat
phat_imputed <- calculate.ptrtspline(fitdat = dat,
                                     sampdat = dat,
                                     num.knot = 5,
                                     propen.model = ptrtmod,
                                     trtname=treatment_name)
dat$phat <- phat_imputed$ptrt_fitdat
# trt0_phatrange <- range(dat$phat[which(dat[,treatment_name]==0)])
# trt1_phatrange <- range(dat$phat[which(dat[,treatment_name]==1)])
### ALL OVERLAPPING VALUES
# overlap_samp <- c(max(trt0_phatrange[1], trt1_phatrange[1]),
#                   min(trt0_phatrange[2], trt1_phatrange[2]))
### MORE STRINGENT CUTOFF CRITERIA
overlap_samp <- c(-3, 3)
### the sample should only contain observations that have estimated phat in the overlapping region
dat <- dat %>% filter(between(phat, overlap_samp[1], overlap_samp[2]))
npop <- sum(dat$wts)
phat_boot <- list(
  pencomp = matrix(nrow = numboot,ncol=nrow(dat)),
  ptrtzonly = matrix(nrow=numboot, ncol = nrow(dat)),
  ptrtspl = matrix(nrow=numboot, ncol=nrow(dat))
)

## Check overlap ##############
#####  Save data as dat_use_faps if first time running
# pdf(file = "results/faps_overlap.pdf",width = 7,height=7)
# ggplot(dat, aes(x=phat, fill=factor(fsbenefit_12mo))) +
#   geom_histogram(alpha=0.2, position="identity") +
#   labs(x="Estimated propensity score", fill="Trt group")+
#   theme_minimal()
# dev.off()
# saveRDS(dat,"data/dat_use_faps.RDS")
###############################
## GET HT estimate
wtdes <- svydesign(ids=~tspsu,strata=~tsstrata,weights=~wts,data=dat)
est_ht <- svyby(~adltfs,by=~fsbenefit_12mo, design=wtdes, FUN=svymean)
est_ht <- est_ht$adltfs[2]-est_ht$adltfs[1]

## GET Naive regression estimate
mod0 <- svyglm(adltfs~pctpov + age + educ + hhsize_scaled + race + wicelig+age:wicelig+fsbenefit_12mo,
            des = wtdes, family = "binomial")
dat$temppred <- predict(mod0,newdata=dat, type="response")
wtdes2 <- svydesign(ids=~tspsu,strata=~tsstrata,weights=~wts,data=dat)
est_naivemod <- svyby(~temppred,by=~fsbenefit_12mo, design=wtdes2, FUN=svymean)
est_naivemod <- est_naivemod$temppred[2]-est_naivemod$temppred[1]

## Begin bootstrap
boot <- 1
bootflagct = 0
# for(boot in 1:numboot){ ## START bootstrap loop
while(boot <= numboot){
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
  for(h in 1:nstrata){
    cluspool <- cluster_stratum_ref$tspsu[cluster_stratum_ref$tsstrata==h]
    tempnumclus <- length(cluspool)
    resamp <- sample(cluspool,tempnumclus-1,replace=T)
    cluster_stratum_ref$resamp[cluster_stratum_ref$tspsu %in% unique(resamp)] <- table(resamp)
  }
  
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
  boot_all <- dat[dat$resamp !=0,]
  boot_allbb <- dat[dat$bbresamp != 0,]
  # # step 1 imputation ===
  imputedlist_linimpute <- calculate.linimpute(outc_model = outcome_mod_linimpute,
                                               y.varname = y_varname,
                                               trt.varname = treatment_name)
  imputedlist_pencomp  <- calculate.pencomp(propen_model=ptrtmod,
                                            outc_model = outcome_mod,
                                            y.varname = y_varname,
                                            trt.varname = treatment_name,
                                            num.knot=10)
  if(sum(c(is.na(imputedlist_pencomp)))>0){bootflagct = bootflagct + 1; next}
  imputedlist_ptrtzonly <- calculate.pencomp.wfpbb(propen_model=ptrtmod,
                                                   bootdf = boot_allbb,
                                                    outc_model = outcome_mod,
                                                    y.varname = y_varname,
                                                    F_draw = wfpbb_fdraw,
                                                    trt.varname = treatment_name,
                                                    num.knot=10)
  if(sum(is.na( imputedlist_ptrtzonly))>0){bootflagct = bootflagct + 1; next}
  imputedlist_ptrtspl <- calculate.pwfpbb.ptrtspl(propen_model=ptrtmod,
                                                  bootdf = boot_allbb,
                                                  outc_model = outcome_mod,
                                                  y.varname = y_varname,
                                                  F_draw = wfpbb_fdraw,
                                                  trt.varname = treatment_name,
                                                  num.knot=10)
  if(sum(is.na(imputedlist_ptrtspl))>0){bootflagct = bootflagct + 1; next}
  pate_synthpop_ptrtzonly_boot[boot] <- imputedlist_ptrtzonly$ate_pop$estimate
  pate_synthpop_ptrtspl_boot[boot] <- imputedlist_ptrtspl$ate_pop$estimate
  
  #if(bootflag == 1){bootflagct = bootflagct + 1; next}
  
  #   # Calculate ATE
  ate_linimpute_boot[boot,] = c(mean(imputedlist_linimpute[[2]]-imputedlist_linimpute[[1]]),
                                var(imputedlist_linimpute[[2]]-imputedlist_linimpute[[1]])/nrow(dat))
  ate_pencomp_boot[boot,]=c(mean(imputedlist_pencomp[[2]]-imputedlist_pencomp[[1]]),
                            var(imputedlist_pencomp[[2]]-imputedlist_pencomp[[1]])/nrow(dat))
  ate_ptrtzonly_boot[boot,]=c(imputedlist_ptrtzonly$ate_samp %>% unlist %>% as.numeric)
  ate_ptrtspl_boot[boot,]=c(imputedlist_ptrtspl$ate_samp %>% unlist %>% as.numeric)
  
  # step 2 generalize! ====
  wtdres_linimpute <- calculate.wtd(imputedlist = imputedlist_linimpute, bootdf=dat)
  wtdres_pencomp <- calculate.wtd(imputedlist = imputedlist_pencomp, bootdf=dat)
  wtdres_ptrtzonly <- imputedlist_ptrtzonly$wtdest
  wtdres_ptrtspl <- imputedlist_ptrtspl$wtdest
  
  wtdest_linimpute_boot[boot,] <- c(wtdres_linimpute[[1]],
                                    wtdres_linimpute[[2]])
  wtdest_pencomp_boot[boot,] <- c(wtdres_pencomp[[1]],
                                  wtdres_pencomp[[2]])
  wtdest_ptrtzonly_boot[boot,] <- c(wtdres_ptrtzonly[[1]],
                                    wtdres_ptrtzonly[[2]])
  wtdest_ptrtspl_boot[boot,] <- c(wtdres_ptrtspl[[1]],
                                  wtdres_ptrtspl[[2]])
  iptwest_ptrtspl_boot[boot,] <- c(imputedlist_ptrtspl$est_iptw,
                                   as.numeric(imputedlist_ptrtspl$var_iptw))
  
  coef_pencomp_boot <- Map(rbind,coef_pencomp_boot, imputedlist_pencomp$coeflistres)
  coef_ptrtspl_boot <- Map(rbind,coef_ptrtspl_boot, imputedlist_ptrtspl$coeflistres)
  coef_ptrtzonly_boot <- Map(rbind,coef_ptrtzonly_boot, imputedlist_ptrtzonly$coeflistres)
  phat_boot$pencomp[boot,] <- imputedlist_pencomp$phat_samp
  phat_boot$ptrtzonly[boot,] <- imputedlist_ptrtzonly$phat_samp
  phat_boot$ptrtspl[boot,] <- imputedlist_ptrtspl$phat_samp
  print(paste0("Boot is ",boot))
  boot <- boot + 1
  bootflag <- 0
} ## END bootstrap

### Combine
ate_linimpute <- combine.rubin.applied(ate_linimpute_boot)
ate_pencomp <-combine.rubin.applied(ate_pencomp_boot)
ate_ptrtzonly <-combine.rubin.applied(ate_ptrtzonly_boot)
ate_ptrtspl <- combine.rubin.applied(ate_ptrtspl_boot)
wtdest_linimpute <- combine.rubin.applied(wtdest_linimpute_boot)
wtdest_pencomp <-combine.rubin.applied(wtdest_pencomp_boot)
wtdest_ptrtzonly <-combine.rubin.applied(wtdest_ptrtzonly_boot)
wtdest_ptrtspl <- combine.rubin.applied(wtdest_ptrtspl_boot)

names(coef_ptrtzonly_boot) <- c("coef0", "coef1", "knotloc0", "knotloc1")
names(coef_ptrtspl_boot) <- c("coef0", "coef1", "knotloc0", "knotloc1")
names(coef_pencomp_boot) <- c("coef0", "coef1", "knotloc0", "knotloc1")

iptwest_ptrtspl <- combine.rubin.applied(iptwest_ptrtspl_boot)
phat <- Map(colMeans, phat_boot)
res <- list(
  ate_linimpute = ate_linimpute,
  ate_pencomp = ate_pencomp,
  ate_ptrtzonly = ate_ptrtzonly,
  ate_ptrtspl = ate_ptrtspl,
  wtdest_linimpute = wtdest_linimpute,
  wtdest_pencomp = wtdest_pencomp ,
  wtdest_ptrtzonly = wtdest_ptrtzonly,
  wtdest_ptrtspl = wtdest_ptrtspl,
  names_fixedcoef = names_fixedcoef,
  pate_synthpop_ptrtspl_boot = pate_synthpop_ptrtspl_boot,
  pate_synthpop_ptrtzonly_boot=pate_synthpop_ptrtzonly_boot,
  iptwest_ptrtspl = iptwest_ptrtspl,
  coef_ptrtzonly = coef_ptrtzonly_boot,
  coef_pencomp = coef_pencomp_boot,
  coef_ptrtspl = coef_ptrtspl_boot,
  phat= phat,
  bootflagct = bootflagct,
  est_ht = est_ht,
  est_naivemod = est_naivemod
)
saveRDS(res,"results/step1_res_faps_wic2.RDS")

lapply(res$coef_ptrtzonly, colMeans)
