## load libraries
library(dplyr)
library(ggplot2)
library(survey)
library(cmdstanr)
library(LaplacesDemon)

#################################################################################
## Read in code
pathtocode <- "./"
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
dat <- readRDS("dat_use_faps.RDS") 
dat$tsstrata <- dplyr::dense_rank(dat$tsstrata)
dat$tspsu <- dplyr::dense_rank(dat$tspsu)

## Arrange
nstrata <- dat$tsstrata %>% unique() %>% length()
cluster_stratum_ref <- dat %>% group_by(tsstrata, tspsu) %>% summarize(cts=n()) 
cluster_stratum_ref <- cluster_stratum_ref %>% group_by(tsstrata) %>% mutate(nclus = length(unique(tspsu)),
                                                                             mstar = nclus / (nclus-1)) %>% 
  select(-cts)
ptrtvars <- readRDS("ptrtmodelvars_bic.RDS")
# ptrtmod <- "fsbenefit_12mo ~ pctpov + age + race + educ + hhsize_scaled" %>% as.formula
ptrtmod <- ptrtvars %>% as.formula

model.matrix(ptrtmod, dat) %>% colnames
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

est_naivemod_boot <- rep(NA, numboot)
est_naiveht_boot <- rep(NA, numboot)


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
# phat_imputed <- calculate.ptrtspline(fitdat = dat, 
#                                      sampdat= dat,
#                                      num.knot = 10,
#                                      propen.model = ptrtmod,
#                                      trtname=treatment_name)
# dat$phat <- phat_imputed$ptrt_fitdat
# 
# # trt0_phatrange <- range(dat$phat[which(dat[,treatment_name]==0)])
# # trt1_phatrange <- range(dat$phat[which(dat[,treatment_name]==1)])
# ### ALL OVERLAPPING VALUES
# # overlap_samp <- c(max(trt0_phatrange[1], trt1_phatrange[1]),
# #                   min(trt0_phatrange[2], trt1_phatrange[2]))
# ### MORE STRINGENT CUTOFF CRITERIA
# overlap_samp <- c(-3, 3)
# # pdf(file = "results/faps_overlap.pdf",width = 7,height=7)s
# ggplot(dat, aes(x=phat, fill=factor(fsbenefit_12mo))) +
#   geom_histogram(alpha=0.2, position="identity") +
#   labs(x="Estimated propensity score", fill="Trt group")+
#   theme_minimal()
# # dev.off()
# 
# ### the sample should only contain observations that have estimated phat in the overlapping region
# dat <- dat %>% filter(between(phat, overlap_samp[1], overlap_samp[2]))
npop <- sum(dat$wts)
phat_boot <- list(
  pencomp = matrix(nrow = numboot,ncol=nrow(dat)),
  ptrtzonly = matrix(nrow=numboot, ncol = nrow(dat)),
  ptrtspl = matrix(nrow=numboot, ncol=nrow(dat))
)

## Check overlap
# saveRDS(dat,"data/dat_use_faps.RDS")

## GET HT estimate
wtdes <- svydesign(ids=~tspsu,strata=~tsstrata,weights=~wts,data=dat)
wtdes_rep <- as.svrepdesign(wtdes,type = "subbootstrap")
est_ht <- svyby(~adltfs,by=~fsbenefit_12mo, design=wtdes, FUN=svymean)
est_ht <- est_ht$adltfs[2]-est_ht$adltfs[1]

## GET Naive regression estimate
mod0 <- glm(adltfs~pctpov + age + educ + hhsize_scaled + race + wicelig+age:wicelig+fsbenefit_12mo, data = dat,family = "binomial")
# mod0 <- svyglm(formula = "adltfs~pctpov + age + educ + hhsize_scaled + race + wicelig+age:wicelig+fsbenefit_12mo",design = wtdes,family = "binomial")
dat$temppred <- predict(mod0,newdata=dat, type="response")
wtdes_mod2 <- svydesign(ids=~0,data=dat)
est_naivemod <- svyby(formula=~temppred,by=~fsbenefit_12mo, design = wtdes_mod2,FUN=svymean)
est_naivemod <- est_naivemod$temppred[2]-est_naivemod$temppred[1]
## Begin bootstrap
boot <- 1
# boot <- Sys.getenv("SLURM_ARRAY_TASK_ID")
boot <- as.numeric(boot)
bootflagct = 0
# for(boot in 1:numboot){ ## START bootstrap loop
for(boot in 1:numboot){
  # (1) Bootstrap the sample (stratified by treatment group)
  set.seed(boot)
  
  bootflag <- 0
  # cluster_stratum_ref$resamp <- 0
  # cluster_stratum_ref$bbresamp <- 0
  # cluster_stratum_ref$clusmultiplier <- 0
  # cluster_stratum_ref$bbclusmultiplier <- 0
  # 
  # dat$clusmultiplier <- NULL
  # dat$resamp <- NULL
  # dat$bbclusmultiplier <- NULL
  # dat$bbresamp <- NULL
  # ## with clusters and strata, resample clusters
  # for(h in 1:nstrata){
  #   cluspool <- cluster_stratum_ref$tspsu[cluster_stratum_ref$tsstrata==h]
  #   tempnumclus <- length(cluspool)
  #   resamp <- sample(cluspool,tempnumclus-1,replace=T)
  #   cluster_stratum_ref$resamp[cluster_stratum_ref$tspsu %in% unique(resamp)] <- table(resamp)
  # }
  # 
  # for(k in 1:nstrata){
  #   cluspool <- cluster_stratum_ref$tspsu[cluster_stratum_ref$tsstrata==k]
  #   tempnumclus <- length(cluspool)
  #   bootstrap <- BayesianBootstrap(cluspool, 1)
  #   resamp <- rmultinom(1, tempnumclus - 1, bootstrap)
  #   cluster_stratum_ref$bbresamp[cluster_stratum_ref$tspsu %in% cluspool] <- resamp
  # }
  # cluster_stratum_ref$clusmultiplier <- cluster_stratum_ref$resamp*cluster_stratum_ref$mstar
  # cluster_stratum_ref$bbclusmultiplier <- cluster_stratum_ref$bbresamp*cluster_stratum_ref$mstar
  # 
  # dat <- dat %>% left_join(cluster_stratum_ref %>% select(tspsu, tsstrata,clusmultiplier,bbclusmultiplier,
  #                                                         resamp, bbresamp))
  ## with clusters and strata, resample clusters
  bootind_cnt <- sample(ind_cnt, size=size_cnt, replace=T)
  bootind_trt <- sample(ind_trt, size=size_trt, replace=T)
  boot_all <- dat[c(bootind_cnt,bootind_trt),]
  
  # boot_all <- dat[rep(1:nrow(dat), dat$resamp),]
  boot_all$wts <- sum(dat$wts)*normalize(boot_all$wts)
  ## Get naive model estimator
  wtdes2 <- svydesign(ids=~0, weights=~wts,data=boot_all)
  mod0 <- glm(adltfs~pctpov + age + educ + hhsize_scaled + race + wicelig+age:wicelig+fsbenefit_12mo,data = boot_all, family = "binomial")
  boot_all$temppred <- predict(mod0,newdata=boot_all, type="response")
  wtdes2_mod<- svydesign(ids=~0, data=boot_all)
  est_naivemodtemp <- svyby(formula=~temppred,by=~fsbenefit_12mo, design = wtdes2_mod,FUN=svymean)
  est_naivemodtemp <- est_naivemodtemp$temppred[2]-est_naivemodtemp$temppred[1]
  
  est_naivemod_boot[boot] <-  est_naivemodtemp
  est_httemp <- svyby(~adltfs,by=~fsbenefit_12mo, design=wtdes2, FUN=svymean)
  est_naiveht_boot[boot] <- est_httemp$adltfs[2]-est_httemp$adltfs[1]
  
  
  # # step 1 imputation ===
  imputedlist_linimpute <- calculate.linimpute(outc_model = outcome_mod_linimpute,
                                               y.varname = y_varname,
                                               trt.varname = treatment_name)

    #if(bootflag == 1){bootflagct = bootflagct + 1; next}
  
  #   # Calculate ATE
  ate_linimpute_boot[boot,] = c(mean(imputedlist_linimpute[[2]]-imputedlist_linimpute[[1]]),
                                var(imputedlist_linimpute[[2]]-imputedlist_linimpute[[1]])/nrow(dat))
  # ate_pencomp_boot[1,]=c(mean(imputedlist_pencomp[[2]]-imputedlist_pencomp[[1]]),
  #                           var(imputedlist_pencomp[[2]]-imputedlist_pencomp[[1]])/nrow(dat))
  # ate_ptrtzonly_boot[1,]=c(mean(imputedlist_ptrtzonly[[2]]-imputedlist_ptrtzonly[[1]]), 
  #                             var(imputedlist_ptrtzonly[[2]]-imputedlist_ptrtzonly[[1]])/nrow(dat))
  # ate_ptrtspl_boot[1,]=c(mean(imputedlist_ptrtspl[[2]]-imputedlist_ptrtspl[[1]]),
  #                           var(imputedlist_ptrtspl[[2]]-imputedlist_ptrtspl[[1]])/nrow(dat))

  # ate_ptrtzonly_boot[1,]=c(imputedlist_ptrtzonly$ate_samp %>% unlist %>% as.numeric)
  # ate_ptrtspl_boot[1,]=c(imputedlist_ptrtspl$ate_samp %>% unlist %>% as.numeric)
  
  # step 2 generalize! ====
  wtdres_linimpute <- calculate.wtd(imputedlist = imputedlist_linimpute, bootdf=dat)
  # wtdres_pencomp <- calculate.wtd(imputedlist = imputedlist_pencomp, bootdf=dat)
  # wtdres_ptrtzonly <- calculate.wtd(imputedlist = imputedlist_ptrtzonly, bootdf=dat)
  # wtdres_ptrtspl <- calculate.wtd(imputedlist = imputedlist_ptrtspl, bootdf=dat)
  
  wtdest_linimpute_boot[boot,] <- c(wtdres_linimpute[[1]],
                                    wtdres_linimpute[[2]])
  # wtdest_pencomp_boot[1,] <- c(wtdres_pencomp[[1]],
  #                                 wtdres_pencomp[[2]])
  # wtdest_ptrtzonly_boot[1,] <- c(imputedlist_ptrtzonly$wtdest[[1]], imputedlist_ptrtzonly$wtdest[[2]])
  # wtdest_ptrtspl_boot[1,] <- c(imputedlist_ptrtspl$wtdest[[1]], imputedlist_ptrtspl$wtdest[[2]])
  # iptwest_ptrtspl_boot[1,] <- c(imputedlist_ptrtspl$est_iptw,
  #                                  as.numeric(imputedlist_ptrtspl$var_iptw))
  
  # coef_pencomp_boot <- Map(rbind,coef_pencomp_boot, imputedlist_pencomp$coeflistres)
  # coef_ptrtspl_boot <- Map(rbind,coef_ptrtspl_boot, imputedlist_ptrtspl$coeflistres)
  # coef_ptrtzonly_boot <- Map(rbind,coef_ptrtzonly_boot, imputedlist_ptrtzonly$coeflistres)
  # phat_boot$pencomp[1,] <- imputedlist_pencomp$phat_samp
  # phat_boot$ptrtzonly[1,] <- imputedlist_ptrtzonly$phat_samp
  # phat_boot$ptrtspl[1,] <- imputedlist_ptrtspl$phat_samp
  print(paste0("Boot is ",boot))
  bootflag <- 0
} ## END bootstrap

### Combine
ate_linimpute <- combine.rubin.applied(ate_linimpute_boot)
# ate_pencomp <-combine.rubin.applied(ate_pencomp_boot)
# ate_ptrtzonly <-combine.rubin.applied(ate_ptrtzonly_boot)
# ate_ptrtspl <- combine.rubin.applied(ate_ptrtspl_boot)
wtdest_linimpute <- combine.rubin.applied(wtdest_linimpute_boot)
# wtdest_pencomp <-combine.rubin.applied(wtdest_pencomp_boot)
# wtdest_ptrtzonly <-combine.rubin.applied(wtdest_ptrtzonly_boot)
# wtdest_ptrtspl <- combine.rubin.applied(wtdest_ptrtspl_boot)

# names(coef_ptrtzonly_boot) <- c("coef0", "coef1", "knotloc0", "knotloc1")
# names(coef_ptrtspl_boot) <- c("names","coef0", "coef1", "knotloc0", "knotloc1")
# names(coef_pencomp_boot) <- c("fixed0", "fixed1","rand0", "rand1","knotloc0", "knotloc1")

# iptwest_ptrtspl <- combine.rubin.applied(iptwest_ptrtspl_boot)
# phat <- Map(colMeans, phat_boot)
res <- list(
  ate_linimpute = ate_linimpute,
  wtdest_linimpute = wtdest_linimpute,
  est_ht = est_ht,
  est_naivemod = est_naivemod,
  se_ht = sqrt(var(est_naiveht_boot)),
  se_mod = sqrt(var(est_naivemod_boot))
)
res$cilower_ht <- res$est_ht - qt(0.975, min(size_cnt-1, size_trt-1))*res$se_ht
res$ciupper_ht <- res$est_ht + qt(0.975, min(size_cnt-1, size_trt-1))*res$se_ht
res$cilower_naivemod <- res$est_naivemod - qt(0.975, min(size_cnt-1, size_trt-1))*res$se_mod
res$ciupper_naivemod <- res$est_naivemod + qt(0.975, min(size_cnt-1, size_trt-1))*res$se_mod

saveRDS(res,paste0("results-clean/step1_res_faps_naive.RDS"))

