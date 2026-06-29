## load libraries
library(dplyr)
library(ggplot2)
library(survey)
library(cmdstanr)
library(rootSolve)
cluster = T
#################################################################################
## Read in code
pathtocode = NA
dat = NA
if(cluster){
  pathtocode <- "../causal-bin-sim-new/CODE/"
  dat <- readRDS("dat_use_faps.RDS")
}else{
  pathtocode <- "../causal-bin-sim-new/CODE/"
  dat <- readRDS("data/dat_use_faps.RDS")
}
## AIPTW
# source(paste0(pathtocode, "gcomputeFunc_application.R"))
# source(paste0(pathtocode, "run_AIPTW_application.R"))
source(paste0(pathtocode, "run_AIPTW_new.R"))
source(paste0(pathtocode, "helpers.R"))


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
# ptrtvars <- readRDS("ptrtmodelvars_bic.RDS")
# ptrtmod <- "fsbenefit_12mo ~ pctpov + age + race + educ + hhsize_scaled" %>% as.formula
# ptrtmod <- ptrtvars %>% as.formula

# model.matrix(ptrtmod, dat) %>% colnames
ptrtmod_aiptw <- "fsbenefit_12mo ~ pctpov + age + wicelig + race + educ + as.vector(hhsize_scaled) + age:wicelig"
# outcome_mod <- readRDS("outmodeltxt_faps.RDS") %>% as.formula
outcome_mod <- "adltfs ~ pctpov + age + wicelig + hhsize_scaled + educ + race" %>% as.formula
# outcome_mod <- "adltfs ~ pctpov + age + educ + race" %>% as.formula

names_fixedcoef <- model.matrix(outcome_mod, dat) %>% colnames

outcome_mod_linimpute = outcome_mod
## parameters
treatment_name <- "fsbenefit_12mo"
y_varname <- "adltfs"
numboot = 1000


## Begin analysis
#- Methods we are looking at: 
# classical PENCOMP
# WFPBB-PENCOMP with $Z$ variables only (wtdest)
# WFPBB-PENCOMP using the spline for sampling weights (wtdest)
res_aipw <- list(
  est_wtd_pt = NA,
  est_unwtd_pt = NA,
  est = rep(NA, numboot),
  est_wtd = rep(NA, numboot),
  time = rep(NA, numboot)
)
tempres <- AIPW.application(dataInput = dat,
                               propen.mod = ptrtmod_aiptw,
                               out.name = y_varname,
                               trt.name = treatment_name,
                               modely1 = outcome_mod,
                               modely0 = outcome_mod,
                               outcomefamily = "binomial")
res_aipw$est_wtd_pt <- tempres$pate_wtd
res_aipw$est_unwtd_pt <- tempres$pate_unwtd

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

# trt0_phatrange <- range(dat$phat[which(dat[,treatment_name]==0)])
# trt1_phatrange <- range(dat$phat[which(dat[,treatment_name]==1)])
### ALL OVERLAPPING VALUES
# overlap_samp <- c(max(trt0_phatrange[1], trt1_phatrange[1]),
#                   min(trt0_phatrange[2], trt1_phatrange[2]))
### MORE STRINGENT CUTOFF CRITERIA
# overlap_samp <- c(-3, 3)
# # pdf(file = "results/faps_overlap.pdf",width = 7,height=7)s
# ggplot(dat, aes(x=phat, fill=factor(fsbenefit_12mo))) +
#   geom_histogram(alpha=0.2, position="identity") +
#   labs(x="Estimated propensity score", fill="Trt group")+
#   theme_minimal()
# dev.off()

### the sample should only contain observations that have estimated phat in the overlapping region
# dat <- dat %>% filter(between(phat, overlap_samp[1], overlap_samp[2]))
npop <- sum(dat$wts)
phat_boot <- list(
  pencomp = matrix(nrow = numboot,ncol=nrow(dat)),
  ptrtzonly = matrix(nrow=numboot, ncol = nrow(dat)),
  ptrtspl = matrix(nrow=numboot, ncol=nrow(dat))
)

## Check overlap
# saveRDS(dat,"dat_use_faps_trimmed.RDS")

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


# boot <- Sys.getenv("SLURM_ARRAY_TASK_ID")
# boot <- as.numeric(boot)
bootflagct = 0
# for(boot in 1:numboot){ ## START bootstrap loop
#while(boot <= numboot){
# (1) Bootstrap the sample (stratified by treatment group)
# set.seed(boot)
for( boot in 1:numboot){
bootflag <- 0
# cluster_stratum_ref$resamp <- 0
# cluster_stratum_ref$bbresamp <- 0
# cluster_stratum_ref$clusmultiplier <- 0
# cluster_stratum_ref$bbclusmultiplier <- 0

# dat$clusmultiplier <- NULL
# dat$resamp <- NULL
# dat$bbclusmultiplier <- NULL
# dat$bbresamp <- NULL
# ## with clusters and strata, resample clusters
bootind_cnt <- sample(ind_cnt, size=size_cnt, replace=T)
bootind_trt <- sample(ind_trt, size=size_trt, replace=T)
boot_all <- dat[c(bootind_cnt,bootind_trt),]
boot_all$wts <-  nrow(boot_all) *normalize(boot_all$wts)
boot_all$wts_scaled <- nrow(boot_all) * normalize(boot_all$wts)

temp <- AIPW.application(dataInput = boot_all,
                 propen.mod = ptrtmod_aiptw,
                 out.name = y_varname,
                 trt.name = treatment_name,
                 modely1 = outcome_mod,
                 modely0 = outcome_mod,
                 outcomefamily = "binomial")

res_aipw$est[boot] <- temp$pate_unwtd
res_aipw$est_wtd[boot] <- temp$pate_wtd
res_aipw$time[boot] <- temp$time

print(paste0("Boot is ",boot))
bootflag <- 0
} ## END bootstrap


### Combine
# res <- list(est = aiptw_boot, time = time1-time0)
saveRDS(res_aipw,paste0("results/step1_res_faps_aipw_new.RDS"))

