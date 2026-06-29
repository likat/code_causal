library(nlme)
library(dplyr)
library(survey)
library(ggplot2)
library(sampling)
library(cmdstanr)
library("rootSolve")
library(LaplacesDemon)
library(parallel)

type_y <- "bin"
type_outcome = "gb"
type_overlap = "notshared"
  pathtocode = "CODE_final/"
  source(paste0(pathtocode,"helpers.R"))
  source(paste0(pathtocode, "calculate_sate.R"))
  source(paste0(pathtocode, "calculate_sate_x.R"))
  source(paste0(pathtocode, "calculate_pencomp_bin_sim.R"))
  source(paste0(pathtocode, "calculate_wfpbb_pencomp_bin_sim_mclapply.R"))
  source(paste0(pathtocode, "imputeF_application.R"))
  source(paste0(pathtocode, "calculate_wtd.R"))
  source(paste0(pathtocode, "calculate_pencomp_nospline_application.R"))
  
  source(paste0(pathtocode, "STEP0-poptingting-bin.R"))
npop <- nrow(pop)
outcome_mod <- NA
propensity_mod <- NA
x_varnames = NA
y_varname = "Ysamp"
catecat_varname = "catecat"
x_name <- c("X1")
x_name_mod <- c("X1")
z_name <- "Z"

### SIMULATION PARS
sim = 200
numboot=200
set.seed(42)
nknot = 15
numclussamp_perstrat <- numclussamp / nstrat

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

cluster_stratum_ref_pop <- pop %>% group_by(stratum, cluster) %>% summarize(cts=n()) 
stratum_ref_pop <- cluster_stratum_ref_pop %>% group_by(stratum) %>% summarize(cts=n()) 
pop$fpc1 <- numclus / nstrat
pop$fpc2 <- npop/numclus
## PROPENSITY GENERATION
if(type_overlap == "shared"){
  propensity_mod <- as.formula("A ~ X1 + X2 + Z + X1:X2")
  pop$ptrt <-  expit(ax1*pop$X1 + ax2*pop$X2 + ax1x2*pop$X1*pop$X2 +
                       az*pop$Z)
  
}else if(type_overlap == "notshared"){
  propensity_mod <- as.formula("A ~ X1 + X2 + X1:X2")
  pop$ptrt <-  expit(ax1*pop$X1 + ax2*pop$X2 + ax1x2*pop$X1*pop$X2)
}

x_varnames <- c(x_name_mod, z_name, paste0(x_name_mod,":",z_name))

varnames_outcome <-  c("X2", "Z")
varnames_outcome_generalize <- c(x_name_mod, z_name, paste0(x_name_mod, ":", z_name))
outcome_mod <- as.formula("Ysamp ~ X2 + Z")
outcome_mod_nospline <- as.formula("Ysamp ~ X1 + X2 + Z")
outcome_mod_generalize <- formulaF(varList=varnames_outcome_generalize, y.name=y_varname)
pate_true <- mean(pop$Y1 - pop$Y0)

### ESTIMATE CONTAINERS
sate <- rep(NA, sim)
resnames <- c("est","withinvar","btwnvar","mivar","ci_tlength","ci_emplength","cicov_t","cicov_emp") # must match names from rubin function
wtdest_nospline <- vector(mode="list", length=8)
ate_nospline <- vector(mode="list", length=8)

names(wtdest_nospline) <- resnames
names(ate_nospline) <- resnames

bootflag <- 0
bootflagct <- 0
# pop$popID <- 1:nrow(pop)
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
# i=1
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


wtdest_nospline_boot <- matrix(NA,nrow=numboot, ncol=2)
ate_nospline_boot <- matrix(NA,nrow=numboot, ncol=2)

# BOOTSTRAP FOR FIRST STAGE IMPUTATION


boot <- 1
while(boot <= numboot){ ## START bootstrap loop
  
  # (1a) Bootstrap the sample (simple stratified bootstrap for non-WFPBB methods)
  bootind_cnt <- sample(ind_cnt, size=size_cnt, replace=T)
  bootind_trt <- sample(ind_trt, size=size_trt, replace=T)
  boot_all <- samp[c(bootind_cnt,bootind_trt),]

  # # step 1 imputation ===
  imputedlist_nospline <- calculate.linimpute(outc_model = outcome_mod,
                                              y.varname = "Ysamp",
                                              trt.varname = "A") 

  ate_nospline_boot[boot,]=c(mean(imputedlist_nospline[[2]]-imputedlist_nospline[[1]]),
                             var(imputedlist_nospline[[2]]-imputedlist_nospline[[1]])/nrow(samp))## variance computed as though SRS

  wtdest_temp <- calculate.wtd(imputedlist = imputedlist_nospline, bootdf=samp)
  
  wtdest_nospline_boot[boot,]  <- c(wtdest_temp[[1]],
                                    wtdest_temp[[2]])

  print(paste0("Boot is ",boot))
  boot <- boot + 1
  bootflag <- 0
} ## END bootstrap


# step 3 rubin's combining rules and store results====
ate_nospline <- Map(c, ate_nospline, combine.rubin(ate_nospline_boot))
wtdest_nospline <- Map(c, wtdest_nospline, combine.rubin(wtdest_nospline_boot))

print(paste0("I is ",i))

# }
res <- list(patetrue = pate_true,
            ate_nospline = ate_nospline,
            wtdest_nospline = wtdest_nospline,
            bootflagct = bootflagct
)
saveRDS(res, file = paste0("results-raw-2stage-bin/linimpute-wrong/results_tingting_noaiptw", type_outcome, type_overlap,i, ".RDS"))

