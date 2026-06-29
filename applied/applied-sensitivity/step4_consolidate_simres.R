## Read in data
library(dplyr)
library(ggplot2)
library(survey)
library(cmdstanr)
require(LaplacesDemon)
require(polyapost)
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
  trtment = rep(trtment_u, each = 3),
  est_pate = NA,
  se_pate = NA
)
case_df$name <- paste0(abs(case_df$outcome)*10, "_",abs(case_df$trtment)*10)

for(i in 1:nrow(case_df)){
  case <- case_df[i,]
  
  res <- readRDS(paste0("results/step3_simulation_", case$name,".RDS"))
  case_df$est_pate[i] <- mean(res$pate_synthpop_ptrtzonly_boot, na.rm=T)
  case_df$se_pate[i] <- sqrt((1+1/200)*var(res$pate_synthpop_ptrtzonly_boot, na.rm=T))
}
case_df <- readRDS("results/step3_simulation_consolidated.RDS")
case_df <- case_df %>% mutate(
  cilower = est_pate - qt(0.975,57-25)*se_pate,
  ciupper = est_pate + qt(0.975,57-25)*se_pate
) %>% arrange(outcome)

saveRDS(case_df, "results/step3_simulation_consolidated.RDS")
