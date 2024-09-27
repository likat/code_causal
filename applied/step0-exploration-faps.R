library(haven)
library(survey)
library(dplyr)
library(tidyr)
faps <- read.csv("data/faps_household_puf.csv")
faps_indv <- read.csv("data/faps_individual_puf.csv")
faps_indv <- faps_indv %>% filter(AGE_R != "D") %>% mutate(AGE_R = as.numeric(AGE_R))
faps_wts <- read.csv("data/faps_hhweights.csv")

## Rename some columns for easier ref
# focus on adult food security (child food security probably more related to chars of the adults they are staying w)

## ROUND 1 CLEANUP: renaming variables
datfull <- faps %>% 
  left_join(faps_wts) %>% 
  left_join(faps_indv, by=c("hhnum"="HHNUM")) %>% 
  filter(AGE_R >= 18) %>% 
  filter(RELATION_R ==0) %>% 
  mutate(wts = HHWGT,
         sex = case_when(
           SEX == "2" ~ "female",
           SEX == "1" ~ "male"),
         pctpov = case_when(# 349 missing, probably better to use INDFMPIR (pctpov to poverty guideline, but greater number of NA)
           pctpovguidehh_r <= 100 ~ "below 100%",
           pctpovguidehh_r <= 184 ~ "100-184%",
           pctpovguidehh_r >= 185 ~ "185%+")
         , 
         wicelig = case_when( # eligeble for WIC?
           wiccategelig == 0 ~ "No",
           wiccategelig == 1 ~ "Yes"
         ),
         age = case_when(
          AGE_R >= 18 & AGE_R <= 35 ~ "18-35",
          AGE_R >= 36 & AGE_R <= 65 ~ "35-65",
           AGE_R > 65 ~ "65+"
           ),
         race = case_when(
           RACECAT_R == 1 & HISPANIC == 0 ~ "Non-hispanic White",
           RACECAT_R == 2 & HISPANIC == 0 ~ "Non-hispanic Black",
           RACECAT_R %in% c(3,4,5,6) & HISPANIC == 0 ~ "Non-hispanic Other",
           HISPANIC == 1 ~ "Hispanic"
         ),
         educ = case_when(
           EDUCCAT %in% c(1,2) ~ "No HS diploma",
           EDUCCAT == 3~ "HS grad or GED",
           EDUCCAT == 4 ~ "Some college",
           EDUCCAT %in% c(5,6) ~ "College graduate or above"
         ), ## 13 refused or don't know 
         hhsize = hhsize,
         primstoretraveltime = case_when(
           primstoretraveltime >=0 ~ primstoretraveltime,
           TRUE ~ NA
           ),
         adltfs  = case_when(
           adltfscat%in% c(1,2) ~ 0, ## food secure,
           adltfscat %in% c(3,4) ~ 1 ## food insecure
         ), ## 325 missing, missingness correlates with pctpov missingness
         fsbenefit_12mo = case_when(
           snapnowhh == 1 ~ 1,
           snapever==1 & snap12mos==1 ~ 1,
           usdafoods==1 ~ 1,
           wichh==1 ~ 1,
           mealfacility==1 ~ 1,
           mealdelivery==1 ~1,
           snapnowhh == 0 ~ 0,
           snapever==0 & snap12mos==0 ~ 0,
           usdafoods==0 ~ 0,
           wichh==0 ~ 0,
           mealfacility==0 ~ 0,
           mealdelivery==0 ~0
         ) ## Have you received a fs benefit in the past 12 mo?
    )



table(datfull$pctpov, useNA= "always")
table(datfull$age, useNA= "always")
table(datfull$race, useNA= "always")
table(datfull$educ, useNA= "always")
table(datfull$hhsize, useNA= "always")
datfull$wts_scaled <- nrow(datfull) * datfull$wts/sum(datfull$wts)
datpool <- datfull %>% dplyr::select(
  wts, pctpov, age, sex,race, wicelig, educ, hhsize, primstoretraveltime, adltfs, fsbenefit_12mo,
  tspsu, tsstrata, hhwgt1:hhwgt57
) %>% drop_na
datpool$sex <- factor(datpool$sex, levels = c("male","female"))
datpool$traveltime_scaled <- scale(datpool$primstoretraveltime)
datpool$hhsize_scaled <- scale(datpool$hhsize)
datpool$wts_scaled <- nrow(datpool)* datpool$wts/sum(datpool$wts)
## BIC model selection
modelbase <- glm(fsbenefit_12mo ~pctpov + age +sex+sex*age+ wicelig + race + educ + hhsize_scaled + traveltime_scaled, datpool,
                 family = "binomial")
modelbase_outcome <- glm(adltfs~pctpov + age + wicelig +sex+sex*age+ race + educ + hhsize_scaled  + traveltime_scaled + fsbenefit_12mo,
                         family = "binomial",
                         data = datpool)
modeltrt.bic <- MASS::stepAIC(modelbase,
                        scope = list(
                          upper = ~(pctpov + age  +sex+sex*age+ wicelig+ race + educ + hhsize_scaled + traveltime_scaled)^2,
                          lower = ~1), k=log(nrow(datpool)))
modelout.bic <- MASS::stepAIC(modelbase_outcome,
                              scope = list(
                                upper = ~(pctpov + age+ wicelig  +sex+sex*age+ race + educ + hhsize_scaled + traveltime_scaled + fsbenefit_12mo)^2,
                                lower = ~fsbenefit_12mo
                                ), k=log(nrow(datpool)))
# modelout.aic <- MASS::stepAIC(modelbase_outcome,
#                               scope = list(
#                                 upper = ~(pctpov + age  +sex+sex*age+ race + educ + hhsize_scaled + traveltime_scaled + fsbenefit_12mo)^2,
#                                 lower = ~fsbenefit_12mo
#                               ))
# library(glmnet)
# y <- datpool$adltfs %>% as.numeric
# x <- model.matrix(~(pctpov + age + race + educ + fsbenefit_12mo + hhsize_scaled)^2-1, data = datpool)
# lassofit <- glmnet(x=x, y=y, family = "binomial")
# cvfit <- cv.glmnet(x=x, y=y, family = "binomial")
# plot(cvfit)
# lambdamin <- cvfit$lambda.1se
# coef(cvfit,s=lambdamin )
saveRDS(modeltrt.bic, "ptrtmodelvars_bic.RDS")
saveRDS(modelout.bic, "outmodelvars_bic.RDS")
# saves interaction terms for education, age:householdsize, and pctpov:householdsize
## Tableone

## CHECKS ------
## (1) Is there treatment heterogeneity?
svydes <- svydesign(ids=~0, weights=~wts,
                    data=datpool)
svydes_mod <- svydesign(ids=~0, weights=~wts_scaled,
                        data=datpool)
svydesunwtd <- svydesign(ids=~0,data = datpool)
model0_text <- 
  "adltfs~pctpov + age + sex + educ + fsbenefit_12mo + hhsize_scaled + race + age*sex + race*fsbenefit_12mo" %>% 
  as.formula
model0 <- svyglm(model0_text, 
       design = svydes_mod, family = "binomial")
summary(model0)
model0unwtd <- svyglm(model0_text, 
                 design = svydesunwtd, family = "binomial")
summary(model0unwtd)

svymean(~adltfs, design = svydes)
svyby(~adltfs, by=~fsbenefit_12mo,FUN = svymean, design = svydes)
svymean(~adltfs, design = svydesunwtd)
svyby(~adltfs, by=~fsbenefit_12mo,FUN = svymean, design = svydesunwtd)

## (2) Is the sampling informative?
modelwtsy <- svyglm(adltfs~ wts_scaled, design = svydes_mod, family = "binomial")
summary(modelwtsy)

## (3) Is there overlap between sampling and treatment variables?
modelwtsz <- svyglm(wts_scaled ~ race + pctpov + hhsize_scaled, 
                  design = svydes)
summary(modelwtsz)
modelaz <- svyglm(fsbenefit_12mo ~ race + pctpov + hhsize_scaled, 
                  design = svydes, family = "binomial")
summary(modelaz)

## Save cleaned dataset
saveRDS(datpool, "data/faps_cleaned.RDS")
model0_textpencomp <- 
  "adltfs~pctpov + age  + educ + hhsize_scaled + race" %>% 
  as.formula
saveRDS(model0_textpencomp, "outmodeltxt_faps.RDS")
