iter = 200
library(dplyr)
# ate_linimpute = ate_linimpute,
# ate_pencomp = ate_pencomp,
# ate_ptrtzonly = ate_ptrtzonly,
# ate_ptrtspl = ate_ptrtspl,
# wtdest_linimpute = wtdest_linimpute,
# wtdest_pencomp = wtdest_pencomp ,
# wtdest_ptrtzonly = wtdest_ptrtzonly,
# wtdest_ptrtspl = wtdest_ptrtspl,
# names_fixedcoef = names_fixedcoef,
# pate_synthpop_ptrtspl_boot = pate_synthpop_ptrtspl_boot,
# pate_synthpop_ptrtzonly_boot=pate_synthpop_ptrtzonly_boot,
# iptwest_ptrtspl = iptwest_ptrtspl,
# coef_ptrtzonly = coef_ptrtzonly_boot,
# coef_pencomp = coef_pencomp_boot,
# coef_ptrtspl = coef_ptrtspl_boot,
resnames <- c("est","withinvar","btwnvar","mivar","ci_tlength","ci_emplength") # must match names from rubin function

helperfun <- function(tbl, wfpbb = FALSE){
  ## confidence intervals require uniform prior 
  res <- vector(mode = "list", length = 5)
  names(res) <- c("est", "withinvar", "btwnvar", "cilower", "ciupper")
  res$est <- mean(tbl$est)
  res$withinvar <- mean(tbl$withinvar)
  res$btwnvar <- var(tbl$est)
  tdf <- 57-25
  if(wfpbb == T){
    mivar <- (1+1/(iter-1))*res$btwnvar

  }else{
    mivar <- res$withinvar + (1+1/iter)*res$btwnvar
  }
  
  cilower <- res$est - qt(0.975,tdf)*sqrt(mivar)
  ciupper <- res$est + qt(0.975,tdf)*sqrt(mivar)
  # ci_all <- ciunifun(estvar = mivar, estprop = res$est)
  
  res <- list(est=res$est, withinvar=res$withinvar,btwnvar=res$btwnvar, mivar = mivar,
              cilower =  cilower, ciupper =ciupper)
  return(res)
}
  ate_linimpute <- vector(mode = "list", length = 6)
  ate_pencomp <- vector(mode="list", length=6)
  ate_ptrtzonly <- vector(mode="list", length=6)
  ate_ptrtspl <- vector(mode="list", length=6)
  wtdest_linimpute <- vector(mode = "list", length = 6)
  wtdest_pencomp <- vector(mode="list", length=6)
  wtdest_ptrtzonly <- vector(mode="list", length=6)
  wtdest_ptrtspl<- vector(mode="list", length=6)
  synthpop_ptrtzonly <- rep(NA, 220)
  synthpop_ptrtspl <- rep(NA,220)
  synthpop_ptrtzonly_time <- rep(NA, 220)
  synthpop_ptrtspl_time <- rep(NA,220)
  
  itersum <- 0
  skipnum = 0
  for(i in 1:220){
    temp <- tryCatch ( {
      res <- readRDS(paste0("results-raw/step1_res_faps",i,".RDS"))
    }, error=function(e) return(NA) )
    if(sum(is.na(temp))>0){skipnum = skipnum + 1; next}
    res <- readRDS(paste0("results-raw/step1_res_faps",i,".RDS"))
    
    
    ate_linimpute <- Map(c, ate_linimpute, res$ate_linimpute)
    ate_pencomp <- Map(c, ate_pencomp, res$ate_pencomp)
    ate_ptrtzonly <- Map(c, ate_ptrtzonly, res$ate_ptrtzonly)
    ate_ptrtspl <- Map(c, ate_ptrtspl, res$ate_ptrtspl)
    wtdest_linimpute <- Map(c, wtdest_linimpute, res$wtdest_linimpute)
    wtdest_pencomp <- Map(c, wtdest_pencomp, res$wtdest_pencomp)
    wtdest_ptrtzonly <- Map(c, wtdest_ptrtzonly, res$wtdest_ptrtzonly)
    wtdest_ptrtspl<- Map(c, wtdest_ptrtspl, res$wtdest_ptrtspl)
    synthpop_ptrtspl[i] <- res$pate_synthpop_ptrtspl_boot
    synthpop_ptrtzonly[i] <- res$pate_synthpop_ptrtzonly_boot
    synthpop_ptrtzonly_time[i] <- res$time_wfpbbzonly
    synthpop_ptrtspl_time[i] <- res$time_wfpbbspl
    itersum = itersum + 1
    if(itersum >= 200){break}
  }
  res_synthpop_ptrtspl <- data.frame(
    est = synthpop_ptrtspl[!is.na(synthpop_ptrtspl)],
    withinvar = NA
  )
  res_synthpop_ptrtzonly <- data.frame(
    est =synthpop_ptrtzonly[!is.na(synthpop_ptrtzonly)],
    withinvar = NA
  )
  
  names(ate_linimpute) <- resnames
  names(ate_pencomp) <- resnames
  names(ate_ptrtzonly) <- resnames
  names(ate_ptrtspl) <- resnames
  names(wtdest_linimpute) <- resnames
  names(wtdest_pencomp) <- resnames
  names(wtdest_ptrtzonly) <- resnames
  names(wtdest_ptrtspl) <- resnames
  # helperfun(wtdest_ptrtspl, wfpbb=T)
  # helperfun(wtdest_pencomp, wfpbb = F)
  # helperfun(wtdest_ptrtzonly, wfpbb = T)
  reswrite <- list(
    ate_linimipute = ate_linimpute %>% helperfun,
    ate_pencomp= ate_pencomp %>% helperfun,
    ate_ptrtzonly = ate_ptrtzonly %>% helperfun(wfpbb=T),
    ate_ptrtspl = ate_ptrtspl %>% helperfun(wfpbb = T),
    wtd_linimpute = wtdest_linimpute %>% helperfun,
    wtd_pencomp= wtdest_pencomp %>% helperfun,
    wtd_ptrtzonly = wtdest_ptrtzonly %>% helperfun(wfpbb=T),
    wtd_ptrtspl = wtdest_ptrtspl %>% helperfun(wfpbb = T),
    synthpop_ptrtspl = res_synthpop_ptrtspl %>% helperfun(wfpbb=T),
    synthpop_ptrtzonly = res_synthpop_ptrtzonly %>% helperfun(wfpbb=T)
  )
saveRDS(reswrite, "results-clean/res_noaiptw_application.RDS")

temp <- readRDS("results-clean/res_noaiptw_application.RDS")
