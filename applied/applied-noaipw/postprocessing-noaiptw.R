iter = 200
library(dplyr)
aiptwres.new <- readRDS("../causal-application/results-clean/res_aipw_new.RDS")
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
    mivar <- (1+1/(iter))*res$btwnvar

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
helperfun.coef <- function(tbl, wfpbb = FALSE){
  ## NOTE: IF REPORTING RESULTS FOR NONWFPBB MODELS, THEN NEED TO INCORPORATE WITHIN VAR
  ## confidence intervals require uniform prior 
  res <- vector(mode = "list", length = 4)
  names(res) <- c("est", "btwnvar", "cilower", "ciupper")
  res$est <- colMeans(tbl)
  res$btwnvar <- apply(tbl, 2, var)
  tdf <- 57-25
    mivar <- (1+1/(iter))*res$btwnvar
  cilower <- res$est - qt(0.975,tdf)*sqrt(mivar)
  ciupper <- res$est + qt(0.975,tdf)*sqrt(mivar)
  # ci_all <- ciunifun(estvar = mivar, estprop = res$est)
  
  res <- list(est=res$est, btwnvar=res$btwnvar, SE = sqrt(mivar),
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
  coef_pencomp <- list(
    modtrt = matrix(nrow = iter, ncol = 15),
    fixed0 = matrix(nrow = iter, ncol = 14),
    fixed1 = matrix(nrow = iter, ncol = 14),
    rand0 = matrix(nrow = iter, ncol = 10),
    rand1 = matrix(nrow = iter, ncol = 10)
  )
  coef_linimpute <- list(
    fixed0 = matrix(nrow = iter, ncol = 13),
    fixed1 = matrix(nrow = iter, ncol = 13)  )
  coef_ptrtzonly <- list(
    modtrt = matrix(nrow = iter, ncol = 15),
    fixed0 = matrix(nrow = iter, ncol = 14),
    fixed1 = matrix(nrow = iter, ncol = 14),
    rand0 = matrix(nrow = iter, ncol = 10),
    rand1 = matrix(nrow = iter, ncol = 10),
    sigma_spline0 = rep(NA, iter),
    sigma_spline1 = rep(NA, iter)
  )
  
  synthpop_ptrtzonly <- rep(NA, 220)
  synthpop_ptrtspl <- rep(NA,220)
  synthpop_ptrtzonly_time <- rep(NA, 220)
  synthpop_ptrtspl_time <- rep(NA,220)
  
  itersum <- 0
  skipnum = 0
  for(i in 1:220){
    temp <- tryCatch ( {
      res <- readRDS(paste0("results-raw/step1_res_faps_noptrtspl",i,".RDS"))
    }, error=function(e) return(NA) )
    if(sum(is.na(temp))>0){skipnum = skipnum + 1; next}
    res <- readRDS(paste0("results-raw/step1_res_faps_noptrtspl",i,".RDS"))
    # res <- readRDS(paste0("results-raw/step1_res_faps_noptrtspl", i, ".RDS"))
    # res_sigspl <- readRDS(paste0("results-raw-wfpbbpencomponly/step1_res_faps_noptrtspl_numknot10_", i, ".RDS"))
    ate_linimpute <- Map(c, ate_linimpute, res$ate_linimpute)
    ate_pencomp <- Map(c, ate_pencomp, res$ate_pencomp)
    ate_ptrtzonly <- Map(c, ate_ptrtzonly, res$ate_ptrtzonly)
    ate_ptrtspl <- Map(c, ate_ptrtspl, res$ate_ptrtspl)
    coef_linimpute$fixed0[itersum+1,] <- res$coef_linimpute[[1]]
    coef_linimpute$fixed1[itersum+1,] <- res$coef_linimpute[[2]]

    coef_pencomp$modtrt[itersum+1,] <- res$coef_pencomp[[1]]
    coef_pencomp$fixed0[itersum+1,] <- res$coef_pencomp[[2]]
    coef_pencomp$fixed1[itersum+1,] <- res$coef_pencomp[[3]]
    coef_pencomp$rand0[itersum+1,] <- res$coef_pencomp[[4]]
    coef_pencomp$rand1[itersum+1,] <-res$coef_pencomp[[5]]
    
    coef_ptrtzonly$modtrt[itersum+1,] <- res$coef_ptrtzonly[[1]]
    coef_ptrtzonly$fixed0[itersum+1,] <- res$coef_ptrtzonly[[2]]
    coef_ptrtzonly$fixed1[itersum+1,] <- res$coef_ptrtzonly[[3]]
    coef_ptrtzonly$rand0[itersum+1,] <- res$coef_ptrtzonly[[4]]
    coef_ptrtzonly$rand1[itersum+1,] <-res$coef_ptrtzonly[[5]]
    coef_ptrtzonly$sigma_spline0[itersum+1] <- res$coef_ptrtzonly[[6]][1,1] %>% unlist %>% as.numeric()
    coef_ptrtzonly$sigma_spline1[itersum+1] <- res$coef_ptrtzonly[[6]][1,2]%>% unlist %>% as.numeric()
    
    wtdest_linimpute <- Map(c, wtdest_linimpute, res$wtdest_linimpute)
    wtdest_pencomp <- Map(c, wtdest_pencomp, res$wtdest_pencomp)
    wtdest_ptrtzonly <- Map(c, wtdest_ptrtzonly, res$wtdest_ptrtzonly)
    # wtdest_ptrtspl<- Map(c, wtdest_ptrtspl, res$wtdest_ptrtspl)
    # synthpop_ptrtspl[i] <- res$pate_synthpop_ptrtspl_boot
    synthpop_ptrtzonly[i] <- res$pate_synthpop_ptrtzonly_boot
    synthpop_ptrtzonly_time[i] <- res$time_wfpbbzonly
    # synthpop_ptrtspl_time[i] <- res$time_wfpbbspl
    itersum = itersum + 1
    if(itersum >= 200){break}
  }
  # res_synthpop_ptrtspl <- data.frame(
  #   est = synthpop_ptrtspl[!is.na(synthpop_ptrtspl)],
  #   withinvar = NA
  # )
  res_synthpop_ptrtzonly <- data.frame(
    est =synthpop_ptrtzonly[!is.na(synthpop_ptrtzonly)],
    withinvar = NA
  )
  names(ate_linimpute) <- resnames
  names(ate_pencomp) <- resnames
  names(ate_ptrtzonly) <- resnames
  # names(ate_ptrtspl) <- resnames
  names(wtdest_linimpute) <- resnames
  names(wtdest_pencomp) <- resnames
  names(wtdest_ptrtzonly) <- resnames
  # names(wtdest_ptrtspl) <- resnames
  # helperfun(wtdest_ptrtspl, wfpbb=T)
  # helperfun(wtdest_pencomp, wfpbb = F)
  # helperfun(wtdest_ptrtzonly, wfpbb = T)
  reswrite <- list(
    ate_linimipute = ate_linimpute %>% helperfun,
    ate_pencomp= ate_pencomp %>% helperfun,
    ate_ptrtzonly = ate_ptrtzonly %>% helperfun(wfpbb=T),
    # ate_ptrtspl = ate_ptrtspl %>% helperfun(wfpbb = T),
    coef_linimpute = lapply(coef_linimpute, helperfun.coef),
    coef_pencomp= lapply(coef_pencomp, helperfun.coef),
    coef_ptrtzonly = lapply(coef_ptrtzonly[1:5], helperfun.coef),
    sigmaspl_ptrtzonly = c(mean(coef_ptrtzonly[[6]]), mean(coef_ptrtzonly[[7]])),
    sigmaspl_ptrtzonly_sd = sqrt(c((1+1/(iter))*var(coef_ptrtzonly[[6]]), (1+1/(iter))*var(coef_ptrtzonly[[7]]))),
    wtd_linimpute = wtdest_linimpute %>% helperfun,
    wtd_pencomp= wtdest_pencomp %>% helperfun,
    wtd_ptrtzonly = wtdest_ptrtzonly %>% helperfun(wfpbb=T),
    # wtd_ptrtspl = wtdest_ptrtspl %>% helperfun(wfpbb = T),
    # synthpop_ptrtspl = res_synthpop_ptrtspl %>% helperfun(wfpbb=T),
    synthpop_ptrtzonly = res_synthpop_ptrtzonly %>% helperfun(wfpbb=T),
    time_synthpop_ptrtzonly = mean(synthpop_ptrtzonly_time, na.rm=T)
  )

saveRDS(reswrite, "results-clean/res_noaiptw_application.RDS")

temp <- readRDS("results-clean/res_noaiptw_application.RDS")




