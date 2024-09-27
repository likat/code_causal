iter = 200
resnames <- c("est","withinvar","btwnvar","mivar","ci_tlength","ci_emplength","cicov_t","cicov_emp") # must match names from rubin function

get.summary <- function(x, truex){
  sampsize <- readRDS(paste0("results-clean/sampsize_",type_outcome, type_overlap,".RDS"))
  fpc <- (N - sampsize) / (N-1)
  
  bias <- mean(x$est - truex)
  rmse <- sqrt(mean((x$est - truex)^2))
  x$mivar <- x$mivar*fpc
  tdf <- round((200-1)*(1+x$withinvar/((200+1)*x$btwnvar)^2))
  cilower_t <- x$est - qt(0.975,tdf)*sqrt(x$mivar)
  ciupper_t <- x$est + qt(0.975,tdf)*sqrt(x$mivar)
  # cilength <- mean(x$ci_tlength)
  cilength <- mean(ciupper_t - cilower_t)
  cicov_t <- mean(truex >= cilower_t & truex <= ciupper_t)
  
    res <- list(
    bias = bias,
    rmse = rmse,
    cilength = cilength,
    cicov_t = cicov_t
    )
return(res)
}

getres <- function(type_outcome, type_overlap){
  wtdest_pencomp <- vector(mode="list", length=8)
  wtdest_pencompwfpbb <- vector(mode="list", length=8)
  wtdest_nospline <- vector(mode="list", length=8)
  twostage_pencomp <- vector(mode="list", length=8)
  twostage_pencompwfpbb <- vector(mode="list", length=8)
  twostage_nospline <- vector(mode="list", length=8)
  synthpop_pencompwfpbb <- vector(mode = "list", length = 8)
  ate_pencomp <- vector(mode="list", length=8)
  ate_pencomp_wfpbb <- vector(mode="list", length=8)
  ate_nospline <- vector(mode="list", length=8)
  coef_pencomp <- vector(mode = "list", length = 4)
  coef_pencompwfpbb <- vector(mode="list", length = 4)
  
  itersum <- 0
  for(i in 1:iter){
	  res <- readRDS(paste0("results-raw/results_tingting_linimpute", type_outcome, type_overlap,i,".RDS"))

    ate_nospline <- Map(c, ate_nospline, res$ate_nospline)
    wtdest_nospline <- Map(c, wtdest_nospline, res$wtdest_nospline)
    twostage_nospline<- Map(c, twostage_nospline, res$twostage_nospline)
    
    # coef_pencomp <- Map(rbind, coef_pencomp,res$coef_pencomp)
    # coef_pencompwfpbb <- Map(rbind, coef_pencompwfpbb,res$coef_pencompwfpbb )
    
  }
  temp <- readRDS(paste0("results-raw/results_tingting_linimpute", type_outcome, type_overlap,"1.RDS"))
  names(ate_nospline) <- resnames
  names(wtdest_nospline) <- resnames
  names(twostage_nospline) <- resnames
 
  reswrite <- list(
    patetrue = temp$patetrue,
    est = twostage_nospline$est,
    summary_ate_nospline=get.summary(ate_nospline, truex = temp$patetrue),
    summary_wtdest_nospline = get.summary(wtdest_nospline, truex = temp$patetrue),
    summary_twostage_nospline = get.summary(twostage_nospline, truex = temp$patetrue)
    # coef_pencompwfpbb= coef_pencompwfpbb
    )
return(reswrite)
}

outcometypes <- c("nogb", "gb")
overlaptypes <- c("notshared", "shared")

for(outcome in 1:length(outcometypes)){
  for(overlap in 1:length(overlaptypes)){
    outcomecurrent <- outcometypes[outcome]
    overlapcurrent <- overlaptypes[overlap]
    getres(type_outcome = outcomecurrent, type_overlap = overlapcurrent)
saveRDS(getres(type_outcome = outcomecurrent, type_overlap = overlapcurrent), file = paste0("results-clean/res_linimpute",outcomecurrent,overlapcurrent,".RDS"))
    
  }
}
