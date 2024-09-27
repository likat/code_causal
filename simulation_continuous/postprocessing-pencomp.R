iter = 200
N=66000
n=6600
resnames <- c("est","withinvar","btwnvar","mivar","ci_tlength","ci_emplength","cicov_t","cicov_emp") # must match names from rubin function

get.summary <- function(x, truex, typewfpbb = F){
  sampsize <- readRDS(paste0("results-clean/sampsize_",type_outcome, type_overlap,".RDS"))
  fpc <- (N - sampsize) / (N-1)
  if(typewfpbb == T){
    bias <- mean(x$est - truex)
    rmse <- sqrt(mean((x$est - truex)^2))
    cilength <- mean(x$ci_tlength)
    cicov_t <- mean(x$cicov_t)
    
  }else{
    bias <- mean(x$est - truex)
    rmse <- sqrt(mean((x$est - truex)^2))
    x$mivar <- x$mivar*fpc
    tdf <- round((200-1)*(1+x$withinvar/((200+1)*x$btwnvar)^2))
    cilower_t <- x$est - qt(0.975,tdf)*sqrt(x$mivar)
    ciupper_t <- x$est + qt(0.975,tdf)*sqrt(x$mivar)
    # cilength <- mean(x$ci_tlength)
    cilength <- mean(ciupper_t - cilower_t)
    cicov_t <- mean(truex >= cilower_t & truex <= ciupper_t)
    # cicov_t <- mean(x$cicov_t)
    
  }
  
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
	  res <- readRDS(paste0("results-raw/results_tingting_noaiptw", type_outcome, type_overlap,i,".RDS"))
    ate_nospline <- Map(c, ate_nospline, res$ate_nospline)
    ate_pencomp <- Map(c,ate_pencomp, res$ate_pencomp)
    wtdest_pencomp <- Map(c, wtdest_pencomp, res$wtdest_pencomp)
    wtdest_pencompwfpbb <- Map(c, wtdest_pencompwfpbb, res$wtdest_pencompwfpbb)
    wtdest_nospline <- Map(c, wtdest_nospline, res$wtdest_nospline)
    synthpop_pencompwfpbb <- Map(c, synthpop_pencompwfpbb, res$synthpop_pencompwfpbb)
    twostage_pencomp <- Map(c, twostage_pencomp, res$twostage_pencomp)
    twostage_pencompwfpbb <- Map(c, twostage_pencompwfpbb, res$twostage_pencompwfpbb)
    twostage_nospline<- Map(c, twostage_nospline, res$twostage_nospline)
    
    coef_pencomp <- Map(rbind, coef_pencomp,res$coef_pencomp)
    coef_pencompwfpbb <- Map(rbind, coef_pencompwfpbb,res$coef_pencompwfpbb )
    
  }
  temp <- readRDS(paste0("results-raw/results_tingting_noaiptw", type_outcome, type_overlap,"1.RDS"))
  names(ate_nospline) <- resnames
  names(ate_pencomp) <- resnames
  names(wtdest_pencomp) <- resnames
  names(wtdest_pencompwfpbb) <- resnames
  names(synthpop_pencompwfpbb) <- resnames
  names(wtdest_nospline) <- resnames
  names(twostage_pencomp) <- resnames
  names(twostage_pencompwfpbb) <- resnames
  names(twostage_nospline) <- resnames
      twostage_pencompwfpbb$cilower <- twostage_pencompwfpbb$est - (1+1/200)*sqrt(twostage_pencompwfpbb$btwnvar)*qt(0.975, 200-1)
  twostage_pencompwfpbb$ciupper <- twostage_pencompwfpbb$est + (1+1/200)*sqrt(twostage_pencompwfpbb$btwnvar)*qt(0.975, 200-1)
  twostage_pencompwfpbb$ci_tlength <- twostage_pencompwfpbb$ciupper-twostage_pencompwfpbb$cilower
  twostage_pencompwfpbb$cicov_t <- (temp$patetrue <= twostage_pencompwfpbb$ciupper & temp$patetrue >= twostage_pencompwfpbb$cilower)
  
   wtdest_pencompwfpbb$cilower <- wtdest_pencompwfpbb$est - (1+1/200)*sqrt(wtdest_pencompwfpbb$btwnvar)*qt(0.975, 200-1)
  wtdest_pencompwfpbb$ciupper <- wtdest_pencompwfpbb$est + (1+1/200)*sqrt(wtdest_pencompwfpbb$btwnvar)*qt(0.975, 200-1)
  wtdest_pencompwfpbb$ci_tlength <- wtdest_pencompwfpbb$ciupper-wtdest_pencompwfpbb$cilower
  wtdest_pencompwfpbb$cicov_t <- (temp$patetrue <= wtdest_pencompwfpbb$ciupper & temp$patetrue >= wtdest_pencompwfpbb$cilower)
  
  synthpop_pencompwfpbb$cilower <- synthpop_pencompwfpbb$est - (1+1/200)*sqrt(synthpop_pencompwfpbb$btwnvar)*qt(0.975, 200-1)
  synthpop_pencompwfpbb$ciupper <- synthpop_pencompwfpbb$est + (1+1/200)*sqrt(synthpop_pencompwfpbb$btwnvar)*qt(0.975, 200-1)
  synthpop_pencompwfpbb$ci_tlength <- synthpop_pencompwfpbb$ciupper-synthpop_pencompwfpbb$cilower
  synthpop_pencompwfpbb$cicov_t <- (temp$patetrue <= synthpop_pencompwfpbb$ciupper & temp$patetrue >= synthpop_pencompwfpbb$cilower)
  reswrite <- list(
    patetrue = temp$patetrue,
    summary_ate_pencomp = get.summary(ate_pencomp, truex = temp$patetrue),
    summary_ate_nospline=get.summary(ate_nospline, truex = temp$patetrue),
    summary_wtdest_pencomp=get.summary(wtdest_pencomp, truex = temp$patetrue),
    summary_wtdest_pencompwfpbb = get.summary(wtdest_pencompwfpbb, truex = temp$patetrue, typewfpbb = T),
    summary_wtdest_nospline = get.summary(wtdest_nospline, truex = temp$patetrue),
    summary_synthpop_pencompwfpbb = get.summary(synthpop_pencompwfpbb, truex = temp$patetrue, typewfpbb = T),
    summary_twostage_pencomp=get.summary(twostage_pencomp, truex = temp$patetrue),
    summary_twostage_pencompwfpbb = get.summary(twostage_pencompwfpbb, truex = temp$patetrue, typewfpbb = T),
    summary_twostage_nospline = get.summary(twostage_nospline, truex = temp$patetrue),
    ate_pencomp=ate_pencomp,
    ate_nospline=ate_nospline,
    wtdest_pencomp=wtdest_pencomp,
    wtdest_pencompwfpbb = wtdest_pencompwfpbb,
    wtdest_nospline = wtdest_nospline,
    synthpop_pencompwfpbb = synthpop_pencompwfpbb,
    twostage_pencomp=twostage_pencomp,
    twostage_pencompwfpbb = twostage_pencompwfpbb,
    twostage_nospline = twostage_nospline,
    coef_pencomp = coef_pencomp,
    coef_pencompwfpbb= coef_pencompwfpbb
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
saveRDS(getres(type_outcome = outcomecurrent, type_overlap = overlapcurrent), file = paste0("results-clean/res_",outcomecurrent,overlapcurrent,".RDS"))
    
  }
}
res <- readRDS(paste0("results-clean/res_",outcomecurrent,overlapcurrent,".RDS"))
