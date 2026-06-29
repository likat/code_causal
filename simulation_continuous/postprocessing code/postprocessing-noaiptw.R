# iter = 250
numboot = 200
nknot = 15
N=66000
n=6600
resnames <- c("est","withinvar","btwnvar","mivar","ci_tlength","ci_emplength","cicov_t","cicov_emp") # must match names from rubin function

get.summary <- function(x, truex, typewfpbb = F){
  # sampsize <- readRDS(paste0("results-clean/sampsize_",type_outcome, type_overlap,".RDS"))
  sampsize <- 6600
  fpc <- (N - sampsize) / (N-1)
  if(typewfpbb == T){
    bias <- mean(x$est - truex)
    rmse <- sqrt(mean((x$est - truex)^2))
    cilength <- mean(x$ci_tlength, na.rm=T)
    cicov_t <- mean(x$cicov_t, na.rm=T)
    
  }else{
    bias <- mean(x$est - truex)
    rmse <- sqrt(mean((x$est - truex)^2))
    x$mivar <- x$mivar*fpc
    tdf <- round((numboot-1)*(1+x$withinvar/((numboot+1)*x$btwnvar)^2))
    cilower_t <- x$est - qt(0.975,tdf)*sqrt(x$mivar)
    ciupper_t <- x$est + qt(0.975,tdf)*sqrt(x$mivar)
    # cilength <- mean(x$ci_tlength)
    cilength <- mean(ciupper_t - cilower_t, na.rm=T)
    cicov_t <- mean(truex >= cilower_t & truex <= ciupper_t, na.rm=T)
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


ciunifun <- function(estvar, estprop, alpha = 0.05){
  deff <- estvar /(estprop*(1-estprop)/sampsize)
  srs_var <- estprop*(1-estprop)
  neff <- sampsize / deff
  kstar = estprop * neff
  cilower <- qbeta(alpha/2, kstar + 1, neff - kstar + 1)
  ciupper <- qbeta(1-alpha/2, kstar + 1, neff - kstar + 1)
  res <- list(cilower = cilower, ciupper = ciupper)
}

getres <- function(type_outcome, type_overlap, iter){
  wtdest_pencomp <- vector(mode="list", length=8)
  wtdest_pencompwfpbb <- vector(mode="list", length=8)
  wtdest_nospline <- vector(mode="list", length=8)
  wtdest_nospline_wrongmod <- vector(mode="list", length=8)
  
  twostage_pencomp <- vector(mode="list", length=8)
  twostage_pencompwfpbb <- vector(mode="list", length=8)
  twostage_nospline <- vector(mode="list", length=8)
  synthpop_pencompwfpbb <- vector(mode = "list", length = 8)
  ate_pencomp <- vector(mode="list", length=8)
  ate_pencomp_wfpbb <- vector(mode="list", length=8)
  ate_nospline <- vector(mode="list", length=8)
  ate_nospline_wrongmod <- vector(mode="list", length=8)
  coef_pencomp <- vector(mode = "list", length = 6)
  coef_pencompwfpbb <- vector(mode="list", length = 4)
  
  time_pencompwfpbb <- 0
  itersum <- 0
  skipnum <- 0
  for(i in 1:iter){
    
    temp <- tryCatch ( {
      res <- readRDS(paste0("../results-raw-2stage/", type_outcome, type_overlap, "/results_tingting_noaiptw", type_outcome, type_overlap,i,".RDS"))
       }, error=function(e) return(NA) )
    if(sum(is.na(temp))>0){skipnum = skipnum + 1; next}
    res <- readRDS(paste0("../results-raw-2stage/", type_outcome, type_overlap, "/results_tingting_noaiptw", type_outcome, type_overlap,i,".RDS"))
    res2 <- readRDS(paste0("../results-raw-2stage/linimpute/", type_outcome, type_overlap,"/results_tingting_noaiptw", type_outcome, type_overlap,i, ".RDS"))
    ate_nospline <- Map(c, ate_nospline, res$ate_nospline)
    ate_pencomp <- Map(c,ate_pencomp, res$ate_pencomp)
    ate_nospline_wrongmod <- Map(c, ate_nospline_wrongmod, res2$ate_nospline)
    wtdest_pencomp <- Map(c, wtdest_pencomp, res$wtdest_pencomp)
    wtdest_pencompwfpbb <- Map(c, wtdest_pencompwfpbb, res$wtdest_pencompwfpbb)
    wtdest_nospline <- Map(c, wtdest_nospline, res$wtdest_nospline)
    wtdest_nospline_wrongmod <- Map(c, wtdest_nospline_wrongmod, res2$wtdest_nospline)
    synthpop_pencompwfpbb <- Map(c, synthpop_pencompwfpbb, res$synthpop_pencompwfpbb)
    time_pencompwfpbb <- time_pencompwfpbb + res$time_pencompwfpbb / iter # in seconds 
    # twostage_pencomp <- Map(c, twostage_pencomp, res$twostage_pencomp)
    # twostage_pencompwfpbb <- Map(c, twostage_pencompwfpbb, res$twostage_pencompwfpbb)
    # twostage_nospline<- Map(c, twostage_nospline, res$twostage_nospline)
    
    # coef_pencomp <- Map(rbind, coef_pencomp,res$coef_pencomp)
    # coef_pencompwfpbb <- Map(rbind, coef_pencompwfpbb,res$coef_pencompwfpbb )
    itersum <- itersum + 1
    # print(itersum)
    
    
    if(itersum == 200){break}
  }
  temp <- readRDS(paste0("../results-raw-2stage/", type_outcome, type_overlap, "/results_tingting_noaiptw", type_outcome, type_overlap,2,".RDS"))
  names(ate_nospline) <- names(ate_nospline_wrongmod) <- resnames
  names(ate_pencomp) <- resnames
  names(wtdest_pencomp) <- resnames
  # names(wtdest_pencompwfpbb) <- resnames
  names(synthpop_pencompwfpbb) <- resnames
  names(wtdest_nospline) <- names(wtdest_nospline_wrongmod) <- resnames
  # names(twostage_pencomp) <- resnames
  # names(twostage_pencompwfpbb) <- resnames
  # names(twostage_nospline) <- resnames
  #     twostage_pencompwfpbb$cilower <- twostage_pencompwfpbb$est - (1+1/numboot)*sqrt(twostage_pencompwfpbb$btwnvar)*qt(0.975, numboot-1)
  # twostage_pencompwfpbb$ciupper <- twostage_pencompwfpbb$est + (1+1/numboot)*sqrt(twostage_pencompwfpbb$btwnvar)*qt(0.975, numboot-1)
  # twostage_pencompwfpbb$ci_tlength <- twostage_pencompwfpbb$ciupper-twostage_pencompwfpbb$cilower
  # twostage_pencompwfpbb$cicov_t <- (temp$patetrue <= twostage_pencompwfpbb$ciupper & temp$patetrue >= twostage_pencompwfpbb$cilower)
  
  #  wtdest_pencompwfpbb$cilower <- wtdest_pencompwfpbb$est - (1+1/numboot)*sqrt(wtdest_pencompwfpbb$btwnvar)*qt(0.975, numboot-1)
  # wtdest_pencompwfpbb$ciupper <- wtdest_pencompwfpbb$est + (1+1/numboot)*sqrt(wtdest_pencompwfpbb$btwnvar)*qt(0.975, numboot-1)
  # wtdest_pencompwfpbb$ci_tlength <- wtdest_pencompwfpbb$ciupper-wtdest_pencompwfpbb$cilower
  # wtdest_pencompwfpbb$cicov_t <- (temp$patetrue <= wtdest_pencompwfpbb$ciupper & temp$patetrue >= wtdest_pencompwfpbb$cilower)
  
  synthpop_pencompwfpbb$cilower <- synthpop_pencompwfpbb$est - (1+1/numboot)*sqrt(synthpop_pencompwfpbb$btwnvar)*qt(0.975, numboot-1)
  synthpop_pencompwfpbb$ciupper <- synthpop_pencompwfpbb$est + (1+1/numboot)*sqrt(synthpop_pencompwfpbb$btwnvar)*qt(0.975, numboot-1)
  synthpop_pencompwfpbb$ci_tlength <- synthpop_pencompwfpbb$ciupper-synthpop_pencompwfpbb$cilower
  synthpop_pencompwfpbb$cicov_t <- (temp$patetrue <= synthpop_pencompwfpbb$ciupper & temp$patetrue >= synthpop_pencompwfpbb$cilower)
  reswrite <- list(
    patetrue = temp$patetrue,
    summary_ate_pencomp = get.summary(ate_pencomp, truex = temp$patetrue),
    summary_ate_nospline=get.summary(ate_nospline, truex = temp$patetrue),
    summary_ate_nospline_wrongmod=get.summary(ate_nospline_wrongmod, truex = temp$patetrue),
    summary_wtdest_pencomp=get.summary(wtdest_pencomp, truex = temp$patetrue),
    summary_wtdest_nospline = get.summary(wtdest_nospline, truex = temp$patetrue),
    summary_wtdest_nospline_wrongmod = get.summary(wtdest_nospline_wrongmod, truex = temp$patetrue),
    summary_synthpop_pencompwfpbb = get.summary(x= synthpop_pencompwfpbb, truex = temp$patetrue, typewfpbb = T),
    ate_pencomp=ate_pencomp,
    ate_nospline=ate_nospline,
    ate_nospline_wrongmod = ate_nospline_wrongmod,
    wtdest_pencomp=wtdest_pencomp,
    wtdest_nospline = wtdest_nospline,
    wtdest_nospline_wrongmod = wtdest_nospline_wrongmod,
    synthpop_pencompwfpbb = synthpop_pencompwfpbb,
    coef_pencomp = coef_pencomp,
    coef_pencompwfpbb= coef_pencompwfpbb,
    gooditer = iter-skipnum,
    skipnum = skipnum,
    time_pencompwfpbb = time_pencompwfpbb
    )
return(reswrite)
}

getres.other <- function(type_outcome, type_overlap, iter,idtag=NULL){
  wtdest_pencomp <- vector(mode="list", length=8)
  wtdest_nospline <- vector(mode="list", length=8)
  synthpop_pencompwfpbb <- vector(mode = "list", length = 8)
  ate_pencomp <- vector(mode="list", length=8)
  ate_pencomp_wfpbb <- vector(mode="list", length=8)
  ate_nospline <- vector(mode="list", length=8)
  coef_pencomp <- vector(mode = "list", length = 6)
  coef_pencompwfpbb <- vector(mode="list", length = 4)
  
  time_pencompwfpbb <- 0
  itersum <- 0
  skipnum <- 0
  for(i in 1:iter){
    
    temp <- tryCatch ( {
      res <- readRDS(paste0("../results-raw-2stage/", type_outcome, type_overlap, "/results_tingting_noaiptw", type_outcome, type_overlap,i,idtag,".RDS"))
    }, error=function(e) return(NA) )
    if(sum(is.na(temp))>0){skipnum = skipnum + 1; next}
    res <- readRDS(paste0("../results-raw-2stage/", type_outcome, type_overlap, "/results_tingting_noaiptw", type_outcome, type_overlap,i,idtag,".RDS"))
    # res2 <- readRDS(paste0("../results-raw-2stage/linimpute/", type_outcome, type_overlap,"/results_tingting_noaiptw", type_outcome, type_overlap,i, ".RDS"))
    ate_nospline <- Map(c, ate_nospline, res$ate_nospline)
    ate_pencomp <- Map(c,ate_pencomp, res$ate_pencomp)
    # ate_nospline_wrongmod <- Map(c, ate_nospline_wrongmod, res2$ate_nospline)
    wtdest_pencomp <- Map(c, wtdest_pencomp, res$wtdest_pencomp)
    # wtdest_pencompwfpbb <- Map(c, wtdest_pencompwfpbb, res$wtdest_pencompwfpbb)
    wtdest_nospline <- Map(c, wtdest_nospline, res$wtdest_nospline)
    # wtdest_nospline_wrongmod <- Map(c, wtdest_nospline_wrongmod, res2$wtdest_nospline)
    synthpop_pencompwfpbb <- Map(c, synthpop_pencompwfpbb, res$synthpop_pencompwfpbb)
    time_pencompwfpbb <- time_pencompwfpbb + res$time_pencompwfpbb / iter # in seconds 
    itersum <- itersum + 1
    # print(itersum)
    
    
    if(itersum == 200){break}
  }
  temp <- readRDS(paste0("../results-raw-2stage/", type_outcome, type_overlap, "/results_tingting_noaiptw", type_outcome, type_overlap,2,".RDS"))
  names(ate_nospline) <- resnames
  names(ate_pencomp) <- resnames
  names(wtdest_pencomp) <- resnames
  # names(wtdest_pencompwfpbb) <- resnames
  names(synthpop_pencompwfpbb) <- resnames
  names(wtdest_nospline) <- resnames
  synthpop_pencompwfpbb$cilower <- synthpop_pencompwfpbb$est - (1+1/numboot)*sqrt(synthpop_pencompwfpbb$btwnvar)*qt(0.975, numboot-1)
  synthpop_pencompwfpbb$ciupper <- synthpop_pencompwfpbb$est + (1+1/numboot)*sqrt(synthpop_pencompwfpbb$btwnvar)*qt(0.975, numboot-1)
  synthpop_pencompwfpbb$ci_tlength <- synthpop_pencompwfpbb$ciupper-synthpop_pencompwfpbb$cilower
  synthpop_pencompwfpbb$cicov_t <- (temp$patetrue <= synthpop_pencompwfpbb$ciupper & temp$patetrue >= synthpop_pencompwfpbb$cilower)
  reswrite <- list(
    patetrue = temp$patetrue,
    summary_ate_pencomp = get.summary(ate_pencomp, truex = temp$patetrue),
    summary_ate_nospline=get.summary(ate_nospline, truex = temp$patetrue),
    # summary_ate_nospline_wrongmod=get.summary(ate_nospline_wrongmod, truex = temp$patetrue),
    summary_wtdest_pencomp=get.summary(wtdest_pencomp, truex = temp$patetrue),
    summary_wtdest_nospline = get.summary(wtdest_nospline, truex = temp$patetrue),
    # summary_wtdest_nospline_wrongmod = get.summary(wtdest_nospline_wrongmod, truex = temp$patetrue),
    summary_synthpop_pencompwfpbb = get.summary(x= synthpop_pencompwfpbb, truex = temp$patetrue, typewfpbb = T),
    ate_pencomp=ate_pencomp,
    ate_nospline=ate_nospline,
    # ate_nospline_wrongmod = ate_nospline_wrongmod,
    wtdest_pencomp=wtdest_pencomp,
    wtdest_nospline = wtdest_nospline,
    # wtdest_nospline_wrongmod = wtdest_nospline_wrongmod,
    synthpop_pencompwfpbb = synthpop_pencompwfpbb,
    coef_pencomp = coef_pencomp,
    coef_pencompwfpbb= coef_pencompwfpbb,
    gooditer = iter-skipnum,
    skipnum = skipnum,
    time_pencompwfpbb = time_pencompwfpbb
  )
  return(reswrite)
}
restemp <- getres.other(type_outcome = "gb", type_overlap = "shared", iter = 200,idtag = "wrongptrt_rightoutcome")
saveRDS(restemp, file = paste0("./res_gbshared_wrongptrt_rightoutcome.RDS"))
outcometypes <- c("nogb", "gb")
overlaptypes <- c("notshared", "shared")
refdf <- data.frame(
  outcometypes = rep(outcometypes, 2),
  overlaptypes = rep(overlaptypes, each = 2),
  numiter = c(200)
)
for(j in 1:nrow(refdf)){
# for(j in 1:1){
  
  outcomecurrent <- refdf$outcometypes[j]
  overlapcurrent <- refdf$overlaptypes[j]
  restemp <- getres(type_outcome = outcomecurrent, type_overlap = overlapcurrent, iter = refdf$numiter[1])
  saveRDS(restemp, file = paste0("./res_",outcomecurrent,overlapcurrent,".RDS"))
  print(j)
}
# outcomecurrent = 'gb'
# overlapcurrent = "notshared"
# temp <- getres(type_outcome = "gb", type_overlap = 'shared', iter = 100)
# saveRDS(temp, file = paste0("../../results-clean/res_",outcomecurrent,overlapcurrent,"nknot",nknot,".RDS"))
# saveRDS(temp, file = paste0("../../results-clean/res_",outcomecurrent,overlapcurrent,".RDS"))
