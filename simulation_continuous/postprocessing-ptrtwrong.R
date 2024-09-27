iter = 200
N = 66000
n = 6600
require(dplyr)
resnames <- c("est","withinvar","btwnvar","mivar","ci_tlength","ci_emplength","cicov_t","cicov_emp") # must match names from rubin function

get.summary <- function(x, truex, wfpbbstatus = F, L=NULL){
    bias <- mean(x$est - truex)
    rmse <- sqrt(mean((x$est - truex)^2))
    cilength = NA
    cicov_t = NA
    if(wfpbbstatus==T){
      cilower <- x$est - sqrt((1+1/L)*x$btwnvar)*qt(0.975, L-1)
      ciupper = x$est + sqrt((1+1/L)*x$btwnvar)*qt(0.975, L-1)
      cilength <- mean(ciupper - cilower)
      cicov_t <- mean(truex <= ciupper & truex >= cilower)
      
    }else{
    cilength <- mean(x$ci_tlength)
    cicov_t <- mean(x$cicov_t)
    }
    
    res <- list(
    bias = bias,
    rmse = rmse,
    cilength = cilength,
    cicov_t = cicov_t
    )
return(res)
}

get.wfpbb.summary <- function(x, truex, L=50){
  res <- data.frame(cilower = NA,
                    ciupper = NA, 
                    ci_tlength = NA,
                    cicov_t = NA)
  res$cilower <- x$est - sqrt((1+1/L)*x$btwnvar)*qt(0.975, L-1)
  res$ciupper <- x$est + sqrt((1+1/L)*x$btwnvar)*qt(0.975, L-1)
  res$ci_tlength <- res$ciupper - res$cilower
  res$cicov_t <- (truex <= res$ciupper & truex >= res$cilower)

  return(res)
}

getres <- function(type_outcome, type_overlap, type_outmod){
  wtdest_noz <- vector(mode="list", length=8)
  wtdest_correct <- vector(mode="list", length=8)
  wtdest_pselspl <- vector(mode="list", length=8)
  
  synthpop_noz <- vector(mode="list", length = 8)
  synthpop_correct<- vector(mode="list", length = 8)
  synthpop_pselspl<- vector(mode="list", length = 8)
  
  twostage_noz <- vector(mode="list", length=8)
  twostage_correct<- vector(mode="list", length=8)
  twostage_pselspl<- vector(mode="list", length=8)
  
  itersum <- 0
  skipnum = 0
  for(i in 1:iter){
    temp <- tryCatch ( {
      res <- readRDS(paste0("results-raw-wrongptrt/results_wrongptrt_",type_outmod, type_outcome, type_overlap,i,".RDS"))
    }, error=function(e) return(NA) )
    if(sum(is.na(temp))>0){skipnum = skipnum + 1; next}else{
      
    
    
	  res <- readRDS(paste0("results-raw-wrongptrt/results_wrongptrt_",type_outmod, type_outcome, type_overlap,i,".RDS"))

	  wtdest_noz <- Map(c, wtdest_noz, res$wtdest_noz)
	  wtdest_correct <- Map(c, wtdest_correct, res$wtdest_correct)
	  wtdest_pselspl <- Map(c, wtdest_pselspl, res$wtdest_pselspl)
	  
	  twostage_noz <- Map(c, twostage_noz, res$twostage_noz)
	  twostage_correct <- Map(c, twostage_correct, res$twostage_correct)
	  twostage_pselspl <- Map(c, twostage_pselspl, res$twostage_pselspl)
	  
	  synthpop_noz <- Map(c, synthpop_noz, res$synthpop_noz)
	  synthpop_correct <- Map(c, synthpop_correct, res$synthpop_correct)
	  synthpop_pselspl <- Map(c, synthpop_pselspl, res$synthpop_pselspl)
	  itersum <- itersum + 1
	  if(itersum > 100){break}
    # coef_pencomp <- Map(rbind, coef_pencomp,res$coef_pencomp)
    # coef_pencompwfpbb <- Map(rbind, coef_pencompwfpbb,res$coef_pencompwfpbb )
    }
  }
  temp <- readRDS(paste0("results-raw/results_wrongptrt_",type_outmod, type_outcome, type_overlap,"1.RDS"))
  names(wtdest_noz) <- resnames
  names(wtdest_correct) <- resnames
  names(wtdest_pselspl) <- resnames
  
  names(synthpop_noz) <- resnames
  names(synthpop_correct) <- resnames
  names(synthpop_pselspl) <- resnames
  
  names(twostage_noz) <- resnames
  names(twostage_correct) <- resnames
  names(twostage_pselspl) <- resnames
  
  # twostage_pencompwfpbb$cilower <- twostage_pencompwfpbb$est - (1+1/200)*sqrt(twostage_pencompwfpbb$btwnvar)*qt(0.975, 200-1)
  # twostage_pencompwfpbb$ciupper <- twostage_pencompwfpbb$est + (1+1/200)*sqrt(twostage_pencompwfpbb$btwnvar)*qt(0.975, 200-1)
  # twostage_pencompwfpbb$ci_tlength <- twostage_pencompwfpbb$ciupper-twostage_pencompwfpbb$cilower
  # twostage_pencompwfpbb$cicov_t <- (temp$patetrue <= twostage_pencompwfpbb$ciupper & temp$patetrue >= twostage_pencompwfpbb$cilower)
  # 
  #  wtdest_pencompwfpbb$cilower <- wtdest_pencompwfpbb$est - (1+1/200)*sqrt(wtdest_pencompwfpbb$btwnvar)*qt(0.975, 200-1)
  # wtdest_pencompwfpbb$ciupper <- wtdest_pencompwfpbb$est + (1+1/200)*sqrt(wtdest_pencompwfpbb$btwnvar)*qt(0.975, 200-1)
  # wtdest_pencompwfpbb$ci_tlength <- wtdest_pencompwfpbb$ciupper-wtdest_pencompwfpbb$cilower
  # wtdest_pencompwfpbb$cicov_t <- (temp$patetrue <= wtdest_pencompwfpbb$ciupper & temp$patetrue >= wtdest_pencompwfpbb$cilower)
  # 
  # synthpop_pencompwfpbb$cilower <- synthpop_pencompwfpbb$est - (1+1/200)*sqrt(synthpop_pencompwfpbb$btwnvar)*qt(0.975, 200-1)
  # synthpop_pencompwfpbb$ciupper <- synthpop_pencompwfpbb$est + (1+1/200)*sqrt(synthpop_pencompwfpbb$btwnvar)*qt(0.975, 200-1)
  # synthpop_pencompwfpbb$ci_tlength <- synthpop_pencompwfpbb$ciupper-synthpop_pencompwfpbb$cilower
  # synthpop_pencompwfpbb$cicov_t <- (temp$patetrue <= synthpop_pencompwfpbb$ciupper & temp$patetrue >= synthpop_pencompwfpbb$cilower)
  reswrite <- list(
    patetrue = temp$patetrue,
    summary_wtdest_noz = wtdest_noz %>% get.summary(truex = temp$patetrue, wfpbbstatus = T, L = 50),
    summary_wtdest_correct = wtdest_correct%>% get.summary(truex = temp$patetrue, wfpbbstatus = T, L = 50),
    summary_wtdest_pselspl = wtdest_pselspl%>% get.summary(truex = temp$patetrue, wfpbbstatus = T, L = 50),
    summary_synthpop_noz = synthpop_noz%>% get.summary(truex = temp$patetrue, wfpbbstatus = T, L = 50),
    summary_synthpop_correct = synthpop_correct%>% get.summary(truex = temp$patetrue, wfpbbstatus = T, L = 50),
    summary_synthpop_pselspl = synthpop_pselspl%>% get.summary(truex = temp$patetrue, wfpbbstatus = T, L = 50),
    summary_twostage_noz = twostage_noz%>% get.summary(truex = temp$patetrue, wfpbbstatus = T, L = 50),
    summary_twostage_correct= twostage_correct%>% get.summary(truex = temp$patetrue, wfpbbstatus = T, L = 50),
    summary_twostage_pselspl= twostage_pselspl%>% get.summary(truex = temp$patetrue, wfpbbstatus = T, L = 50)
    # coef_pencomp = coef_pencomp,
    # coef_pencompwfpbb= coef_pencompwfpbb
    )
return(reswrite)
}

outcometypes <- c("nogb")
overlaptypes <- c("shared")
outmodtypes <- c("x1x2", "x2z")
for(outcome in 1:length(outcometypes)){
  for(overlap in 1:length(overlaptypes)){
    for(outmod in 1:length(outmodtypes)){
    outcomecurrent <- outcometypes[outcome]
    overlapcurrent <- overlaptypes[overlap]
    outmodcurrent <- outmodtypes[outmod]
    # getres(type_outcome = outcomecurrent, type_overlap = overlapcurrent, type_outmod = outmodcurrent)
saveRDS(getres(type_outcome = outcomecurrent, type_overlap = overlapcurrent,type_outmod =outmodcurrent), file = paste0("results-clean/res_wrongptrt_",outmodcurrent,outcomecurrent,overlapcurrent,".RDS"))
    }
  }
}
temp_noz <- readRDS("results-clean/res_wrongptrt_x1x2gbshared.RDS")
temp_withz <- readRDS("results-clean/res_wrongptrt_x2zgbshared.RDS")
temp_noz$summary_synthpop_noz
temp_noz$summary_synthpop_pselspl
temp_noz$summary_synthpop_correct

temp_withz$summary_synthpop_noz
temp_withz$summary_synthpop_pselspl
temp_withz$summary_synthpop_correct

