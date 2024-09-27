
###final.mod-prediction model
###baselineVar-baseline covariates
###numCarlo-number of monte carlo runs use in G computation
###here I use A0-to denote the first treatment indicator
###treat-treatment indicator (for the treatment that you are interested in estimating the associated potential outcomes)

###treat.varname-variable name denoting first treatment
###in the simulation studies, A0-denotes treatment at first time point

gcomputeFunc=function(final.mod, data,
                      des=mod3_des, 
                      baselineVar, numCarlo, treat.varname, outcome.varname, treat){
  
  ############################################################### 
  ### generate baseline covariates from empirical distributions########
  ## simply bootstrap index to preserve relationships between the variables
  baselineSim=NULL
  baseline_des <- model.matrix(des)
  # draw_g=NULL
  bootstrapind <- sample(1:nrow(baseline_des), size=numCarlo, replace=T)
  baselineSim <- as.data.frame(baseline_des[bootstrapind,])
  # for(ind in 1:length(baselineVar)){
  #   Xcov_current <- baselineVar[ind]
  #   Xcov_prop <- table(Xcov_current)/nrow(Xcov_current)
  #   if(colnames(baselineVar)[ind] == "hhsize_scaled"){
  #     draw_gtemp <- 
  #       t(rmultinom(n=numCarlo,size = 1, prob = Xcov_prop)*as.numeric(names(Xcov_prop))) %>% 
  #       rowSums
  #     draw_g <- cbind(draw_g,draw_gtemp)
  #   }else{
  #   draw_g <- cbind(draw_g, t(rmultinom(n=numCarlo,size = 1, prob = Xcov_prop[-1])))
  #   }
  #   # draw_g=cbind(draw_g, rnorm(numCarlo, mean=mean(data[,baselineVar[ind]]), sd=sd(data[,baselineVar[ind]])) )
  # }
  
  # colnames(draw_g)=paste0(covariateXnames)
  # baselineSim=cbind(baselineSim, draw_g)
  # baselineSim <- as.data.frame(draw_g)

  ####squared terms ###############
  # baselineSim[ , "age35-65:wiceligYes"] = baselineSim[,"35-65"]*baselineSim[,"Yes"]
  # baselineSim[ , "age65+:wiceligYes"] = baselineSim[,"65+"]*baselineSim[,"Yes"]
  # baselineSim[ , "L3_sq"] = baselineSim[, "L3"]^2
  # baselineSim[, "L1L2"] = baselineSim[, "L1"] * baselineSim[, "L2"]
  # includeDesign <- baselineSim %>% as.matrix
  includeDesign <- model.matrix(final.mod,baselineSim) %>% as.matrix
  
  
  ###### generate the final outcome of interest Y 
  ## LEFT OFF HERE 05/23
  # includeVar = NULL
  # includeVar = names(final.mod$coef)[-c(1)]
  # includeDesign=NULL
  # for(ind in 1:length(includeVar)){
  #   includeDesign=cbind(includeDesign, baselineSim[, includeVar[ind]])
  # }
  
  
  # Y_mean=(cbind(rep(1,numCarlo), includeDesign)%*% final.mod$coefficients) 
  Y_mean=(includeDesign%*% final.mod$coefficients) 
  
  baselineSim$A0_g <- treat
  
  names(baselineSim)[which(names(baselineSim) == "A0_g")]=treat.varname
  # Y_g=rnorm(numCarlo, mean=Y_mean, sd=summary(final.mod)$sigma) 
  Y_g= rbinom(numCarlo, size = 1, prob = expit(Y_mean))
  
  baselineSim = data.frame(baselineSim, Y=Y_g)
  names(baselineSim)[which(names(baselineSim) == "Y")] = outcome.varname
  return( baselineSim)
  
  
}




