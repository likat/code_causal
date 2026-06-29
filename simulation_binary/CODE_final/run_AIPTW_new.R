### function to calculate the doubly robust coefficients from the estimating equations 
###dataInput-the dataset
###propenVarlist-the variables in the propensity model
###treat.varname-names of the treatment 
###outcome.varname-name of the outcome variable
##outcomeVarList0-variables included in the outcome model for Y0
##outcomeVarList1-varialbes included in the outcome model for Y1
##original-1 for calculating treatment effect estimate on the original dataset, instead of bootstrap samples (for calculating standard error)
###numCarlo-number of monte carlo runs in g computation and estimating the augmentation term in AIPTW
###baseList-variables at baseline, draw from baseline empirical distributions (in the simulations there are 3 baseline covariates for one time point treatment)

helper.aipw <- function(treatIndicator, selwts, ptrt, y, ypred){
  # res <- 1/sum(selwts)*sum(1/ptrt *selwts* (treatIndicator*(y-ypred) + ptrt*ypred))
  res <- 1/sum(selwts)*sum(selwts* (ypred + treatIndicator/ptrt *(y-ypred)))
  return(res)
}


AIPW = 
  function(
    dataInput=samp, 
    propen.mod = propensity_mod_aiptw,
    propen.mod.mis = propensity_mod_mis,
    modely1 = outcome_mod,
    modely0 = outcome_mod,
    modely1_correct = NULL,
    modely0_correct = NULL,
    benchmark = F,
    outcomefamily = list(c("gaussian", "binomial"))) { 
  
    t0 <- Sys.time()
    if(benchmark){
      
      ## fit propensity model
      modeltrt=glm(propen.mod, data=dataInput, family="binomial") 
      modeltrt.mis <- glm(propen.mod.mis, data=dataInput, family="binomial") 

      ## treatment indicator
      treatInd <- dataInput$A
      
      ## predicted propensity of receiving the assigned treatment
      pred_ptrt=treatInd * modeltrt$fitted.values + (1-treatInd) * (1-modeltrt$fitted.values)
      pred_ptrt_mis = treatInd * modeltrt.mis$fitted.values + (1-treatInd) * (1-modeltrt.mis$fitted.values)
      dat1 = dataInput[dataInput$A==1,]
      dat0 = dataInput[dataInput$A==0,]

      ## unweighted outcome model
      modelout_1_unwtd= glm(modely1, data=dataInput[dataInput$A==1,], family = outcomefamily) 
      modelout_0_unwtd= glm(modely0, data=dataInput[dataInput$A==0,], family = outcomefamily) 
      modelout_1_correct= glm(modely1_correct, data=dataInput[dataInput$A==1,], family = outcomefamily) 
      modelout_0_correct= glm(modely0_correct, data=dataInput[dataInput$A==0,], family = outcomefamily) 
      
      ## predicted Y
      ypred1_unwtd <- predict(modelout_1_unwtd, newdata = dataInput)
      ypred0_unwtd <- predict(modelout_0_unwtd, newdata = dataInput)
      ypred1_correct <- predict(modelout_1_correct, newdata = dataInput)
      ypred0_correct <- predict(modelout_0_correct, newdata = dataInput)
      
      ## AIPW estimator
      
      mu1.correct <- helper.aipw(treatIndicator = dataInput$A==1, ptrt = pred_ptrt,
                         y = dataInput$Ysamp, ypred = ypred1_correct,
                         selwts = dataInput$wts)
      mu0.correct <- helper.aipw(treatIndicator = dataInput$A==0, ptrt = pred_ptrt,
                         y = dataInput$Ysamp, ypred = ypred0_correct,
                         selwts = dataInput$wts)
      mu1.ptrtwrong <- helper.aipw(treatIndicator = dataInput$A==1, ptrt = pred_ptrt_mis,
                                y = dataInput$Ysamp, ypred = ypred1_correct,
                                selwts = dataInput$wts)
      mu0.ptrtwrong <- helper.aipw(treatIndicator = dataInput$A==0, ptrt = pred_ptrt_mis,
                                y = dataInput$Ysamp, ypred = ypred0_correct,
                                selwts = dataInput$wts)
      mu1.outwrong <- helper.aipw(treatIndicator = dataInput$A==1, ptrt = pred_ptrt,
                                   y = dataInput$Ysamp, ypred = ypred1_unwtd,
                                   selwts = dataInput$wts)
      mu0.outwrong <- helper.aipw(treatIndicator = dataInput$A==0, ptrt = pred_ptrt,
                                   y = dataInput$Ysamp, ypred = ypred0_unwtd,
                                   selwts = dataInput$wts)
      t1 <- Sys.time()
      res <- list(pate_correct = mu1.correct - mu0.correct, 
                  pate_ptrtwrong = mu1.ptrtwrong - mu0.ptrtwrong, 
                  pate_outwrong = mu1.outwrong - mu0.outwrong,
                  # coef1 = modelout_1$coefficients, coef0 = modelout_0$coefficients,
                  # coef1_unwtd = modelout_1_unwtd$coefficients, coef0_unwtd = modelout_0_unwtd$coefficients,
                  # trtcoef_unwtd = modeltrt.mis$coefficients,
                  # trtcoef_wtd = modeltrt.miswtd$coefficients,
                  time = t1 - t0)
    }else{ 
      ## fit propensity model
      modeltrt=glm(propen.mod, data=dataInput, family="binomial") 
      modeltrt.mis <- glm(propen.mod.mis, data=dataInput, family="binomial") 
      svydes <- svydesign(ids = ~0, weights =~wts, data = dataInput)
      modeltrt.miswtd <- svyglm(propen.mod.mis, design = svydes, family = "binomial")
      modeltrt.wtd <- svyglm(propen.mod, design = svydes, family = "binomial")
      
      ## treatment indicator
      treatInd <- dataInput$A
      
      ## predicted propensity of receiving the assigned treatment
      pred_ptrt=treatInd * modeltrt$fitted.values + (1-treatInd) * (1-modeltrt$fitted.values)
      pred_ptrt_mis = treatInd * modeltrt.mis$fitted.values + (1-treatInd) * (1-modeltrt.mis$fitted.values)
      pred_ptrt_miswtd = treatInd * modeltrt.miswtd$fitted.values + (1-treatInd) * (1-modeltrt.miswtd$fitted.values)
      pred_ptrt_wtd = treatInd * modeltrt.wtd$fitted.values + (1-treatInd) * (1-modeltrt.wtd$fitted.values)
      
      ## establish survey design for weighted outcome model
      svydes0 <- svydesign(ids = ~0, weights =~wts, data = dataInput[dataInput$A==0,])
      svydes1 <- svydesign(ids = ~0, weights =~wts, data = dataInput[dataInput$A==1,])
      
      dat1 = dataInput[dataInput$A==1,]
      dat0 = dataInput[dataInput$A==0,]
      modelout_1= svyglm(modely1, design = svydes1, family = outcomefamily) 
      modelout_0= svyglm(modely0, design = svydes0, family = outcomefamily) 
      
      ## unweighted outcome model
      modelout_1_unwtd= glm(modely1, data=dataInput[dataInput$A==1,], family = outcomefamily) 
      modelout_0_unwtd= glm(modely0, data=dataInput[dataInput$A==0,], family = outcomefamily) 
      
      ## predicted Y
      ypred1 <- predict(modelout_1, newdata = dataInput)
      ypred0 <- predict(modelout_0, newdata = dataInput)
      ypred1_unwtd <- predict(modelout_1_unwtd, newdata = dataInput)
      ypred0_unwtd <- predict(modelout_0_unwtd, newdata = dataInput)
      
      ## AIPW estimator
      mu1.allwrong <- helper.aipw(treatIndicator = dataInput$A==1, ptrt = pred_ptrt_mis,
                                  y = dataInput$Ysamp, ypred = ypred1_unwtd,
                                  selwts = dataInput$wts)
      mu0.allwrong <- helper.aipw(treatIndicator = dataInput$A==0, ptrt = pred_ptrt_mis,
                                  y = dataInput$Ysamp, ypred = ypred0_unwtd,
                                  selwts = dataInput$wts)
      
      mu1.unwtd.correcttrt <- helper.aipw(treatIndicator = dataInput$A==1, ptrt = pred_ptrt,
                         y = dataInput$Ysamp, ypred = ypred1_unwtd,
                         selwts = dataInput$wts)
      mu0.unwtd.correcttrt <- helper.aipw(treatIndicator = dataInput$A==0, ptrt = pred_ptrt,
                               y = dataInput$Ysamp, ypred = ypred0_unwtd,
                         selwts = dataInput$wts)
      mu1.wtd.wtdtrt <- helper.aipw(treatIndicator = dataInput$A==1, ptrt = pred_ptrt_miswtd,
                                y = dataInput$Ysamp, ypred = ypred1,
                                selwts = dataInput$wts)
      mu0.wtd.wtdtrt <- helper.aipw(treatIndicator = dataInput$A==0, ptrt = pred_ptrt_miswtd,
                                y = dataInput$Ysamp, ypred = ypred0,
                                selwts = dataInput$wts)
      mu1.wtd.mistrt <- helper.aipw(treatIndicator = dataInput$A==1, ptrt = pred_ptrt_mis,
                             y = dataInput$Ysamp, ypred = ypred1,
                             selwts = dataInput$wts)
      mu0.wtd.mistrt <- helper.aipw(treatIndicator = dataInput$A==0, ptrt = pred_ptrt_mis,
                                    y = dataInput$Ysamp, ypred = ypred0,
                             selwts = dataInput$wts)
      mu1.unwtd.wtdtrt <- helper.aipw(treatIndicator = dataInput$A==1, ptrt = pred_ptrt_miswtd,
                                y = dataInput$Ysamp, ypred = ypred1_unwtd,
                                selwts = dataInput$wts)
      mu0.unwtd.wtdtrt <- helper.aipw(treatIndicator = dataInput$A==0, ptrt = pred_ptrt_miswtd,
                                y = dataInput$Ysamp, ypred = ypred0_unwtd,
                                selwts = dataInput$wts)
      t1 <- Sys.time()
      
      res <- list(pate_unwtd_correcttrt = mu1.unwtd.correcttrt-mu0.unwtd.correcttrt, 
                  pate_wtd_wtdtrt = mu1.wtd.wtdtrt - mu0.wtd.wtdtrt, 
                  pate_wtd_mistrt = mu1.wtd.mistrt - mu0.wtd.mistrt,
                  pate_unwtd_wtdtrt = mu1.unwtd.wtdtrt-mu0.unwtd.wtdtrt, 
                  pate_allwrong = mu1.allwrong - mu0.allwrong,
                  coef1 = modelout_1$coefficients, coef0 = modelout_0$coefficients,
                  coef1_unwtd = modelout_1_unwtd$coefficients, coef0_unwtd = modelout_0_unwtd$coefficients,
                  trtcoef_unwtd = modeltrt.mis$coefficients,
                  trtcoef_wtd = modeltrt.miswtd$coefficients,
                  time = t1 - t0)
    }
    return(res)
} 

AIPW.application = 
  function(
    dataInput=samp, 
    propen.mod = propensity_mod_aiptw,
    out.name,
    trt.name,
    # propen.mod.mis = propensity_mod_mis,
    modely1 = outcome_mod,
    modely0 = outcome_mod,
    outcomefamily = list(c("gaussian", "binomial"))) { 
    
    t0 <- Sys.time()
    ## fit propensity model
    # modeltrt=glm(propen.mod, data=dataInput, family="binomial") 
    svydes <- svydesign(ids = ~0, weights =~wts, data = dataInput)
    modeltrt <- svyglm(propen.mod, design = svydes, family = "binomial")
    
    ## treatment indicator
    treatInd <- dataInput[,trt.name]
    
    ## predicted propensity of receiving the assigned treatment
    pred_ptrt=treatInd * modeltrt$fitted.values + (1-treatInd) * (1-modeltrt$fitted.values)
    svydes0 <- svydesign(ids = ~0, weights =~wts, data = dataInput[treatInd==0,])
    svydes1 <- svydesign(ids = ~0, weights =~wts, data = dataInput[treatInd==1,])
    
    dat1 = dataInput[treatInd==1,]
    dat0 = dataInput[treatInd==0,]
    modelout_1= svyglm(modely1, design = svydes1, family = outcomefamily) 
    modelout_0= svyglm(modely0, design = svydes0, family = outcomefamily) 
    
    ## unweighted outcome model
    modelout_1_unwtd= glm(modely1, data=dataInput[treatInd==1,], family = outcomefamily) 
    modelout_0_unwtd= glm(modely0, data=dataInput[treatInd==0,], family = outcomefamily) 
    
    ## predicted Y
    
    ypred1 <- predict(modelout_1, newdata = dataInput, type = "response")
    ypred0 <- predict(modelout_0, newdata = dataInput, type = "response")
    ypred1_unwtd <- predict(modelout_1_unwtd, newdata = dataInput, type = "response")
    ypred0_unwtd <- predict(modelout_0_unwtd, newdata = dataInput, type = "response")

    
    ## AIPW estimator
    mu1 <- helper.aipw(treatIndicator = (treatInd==1), ptrt = pred_ptrt,
                       y = dataInput[,out.name] %>% as.numeric, ypred = ypred1,
                       selwts = dataInput$wts)
    mu0 <- helper.aipw(treatIndicator =( treatInd==0), ptrt = pred_ptrt,y = dataInput[,out.name], ypred = ypred0,
                       selwts = dataInput$wts)
    
    mu1_unwtd <- helper.aipw(treatIndicator = (treatInd==1), ptrt = pred_ptrt,y = dataInput[,out.name], ypred = ypred1_unwtd,
                             selwts = dataInput$wts)
    mu0_unwtd <- helper.aipw(treatIndicator = (treatInd==0), ptrt = pred_ptrt,y = dataInput[,out.name], ypred = ypred0_unwtd,
                             selwts = dataInput$wts)
    t1 <- Sys.time()
    
    res <- list(pate_wtd = mu1-mu0, pate_unwtd = mu1_unwtd - mu0_unwtd,
                coef1 = modelout_1$coefficients, coef0 = modelout_0$coefficients,
                coef1_unwtd = modelout_1_unwtd$coefficients, coef0_unwtd = modelout_0_unwtd$coefficients,
                time = t1 - t0)
    return(res)
  } 

