## FUNCTION FOR PENCOMP
calculate.pencomp.nospline <- function(outc_model,
                                       popdat=pop,
                                       bootind0 = bootind_cnt,
                                       bootind1 = bootind_trt,
                                       sampind0 = ind_cnt,
                                       sampind1 = ind_trt){
  # pop_nonsel <- popdat[popdat$I==0,]
  # (2) Fit propensity score model 
  imputed_ctrl <- rep(NA,nrow(samp))
  imputed_trt <- rep(NA,nrow(samp))
  modcoeff <- c(NA,NA)
  potenout <- c(NA,NA)
  imputed_ctrl_pop <- rep(NA,nrow(popdat))
  imputed_trt_pop <- rep(NA,nrow(popdat))
  
  for(trtgrp in 0:1){
    ## separate fitted probabilities by treatment group
    dfcopy_trtgrp <- NULL
    dfcopy_trtgrp <- boot_all[boot_all$A==trtgrp,]
    ## Fit imputation model
    mod_impute <- lm(outc_model, data=dfcopy_trtgrp)
    modcoeff[trtgrp+1] <- mod_impute$coefficients[1]
    ## save imputations for the other treatment group
    impute_potenout <- predict(mod_impute,newdata = samp[samp$A!=trtgrp,])
    potenout[trtgrp+1] <- mean(impute_potenout)
    if(trtgrp == 0){
      imputed_ctrl[samp$A == 0] <- samp$Ysamp[samp$A == 0]
      imputed_ctrl[samp$A == 1] <- impute_potenout
      imputed_ctrl_pop <- predict(mod_impute, newdata=popdat)
    }
    if(trtgrp == 1){
      imputed_trt[samp$A == 1] <- samp$Ysamp[samp$A == 1]
      imputed_trt[samp$A == 0] <- impute_potenout
      imputed_trt_pop <- predict(mod_impute, newdata=popdat)
    }
    
  }
  
  imputed_ctrl_pop[pop$I==1][samp$A==0] <- samp$Ysamp[samp$A==0]
  imputed_trt_pop[pop$I==1][samp$A==1] <- samp$Ysamp[samp$A==1]
  
  # (4) return imputed values 
  # if(bootflag == 1){print("rebootstrap")}
  # else{
    poptemp <- data.frame(
      ind = pop$catecat,
      imputed0 = imputed_ctrl_pop,
      imputed1 = imputed_trt_pop
    )
    ## boot level ATE
    temptblboot <- data.frame(ind = samp$catecat,
                          y0= imputed_ctrl,
                          y1=imputed_trt)
    temptblboot <- temptblboot %>% group_by(ind) %>% summarise(cate=mean(y1-y0),
                                                       cts = n(),
                                                       catevar=var(y1-y0)/cts)
    ## generalized PATE
      pate_gen <- 
      mean(poptemp$imputed1- poptemp$imputed0)
    varpategen <- var(poptemp$imputed1-poptemp$imputed0)/npop
    ## generalized CATE
    temptbl<-poptemp %>% group_by(ind) %>% summarise(cate=mean(imputed1-imputed0),
                                                                 cts = n(),
                                                                 catevar=var(imputed1-imputed0)/cts)
    # (4) return imputed values 
    return(list(imputed_ctrl, imputed_trt, 
                cate = temptblboot$cate,
                catevar = temptblboot$catevar,
                pate_gen=pate_gen,
                varpategen = varpategen,
                cate_gen = temptbl$cate,
                varcategen = temptbl$catevar
    ))
    
  # }
}




