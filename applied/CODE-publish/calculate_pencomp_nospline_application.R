## FUNCTION FOR PENCOMP
calculate.linimpute <- function(outc_model,
                                       samp = dat,
                                       y.varname,
                                       trt.varname){
  # (2) Fit propensity score model 
  imputed_ctrl <- rep(NA,nrow(samp))
  imputed_trt <- rep(NA,nrow(samp))
  modcoeff <- c(NA,NA)
  potenout <- c(NA,NA)

  for(trtgrp in 0:1){
    ## separate fitted probabilities by treatment group
    dfcopy_trtgrp <- NULL
    dfcopy_trtgrp <- boot_all[boot_all[,trt.varname]==trtgrp,]
    ## Fit imputation model
    mod_impute <- glm(outc_model, data=dfcopy_trtgrp, family = binomial())
    modcoeff[trtgrp+1] <- mod_impute$coefficients[1]
    ## save imputations for the other treatment group
    impute_potenout <- rbinom(nrow(samp), size = 1, predict(mod_impute,newdata = samp, type = "response"))
    potenout[trtgrp+1] <- mean(impute_potenout)
    if(trtgrp == 0){
      imputed_ctrl<- impute_potenout
    }
    if(trtgrp == 1){
      imputed_trt <- impute_potenout
    }
    
  }
  imputed_ctrl[samp[,trt.varname]==0] <- samp[samp[,trt.varname]==0, y.varname]
  imputed_trt[samp[,trt.varname]==1] <- samp[samp[,trt.varname]==1, y.varname]
  

    # (4) return imputed values 
    return(list(imputed_ctrl, imputed_trt
    ))
    
  # }
}




