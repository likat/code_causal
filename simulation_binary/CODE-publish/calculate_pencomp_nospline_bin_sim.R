## FUNCTION FOR PENCOMP
calculate.linimpute <- function(outc_model,popdat = pop,
                                       y.varname,
                                       trt.varname){
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
    dfcopy_trtgrp <- boot_all[boot_all[,trt.varname]==trtgrp,]
    ## Fit imputation model
    mod_impute <- glm(outc_model, data=dfcopy_trtgrp, family = "binomial")
    modcoeff[trtgrp+1] <- mod_impute$coefficients[1]
    ## save imputations for the other treatment group
    impute_potenout <- rbinom(nrow(samp),size=1,predict(mod_impute,newdata = samp, type = "response"))
    potenout[trtgrp+1] <- mean(impute_potenout)
    if(trtgrp == 0){
      imputed_ctrl<- impute_potenout
      imputed_ctrl_pop <- rbinom(nrow(popdat), size= 1, predict(mod_impute, newdata=popdat, type = "response"))
    }
    if(trtgrp == 1){
      imputed_trt <- impute_potenout
      imputed_trt_pop <- rbinom(nrow(popdat),size=1,predict(mod_impute, newdata=popdat, type = "response"))
    }
    
  }
  imputed_ctrl[samp[,trt.varname]==0] <- samp$Ysamp[samp[,trt.varname]==0]
  imputed_trt[samp[,trt.varname]==1] <- samp$Ysamp[samp[,trt.varname]==1]
  
  imputed_ctrl_pop[pop$I==1][samp[,trt.varname]==0] <- samp$Ysamp[samp[,trt.varname]==0]
  imputed_trt_pop[pop$I==1][samp[,trt.varname]==1] <- samp$Ysamp[samp[,trt.varname]==1]
  
  pate_gen <- 
    mean(imputed_trt_pop - imputed_ctrl_pop)
  varpategen <- var(imputed_trt_pop - imputed_ctrl_pop)/npop
    # (4) return imputed values 
  return(list(imputed_ctrl, imputed_trt, 
              pate_gen=pate_gen,
              varpategen = varpategen
  ))
  # }
}




