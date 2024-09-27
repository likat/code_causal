## FUNCTION FOR PENCOMP
calculate.pencomp <- function(propen_model, # string for propensity model, e.g. "A~X+Z"
                              outc_model,
                              samp = dat,
                              y.varname,    # name of the outcome variable   
                              trt.varname,
                              num.knot=5){
  tryCatch ( {
  knotlocations <- vector(mode="list", length = 2)
      # (2) Fit propensity score model 
    ## What if we used the true values for propensity model and outcome model?
      ### it's NOT an issue with the propensity scores;
     ### using the true ptrt biases PATE equally
      mod_trt <- glm(propen_model, data = boot_all, family="binomial")
      phat <- predict(mod_trt, newdata=boot_all,type="link")
      phat_samp <- predict(mod_trt, newdata = samp, type = "link")
      boot_all$phat <- phat
      
      # hist(phat_pop - popdat$ptrt)
      
      ## initialize empty list for imputation storage
      impute_potenout <- vector(mode="list", length=2)
      boot_all$impute <- 0
      imputed_ctrl <- rep(NA,nrow(samp))
      imputed_trt <- rep(NA,nrow(samp))
      estcoef <- vector(mode="list", length=2)
      
      for(trtgrp in 0:1){
        ## separate fitted probabilities by treatment group
        phat_trtgrp <- boot_all$phat[boot_all[,trt.varname]==trtgrp]
        dfcopy_trtgrp <- boot_all[boot_all[,trt.varname]==trtgrp,]
        
        # (3) Fit outcome model with splines
        ## need sensitivity analysis for how many knots
        ## currently borrowing code from zhou 2019
        space=(max(phat_trtgrp)-min(phat_trtgrp))/(num.knot+1)
        knots=(min(phat_trtgrp)+space*(1:num.knot))
        
        ## Construct C2 matrix, truncated linear basis splines
        linearB=NULL
        linearB =outer(phat_trtgrp, knots, "-") ## constructing C2 matrix
        linearB =linearB * (linearB > 0)
        
        ## Build C1 matrix, define outcome
        response=dfcopy_trtgrp[, y.varname]
        covariateX = cbind(model.matrix(outc_model, dfcopy_trtgrp), phat_trtgrp)
        
        ## Fit imputation model
        datastan <- list(
          N = nrow(covariateX),
          num_knot= num.knot,
          numpred = ncol(covariateX),
          y = response,
          pred_mat= covariateX,
          spline_mat= linearB
        )
        psppM <- mod$optimize(
          data = datastan,
          jacobian = F,
          iter= 10000
        )
        estcoef[[trtgrp+1]] <- list(fixed=psppM$summary("beta")$estimate,
                                    random = psppM$summary("gamma")$estimate)
        
        
        impute_potenout <-
          imputeF(
            newdata=samp,
            fixedcoef=psppM$summary("beta")$estimate,
            splinecoef = psppM$summary("gamma")$estimate,
            knotloc = knots,
            outcome.model =outc_model,
            propen.score.new=phat_samp
          )
        
        if(trtgrp == 0){
          imputed_ctrl <- impute_potenout
          knotlocations[[1]] <- knots
        }
        if(trtgrp == 1){
          imputed_trt<- impute_potenout
          knotlocations[[2]] <- knots
        }
        
      }
      coeflistres <- list(
        fixed0 = estcoef[[1]],
        fixed1 = estcoef[[2]],
        rand0 = estcoef[[1]]$random,
        rand1 = estcoef[[2]]$random,
        knotloc0 = knotlocations[[1]],
        knotloc1 = knotlocations[[2]]
      )
      imputed_ctrl[samp[,trt.varname]==0] <- samp[samp[,trt.varname]==0, y.varname]
      imputed_trt[samp[,trt.varname]==1] <- samp[samp[,trt.varname]==1, y.varname]
      
    # (4) return imputed values 
      return(list(imputed_ctrl, imputed_trt, coeflistres=coeflistres, phat_samp = phat_samp
      ))
  }, error=function(e) return(NA) )
}

