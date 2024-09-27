## FUNCTION FOR PENCOMP
calculate.pencomp <- function(propen_model, # string for propensity model, e.g. "A~X+Z"
                              outc_model,
                              popdat=pop,
                              trt.varname = "A",
                              bootind0 = bootind_cnt,
                              bootind1 = bootind_trt,
                              sampind0 = ind_cnt,
                              sampind1 = ind_trt,
                              y.varname,    # name of the outcome variable
                              # x.varnames,     # name of the confounding variables                              
                              num.knot=10){
  # tryCatch ( {
  knotlocations <- vector(mode="list", length = 2)
  # (2) Fit propensity score model 
    ## What if we used the true values for propensity model and outcome model?
    ### it's NOT an issue with the propensity scores;
    ### using the true ptrt biases PATE equally
    mod_trt <- glm(propen_model, data = boot_all, family="binomial")
    phat <- predict(mod_trt, newdata=boot_all,type="link")
    phat_pop <- predict(mod_trt, newdata= popdat, type = "link")
    phat_samp <- predict(mod_trt, newdata = samp, type = "link")
    boot_all$phat <- phat
    
    # hist(phat_pop - popdat$ptrt)
    
    ## initialize empty list for imputation storage
    impute_potenout <- vector(mode="list", length=2)
    boot_all$impute <- 0
    imputed_ctrl <- rep(NA,nrow(samp))
    imputed_trt <- rep(NA,nrow(samp))
    imputed_ctrl_pop <- rep(NA,nrow(popdat))
    imputed_trt_pop <- rep(NA,nrow(popdat))
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
        popmod <- formulaF(varList = varnames_outcome,y.name = "Y0")
        impute_potenoutpop <-
          imputeF(newdata=popdat,
                  fixedcoef=psppM$summary("beta")$estimate,
                  splinecoef = psppM$summary("gamma")$estimate,
                  knotloc = knots,
                  outcome.model =popmod,
                  propen.score.new=phat_pop)
        imputed_ctrl_pop <- impute_potenoutpop
        knotlocations[[1]] <- knots
      }
      if(trtgrp == 1){
        imputed_trt<- impute_potenout
        popmod <- formulaF(varList = varnames_outcome,y.name = "Y1")
        impute_potenoutpop <-
          imputeF(newdata=popdat,
                  fixedcoef=psppM$summary("beta")$estimate,
                  splinecoef = psppM$summary("gamma")$estimate,
                  knotloc = knots,
                  outcome.model =popmod,
                  propen.score.new=phat_pop)
        imputed_trt_pop <- impute_potenoutpop
        knotlocations[[2]] <- knots
      }
      
    }
    coeflistres <- list(
      fixed0 = estcoef[[1]]$fixed,
      fixed1 = estcoef[[2]]$fixed,
      rand0 = estcoef[[1]]$random,
      rand1 = estcoef[[2]]$random,
      knotloc0 = knotlocations[[1]],
      knotloc1 = knotlocations[[2]]
    )
    imputed_ctrl_pop[pop$I==1][samp[,trt.varname]==0] <- samp$Ysamp[samp[,trt.varname]==0]
    imputed_trt_pop[pop$I==1][samp[,trt.varname]==1] <- samp$Ysamp[samp[,trt.varname]==1]
    imputed_ctrl[samp[,trt.varname]==0] <- samp$Ysamp[samp[,trt.varname]==0]
    imputed_trt[samp[,trt.varname]==1] <- samp$Ysamp[samp[,trt.varname]==1]
    
    pate_gen <- 
      mean(imputed_trt_pop - imputed_ctrl_pop)
    varpategen <- var(imputed_trt_pop - imputed_ctrl_pop)/npop
    # (4) return imputed values 
    return(list(imputed0 = imputed_ctrl, imputed1 = imputed_trt, modtrt = mod_trt, 
                pate_gen=pate_gen,
                varpategen = varpategen,
                coeflistres = coeflistres,
                knotlocations= knotlocations,
                phat = phat,
                phat_samp = phat_samp,
                phat_pop = phat_pop
    ))
  # }, error=function(e) return(NA) )
  
}

