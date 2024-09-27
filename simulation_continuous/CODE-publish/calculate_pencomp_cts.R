 ## FUNCTION TO IMPUTE 
# imputeF=function(newdata, 
#                  model, 
#                  # x.varnames, 
#                  knotloc = knots,
#                  outcome.model,
#                  propen.score.new) {
#   
#   # knots=model$knot.loc
#   # space=(max(propen.score.new)-min(propen.score.new))/(num.knot+1)
#   # knots=(min(propen.score.new)+space*(1:num.knot))
#   
#   linearB=NULL
#   linearB =outer(propen.score.new, knotloc, "-")
#   linearB =linearB * (linearB > 0)
#   
#   # designM = cbind(model.matrix(outcome.model, newdata), propen.score.new)
#   designM = cbind(model.matrix(outcome.model, newdata))
#   designM=as.matrix(designM)
#   fittedcoef <- model$coefficients$fixed
#   randomcoef <- model$coefficients$random$all 
#   # predicted = designM %*% model$coefficients$fixed + as.matrix(linearB) %*% as.vector(unlist(model$random)) + rnorm(nrow(newdata), 0, model$sigma)
#   predicted = designM %*% fittedcoef + as.matrix(linearB) %*% as.numeric(randomcoef) + rnorm(nrow(newdata), 0, model$sigma)
#   
#   return(predicted)
#   
# }


  
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
  knotlocations <- vector(mode="list", length = 2)
  tryCatch ( {
      # (2) Fit propensity score model 
    ## What if we used the true values for propensity model and outcome model?
      ### it's NOT an issue with the propensity scores;
     ### using the true ptrt biases PATE equally
      mod_trt <- glm(propen_model, data = boot_all, family="binomial")
      phat <- predict(mod_trt, newdata=boot_all,type="link")
      phat_pop <- predict(mod_trt, newdata= popdat, type = "link")
      phat_samp <- predict(mod_trt, newdata = samp, type = "link")
      
      # phat <- boot_all$ptrt
      # phat_pop <- popdat$ptrt
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
        # covariateX = cbind(model.matrix(outc_model, dfcopy_trtgrp))
        
        ## Fit imputation model
        all=rep(1, nrow(dfcopy_trtgrp))
        ## SINGULARITY ISSUES HERE
        if(Matrix::rankMatrix(covariateX)[1] < ncol(covariateX)){bootflag <<- 1; break}
        # if(sum(table(dfcopy_trtgrp$catecat) < 10) > 0){bootflag <<- 1; break}
        psppM=lme(response ~ covariateX-1, random=list(all=pdIdent(~0+linearB)),
                  control=(maxIter=100))
        # lm(response~covariateX-1)
        estcoef[[trtgrp+1]] <- psppM$coefficients
        
        ## Fit imputation model
        # datastan <- list(
        #   N = nrow(covariateX),
        #   num_knot= num.knot,
        #   numpred = ncol(covariateX),
        #   y = response,
        #   pred_mat= covariateX,
        #   spline_mat= linearB
        # )
        # psppM <- mod$optimize(
        #   data = datastan,
        #   jacobian = F,
        #   iter= 10000
        # )
        # estcoef[[trtgrp+1]] <- list(fixed=psppM$summary("beta")$estimate,
        #                             random = psppM$summary("gamma")$estimate,
        #                             sigma_y = psppM$summary("sigma_y")$estimate)
        
        impute_potenout <-
          imputeF(newdata=samp[samp[,trt.varname]!=trtgrp,],
                  model = psppM,
                  knotloc = knots,
                  outcome.model =outc_model,
                  propen.score.new=phat_samp[samp[,trt.varname]!=trtgrp])
        ## save imputations for the other treatment group
        # impute_potenout <-
        #   imputeF(
        #     newdata=samp[samp[,trt.varname]!=trtgrp,],
        #     fixedcoef=psppM$summary("beta")$estimate,
        #     splinecoef = psppM$summary("gamma")$estimate,
        #       sigmay = psppM$summary("sigma_y")$estimate,
        #     knotloc = knots,
        #     outcome.model =outc_model,
        #     propen.score.new=phat_samp[samp[,trt.varname]!=trtgrp]
        #   )
        
        
        if(trtgrp == 0){
          imputed_ctrl[samp[,trt.varname] == 0] <- samp$Ysamp[samp[,trt.varname] == 0]
          imputed_ctrl[samp[,trt.varname] == 1] <- impute_potenout
          # formulaF=function(varList, y.name){
          #   return ( as.formula(paste(y.name, "~ ", paste(c(varList), collapse = "+"))) )
          # }
          popmod <- formulaF(varList = varnames_outcome,y.name = "Y0")
          impute_potenoutpop <-
            imputeF(newdata=popdat,
                    model = psppM,
                    knotloc = knots,
                    outcome.model =popmod,
                    propen.score.new=phat_pop)
         # impute_potenoutpop <- imputeF(
         #    newdata=popdat,
         #    fixedcoef=psppM$summary("beta")$estimate,
         #    splinecoef = psppM$summary("gamma")$estimate,
         #    sigmay = psppM$summary("sigma_y")$estimate,
         #    knotloc = knots,
         #    outcome.model =popmod,
         #    propen.score.new=phat_pop
         #  )
          
          imputed_ctrl_pop <- impute_potenoutpop
          knotlocations[[1]] <- knots
        }
        if(trtgrp == 1){
          imputed_trt[samp[,trt.varname] == 1] <- samp$Ysamp[samp[,trt.varname] == 1]
          imputed_trt[samp[,trt.varname] == 0] <- impute_potenout
          popmod <- formulaF(varList = varnames_outcome,y.name = "Y1")
          impute_potenoutpop <-
            imputeF(newdata=popdat,
                    model = psppM,
                    knotloc = knots,
                    outcome.model =popmod,
                    propen.score.new=phat_pop)

          # impute_potenoutpop <- imputeF(
          #   newdata=popdat,
          #   fixedcoef=psppM$summary("beta")$estimate,
          #   splinecoef = psppM$summary("gamma")$estimate,
          #   sigmay = psppM$summary("sigma_y")$estimate,
          #   knotloc = knots,
          #   outcome.model =popmod,
          #   propen.score.new=phat_pop
          # )
          imputed_trt_pop <- impute_potenoutpop
          knotlocations[[2]] <- knots
          
        }
        
      }
      # hist((imputed_ctrl_pop - pop$Y1)/mean(pop$Y1))
      imputed_ctrl_pop[pop$I==1][samp[,trt.varname]==0] <- samp$Ysamp[samp[,trt.varname]==0]
      imputed_trt_pop[pop$I==1][samp[,trt.varname]==1] <- samp$Ysamp[samp[,trt.varname]==1]
      
      poptemp <- data.frame(
        ind = pop$catecat,
        imputed0 = imputed_ctrl_pop,
        imputed1 = imputed_trt_pop
      )
      ## boot level ATE
      temptblboot <- data.frame(ind = boot_all$catecat,
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
      coeflistres <- list(
        fixed0 = estcoef[[1]]$fixed,
        fixed1 = estcoef[[2]]$fixed,
        rand0 = estcoef[[1]]$random$all,
        rand1 = estcoef[[2]]$random$all
      )
      # coeflistres <- list(
      #   fixed0 = estcoef[[1]]$fixed,
      #   fixed1 = estcoef[[2]]$fixed,
      #   rand0 = estcoef[[1]]$random,
      #   rand1 = estcoef[[2]]$random,
      #   sigma_y0 = estcoef[[1]]$sigma_y,
      #   sigma_y1 = estcoef[[2]]$sigma_y
      # )
      

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
  }, error=function(e) return(NA) )
}

# ## PRINTING OUT A TABLE OF ESTIMATED COEFFICIENTS FOR each group
# mod_trt <- glm(propensity_mod, data = boot_all, family="binomial")
# phat <- predict(mod_trt, newdata=boot_all,type="link")
# phat_pop <- predict(mod_trt, newdata= popdat, type = "link")
# # phat <- boot_all$ptrt
# # phat_pop <- popdat$ptrt
# boot_all$phat <- phat
# 
# ## SAMPLE COEF
# trtgrp <- 0
# outc_model <- "Ysamp ~ X2 + Z" %>% as.formula()
#     phat_trtgrp <- boot_all$phat[boot_all$A==trtgrp]
#     dfcopy_trtgrp <- boot_all[boot_all$A==trtgrp,]
# 
#     # (3) Fit outcome model with splines
#     ## need sensitivity analysis for how many knots
#     ## currently borrowing code from zhou 2019
#     space=(max(phat_trtgrp)-min(phat_trtgrp))/(num.knot+1)
#     knots=(min(phat_trtgrp)+space*(1:num.knot))
# 
#     ## Construct C2 matrix, truncated linear basis splines
#     linearB=NULL
#     linearB =outer(phat_trtgrp, knots, "-") ## constructing C2 matrix
#     linearB =linearB * (linearB > 0)
# 
#     ## Build C1 matrix, define outcome
#     response=dfcopy_trtgrp[, "Ysamp"]
#     covariateX = cbind(model.matrix(outc_model, dfcopy_trtgrp), phat_trtgrp)
# 
#     ## Fit imputation model
#     all=rep(1, nrow(dfcopy_trtgrp))
#     ## SINGULARITY ISSUES HERE
#     if(Matrix::rankMatrix(covariateX)[1] < ncol(covariateX)){bootflag <<- 1; break}
#     # if(sum(table(dfcopy_trtgrp$catecat) < 10) > 0){bootflag <<- 1; break}
#     psppM=lme(response ~ covariateX-1, random=list(all=pdIdent(~0+linearB)),
#               control=(maxIter=100))
#     print(xtable(summary(psppM)$tTable, digits=3, type = "latex"), file = "coeff_pencomp_nogb_srs_bootY0.tex")
# 
# trtgrp <- 1
#     phat_trtgrp <- boot_all$phat[boot_all$A==trtgrp]
#     dfcopy_trtgrp <- boot_all[boot_all$A==trtgrp,]
# 
#     # (3) Fit outcome model with splines
#     ## need sensitivity analysis for how many knots
#     ## currently borrowing code from zhou 2019
#     space=(max(phat_trtgrp)-min(phat_trtgrp))/(num.knot+1)
#     knots=(min(phat_trtgrp)+space*(1:num.knot))
# 
#     ## Construct C2 matrix, truncated linear basis splines
#     linearB=NULL
#     linearB =outer(phat_trtgrp, knots, "-") ## constructing C2 matrix
#     linearB =linearB * (linearB > 0)
# 
#     ## Build C1 matrix, define outcome
#     response=dfcopy_trtgrp[, "Ysamp"]
#     covariateX = cbind(model.matrix(outc_model, dfcopy_trtgrp), phat_trtgrp)
# 
#     ## Fit imputation model
#     all=rep(1, nrow(dfcopy_trtgrp))
#     ## SINGULARITY ISSUES HERE
#     if(Matrix::rankMatrix(covariateX)[1] < ncol(covariateX)){bootflag <<- 1; break}
#     # if(sum(table(dfcopy_trtgrp$catecat) < 10) > 0){bootflag <<- 1; break}
#     psppM=lme(response ~ covariateX-1, random=list(all=pdIdent(~0+linearB)),
#               control=(maxIter=100))
#     print(xtable(summary(psppM)$tTable, digits=3, type = "latex"), file = "coeff_pencomp_nogb_srs_bootY1.tex")
# 
# ## POP COEF
#     logit <- function(x){return(log(x/(1-x)))}
# phat_trtgrp <- popdat$ptrt %>% logit()
# dfcopy_trtgrp <- popdat
# outc_model = "Y1 ~ X2 + Z" %>% as.formula()
# y.varname = "Y1"
#     # # (3) Fit outcome model with splines
#     # ## need sensitivity analysis for how many knots
#     # ## currently borrowing code from zhou 2019
#     space=(max(phat_trtgrp)-min(phat_trtgrp))/(num.knot+1)
#     knots=(min(phat_trtgrp)+space*(1:num.knot))
# 
#     # ## Construct C2 matrix, truncated linear basis splines
#     linearB=NULL
#     linearB =outer(phat_trtgrp, knots, "-") ## constructing C2 matrix
#     linearB =linearB * (linearB > 0)
# 
#     ## Build C1 matrix, define outcome
#     response=dfcopy_trtgrp[, y.varname]
#     covariateX = cbind(model.matrix(outc_model, dfcopy_trtgrp), phat_trtgrp)
# 
#     ## Fit imputation model
#     all=rep(1, nrow(dfcopy_trtgrp))
#     ## SINGULARITY ISSUES HERE
#     psppM=lme(response ~ covariateX-1, random=list(all=pdIdent(~0+linearB)),
#               control=(maxIter=100))
#     print(xtable(summary(psppM)$tTable, digits=3, type = "latex"), file = "coeff_pencomp_nogb_overlap_popY1.tex")
# 
# outc_model = "Y0 ~ X2 + Z" %>% as.formula()
# y.varname = "Y0"
#     # (3) Fit outcome model with splines
#     ## need sensitivity analysis for how many knots
#     ## currently borrowing code from zhou 2019
#     space=(max(phat_trtgrp)-min(phat_trtgrp))/(num.knot+1)
#     knots=(min(phat_trtgrp)+space*(1:num.knot))
# 
#     ## Construct C2 matrix, truncated linear basis splines
#     linearB=NULL
#     linearB =outer(phat_trtgrp, knots, "-") ## constructing C2 matrix
#     linearB =linearB * (linearB > 0)
# 
#     ## Build C1 matrix, define outcome
#     response=dfcopy_trtgrp[, y.varname]
#     covariateX = cbind(model.matrix(outc_model, dfcopy_trtgrp), phat_trtgrp)
# 
#     ## Fit imputation model
#     all=rep(1, nrow(dfcopy_trtgrp))
#     ## SINGULARITY ISSUES HERE
#     psppM=lme(response ~ covariateX-1, random=list(all=pdIdent(~0+linearB)),
#               control=(maxIter=100))
#     print(xtable(summary(psppM)$tTable, digits=3, type = "latex"), file = "coeff_pencomp_nogb_overlap_popY0.tex")
# 

