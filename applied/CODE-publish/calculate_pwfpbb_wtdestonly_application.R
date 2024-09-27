## FUNCTION FOR PENCOMP
calculate.pencomp.wfpbb <- function(propen_model, # string for propensity model, e.g. "A~X+Z"
                              outc_model,
                              bootdf,
                              samp = dat,
                              y.varname,    # name of the outcome variable
                              x.varnames,     # name of the confounding variables                              
                              num.knot=5,
                              F_draw=1,
                              trt.varname,
                              T_pop = Tfact){
  
    
      ## WFPBB generalizing step
      bootwts <- T_pop*nrow(bootdf) * normalize(bootdf$wts*bootdf$bbclusmultiplier)
      coefflist0 <- vector(mode="list", length=F_draw)
      coefflist1 <- vector(mode="list", length=F_draw)

      ## initialize empty list for imputation storage
      impute_potenout <- vector(mode="list", length=2)
      bootdf$impute <- 0
      ate_pop <- list(estimate = rep(NA,F_draw))
      ate_samp <-  list(estimate = rep(NA,F_draw),
                        var = rep(NA, F_draw))
      wtd_est <- list(estimate = rep(NA,F_draw),
                      var = rep(NA, F_draw))
      knotlocations <- vector(mode="list", length = 2)
      
      for(f in 1:F_draw){
        imputed_ctrl <- rep(NA,nrow(samp))
        imputed_trt <- rep(NA,nrow(samp))
        imputed_ctrlpop <- rep(NA,T_pop*nrow(bootdf))
        imputed_trtpop <- rep(NA, T_pop*nrow(bootdf))
        
        tryCatch ( {        
          
          temp <- 
            wtpolyap(ysamp = 1:nrow(bootdf), wts = bootwts, k = T_pop*nrow(bootdf) - nrow(bootdf) )
          popbootpolya <- 
            bootdf[temp,]
          mod_trt <- 
            glm(propen_model, data = popbootpolya, family="binomial")
          popbootpolya$phat <-
            predict(mod_trt, newdata=popbootpolya,type="link")
          phat_samp <- 
            predict(mod_trt, newdata = samp, type = "link")

          for(trtgrp in 0:1){
            ## separate fitted probabilities by treatment group
            phat_trtgrp <- popbootpolya$phat[popbootpolya[,trt.varname]==trtgrp]
            dfcopy_trtgrp <- popbootpolya[popbootpolya[,trt.varname]==trtgrp,]
            
            # (3) Fit outcome model with splines
            ## currently borrowing code from zhou 2019
            space=(max(phat_trtgrp)-min(phat_trtgrp))/(num.knot+1)
            knots=(min(phat_trtgrp)+space*(1:num.knot))
            
            ## Construct C2 matrix, truncated linear basis splines
            linearB=NULL
            linearB =outer(phat_trtgrp, knots, "-") ## constructing C2 matrix
            linearB =linearB * (linearB > 0)
            
            ## Build C1 matrix, define outcome
            # outc_model <- as.formula("adltfs ~ pctpov +age+ educ + race")
            response=dfcopy_trtgrp[, y.varname] %>% unlist %>% as.numeric
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
              iter = 10000
            )
            
            if(trtgrp==0){
              coefflist0[[f]] <- c(psppM$summary("beta")$estimate, 
                                   psppM$summary("gamma")$estimate)
              
            }else{
              coefflist1[[f]] <- c(psppM$summary("beta")$estimate, 
                                   psppM$summary("gamma")$estimate)
            }
            
            ## mean function for all values (observed and nonobserved)
            impute_potenout <- 
              imputeF(
                newdata=samp,
                  # model = psppM,
                fixedcoef=psppM$summary("beta")$estimate,
                splinecoef = psppM$summary("gamma")$estimate,
                knotloc = knots,
                outcome.model =outc_model,
                propen.score.new=phat_samp
              )
            
            ## imputations for synthetic population
            impute_potenoutpop <- 
              imputeF(
                newdata= popbootpolya,
                # model = psppM,
                fixedcoef=psppM$summary("beta")$estimate,
                splinecoef = psppM$summary("gamma")$estimate,
                knotloc = knots,
                outcome.model =outc_model,
                propen.score.new=popbootpolya$phat
              )
            
            if(trtgrp == 0){
              # imputed_ctrl[samp[,trt.varname] == 0] <- 
              #   samp[samp[,trt.varname] == 0, y.varname] %>% 
              #   unlist %>% as.numeric
              # imputed_ctrl[samp[,trt.varname] == 1] <- impute_potenout
              # imputed_ctrlpop[popbootpolya[,trt.varname] == 0] <- 
              #   popbootpolya[popbootpolya[,trt.varname] == 0, y.varname] %>% 
              #   unlist %>% as.numeric
              # imputed_ctrlpop[popbootpolya[,trt.varname] == 1] <- impute_potenoutpop
              imputed_ctrl <- impute_potenout
              imputed_ctrlpop <- impute_potenoutpop
              knotlocations[[1]] <- knots
            }
            if(trtgrp == 1){
              # imputed_trt[samp[,trt.varname] == 1] <- samp[samp[,trt.varname] == 1, y.varname]%>% 
              #   unlist %>% as.numeric
              # imputed_trt[samp[,trt.varname] == 0] <- impute_potenout
              # 
              # imputed_trtpop[popbootpolya[,trt.varname] == 1] <- 
              #   popbootpolya[popbootpolya[,trt.varname] == 1, y.varname] %>% 
              #   unlist %>% as.numeric
              # imputed_trtpop[popbootpolya[,trt.varname] == 0] <- impute_potenoutpop
              imputed_trt <- impute_potenout
              imputed_trtpop <- impute_potenoutpop
              knotlocations[[2]] <- knots
            }
          }
          imputed_ctrl[samp[,trt.varname]==0] <- samp[samp[,trt.varname]==0, y.varname]
          imputed_trt[samp[,trt.varname]==1] <- samp[samp[,trt.varname]==1, y.varname]
          
          imputed_ctrlpop[popbootpolya[,trt.varname]==0] <- popbootpolya[popbootpolya[,trt.varname]==0, y.varname]
          imputed_trtpop[popbootpolya[,trt.varname]==1] <- popbootpolya[popbootpolya[,trt.varname]==1, y.varname]
          
          ate_samp$estimate[f] <- mean(imputed_trt - imputed_ctrl)
          ate_samp$var[f] <- var(imputed_trt - imputed_ctrl)/nrow(samp)
          ate_pop$estimate[f] <- mean(imputed_trtpop-imputed_ctrlpop)
          wtdtemp <- calculate.wtd(imputedlist = list(imputed_ctrl, imputed_trt),
                                   bootdf = samp)
          wtd_est$estimate[f] <- wtdtemp[[1]]
          wtd_est$var[f] <- wtdtemp[[2]]
        }, error=function(e) return(NA) )
        
      }
      coeff1 <- rep(0,ncol(covariateX)+num.knot)
      coeff0 <- rep(0, ncol(covariateX) + num.knot)
      randeff1 <- rep(0, num.knot)
      randeff0 <- rep(0, num.knot)
      indskip <- 0
      for(ind in 1:F_draw){
        if(length(coefflist1[[ind]])==0 |length(coefflist0[[ind]])==0 ){
          indskip = indskip + 1; next
        }
        coeff0 <- coeff0 + coefflist0[[ind]]
        coeff1 <- coeff1 + coefflist1[[ind]]
      }
      coeff0 <- coeff0/(F_draw-indskip)
      coeff1 <- coeff1/(F_draw-indskip)

      coeflistres <- list(
        coef0 = coeff0,
        coef1 = coeff1,
        knotloc0 = knotlocations[[1]],
        knotloc1 = knotlocations[[2]]
      )
      
      # (4) return imputed values 
      return(list(imputed0 = imputed_ctrl, imputed1 = imputed_trt, 
                  ate_samp=lapply(ate_samp,mean), ate_pop=lapply(ate_pop,mean),
                  wtdest = lapply(wtd_est, mean),
                  phat_samp = phat_samp,
                  phat_fitdat = popbootpolya$phat,
                  coeflistres=coeflistres
      ))
}
  
  
  



