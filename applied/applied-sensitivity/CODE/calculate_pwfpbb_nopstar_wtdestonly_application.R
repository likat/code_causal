## FUNCTION FOR PENCOMP
calculate.pencomp.wfpbb <- function(propen_model, # string for propensity model, e.g. "A~X+Z"
                              outc_model,
                              samp = dat,
                              y.varname,    # name of the outcome variable
                              x.varnames,     # name of the confounding variables                              
                              num.knot=5,
                              F_draw=1,
                              trt.varname,
                              T_pop = Tfact){
  knotlocations <- vector(mode="list", length = 2)
    
      ## initialize empty list for imputation storage
      impute_potenout <- vector(mode="list", length=2)
      boot_all$impute <- 0
      imputed_ctrl <- rep(NA,nrow(samp))
      imputed_trt <- rep(NA,nrow(samp))
    
      ## WFPBB generalizing step
      bootwts <- T_pop*nrow(boot_all) * normalize(boot_all$wts*boot_all$clusmultiplier)
      coefflist0 <- vector(mode="list", length=F_draw)
      coefflist1 <- vector(mode="list", length=F_draw)

      for(f in 1:F_draw){
        tryCatch ( {        
          # temp <- wtpolyap(ysamp = 1:nrow(boot_all), wts = bootwts, k = npop - nrow(boot_all) )
          temp <- wtpolyap(ysamp = 1:nrow(boot_all), wts = bootwts, k = T_pop*nrow(boot_all) - nrow(boot_all) )
          popbootpolya <- boot_all[temp,]
          mod_trt <- glm(propen_model, data = popbootpolya, family="binomial")
          popbootpolya$phat <- predict(mod_trt, newdata=popbootpolya,type="link")
          phat_samp <- predict(mod_trt, newdata = samp, type = "link")

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
            response=dfcopy_trtgrp[, y.varname] %>% unlist %>% as.numeric
            covariateX = cbind(model.matrix(outc_model, dfcopy_trtgrp), linearB)
            
            ## Fit imputation model
            all=rep(1, nrow(dfcopy_trtgrp))
            ## SINGULARITY ISSUES HERE
            if(Matrix::rankMatrix(covariateX)[1] < ncol(covariateX)){bootflag <<- 1; break}
            psppM <- glm(response~covariateX-1, family = "binomial")
            if(trtgrp==0){
              coefflist0[[f]] <- psppM$coefficients
            }else{
              coefflist1[[f]] <- psppM$coefficients
            }
            
            ## save imputations for the other treatment group
            impute_potenout <- 
              imputeF(newdata=samp[samp[,trt.varname]!=trtgrp,],
                      model = psppM,
                      knotloc = knots,
                      outcome.model =outc_model,
                      propen.score.new=phat_samp[samp[,trt.varname]!=trtgrp])
            if(trtgrp == 0){
              imputed_ctrl[samp[,trt.varname] == 0] <- samp[samp[,trt.varname] == 0, y.varname] %>% 
                unlist %>% as.numeric
              imputed_ctrl[samp[,trt.varname] == 1] <- impute_potenout
              knotlocations[[1]] <- knots
            }
            if(trtgrp == 1){
              imputed_trt[samp[,trt.varname] == 1] <- samp[samp[,trt.varname] == 1, y.varname]%>% 
                unlist %>% as.numeric
              imputed_trt[samp[,trt.varname] == 0] <- impute_potenout
              knotlocations[[2]] <- knots
            }
          }
        }, error=function(e) return(NA) )
        
      }
      coeff1 <- rep(0,ncol(covariateX))
      coeff0 <- rep(0, ncol(covariateX))
      randeff1 <- rep(0, num.knot)
      randeff0 <- rep(0, num.knot)
      indskip <- 0
      for(ind in 1:F_draw){
        if(length(coefflist1[[ind]])==0 |length(coefflist0[[ind]])==0 ){
          indskip = indskip + 1; next
        }
        # coeff0 <- coeff0 + coefflist0[[ind]]$fixed
        # coeff1 <- coeff1 + coefflist1[[ind]]$fixed
        # randeff0 <- randeff0 + coefflist0[[ind]]$random$all
        # randeff1 <- randeff1 + coefflist1[[ind]]$random$all
        coeff0 <- coeff0 + coefflist0[[ind]]
        coeff1 <- coeff1 + coefflist1[[ind]]
      }
      coeff0 <- coeff0/(F_draw-indskip)
      coeff1 <- coeff1/(F_draw-indskip)
      # randeff0 <- randeff0/(F_draw-indskip)
      # randeff1 <- randeff1/(F_draw-indskip)

      coeflistres <- list(
        fixed0 = coeff0,
        fixed1 = coeff1,
        # rand0 = randeff0,
        # rand1 = randeff1,
        knotloc0 = knotlocations[[1]],
        knotloc1 = knotlocations[[2]]
      )
      
      # (4) return imputed values 
      return(list(imputed0 = imputed_ctrl, imputed1= imputed_trt,
                  phat_samp = phat_samp,
                  phat_fitdat = popbootpolya$phat,
                  coeflistres=coeflistres
      ))
}
  
  
  



