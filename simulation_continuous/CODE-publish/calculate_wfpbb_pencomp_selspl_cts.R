## FUNCTION FOR PENCOMP
calculate.pwfpbb.pselspl <- function(propen_model, # string for propensity model, e.g. "A~X+Z"
                                    outc_model,
                                    popdat=pop,
                                    bootdf = boot_allbb,
                                    uniqueind,
                                    polyawts,
                                    # covariateXnames,
                                    # propensity_mod_aiptw,
                                    # outcome_mod_aiptw,
                                    bootind0 = bootind_cnt,
                                    bootind1 = bootind_trt,
                                    sampind0 = ind_cnt,
                                    sampind1 = ind_trt,
                                    trt.varname="A",
                                    y.varname,    # name of the outcome variable
                                    x.varnames,     # name of the confounding variables                              
                                    num.knot=10,
                                    F_draw=1
                                    # nCarlo = 2000
                                    ){
  knotlocations <- vector(mode="list", length = 2)
  # (2) Fit propensity score model 
  ## What if we used the true values for propensity model and outcome model?
  ### it's NOT an issue with the propensity scores;
  ### using the true ptrt biases PATE equally
  # hist(phat_pop - popdat$ptrt)
  
  ## initialize empty list for imputation storage
  impute_potenout <- vector(mode="list", length=2)
  bootdf <- boot_allbb
  bootdf$id <- 1:nrow(bootdf)
  # table(bootdf$id)
  bootdf$impute <- 0
  
  ## WFPBB generalizing step
  # bootwts <- npop * normalize(bootdf$wts)
  bootwts <- npop * normalize(bootdf$wts)
  pategen_f <- rep(NA,F_draw)
  pategen_synthpop_f <- rep(NA, F_draw)
  wtdest_f <- matrix(ncol=2, nrow= F_draw)
  varpategen_f <- rep(NA,F_draw)
  varpategen_synthpop_f <- rep(NA, F_draw)
  coefflist0 <- vector(mode="list", length=F_draw)
  coefflist1 <- vector(mode="list", length=F_draw)
  knotloc0 <- vector(mode="list", length = F_draw)
  knotloc1 <- vector(mode="list", length = F_draw)
  pategen_direct <- rep(NA, F_draw)
  time_pencompwfpbb <- rep(NA,F_draw)
  # covariateX <- NA
  for(f in 1:F_draw){
    tryCatch ( {        
      # temp <- wtpolyap(ysamp = 1:nrow(bootdf), wts = bootwts, k = npop - nrow(bootdf) )
      temp <- wtpolyap(ysamp = uniqueind, wts = polyawts, k = npop - length(uniqueind) )
      
      # popbootpolya <- bootdf[temp,]
      popbootpolya <- samp[temp,]
      pategen_direct[f] <- mean(popbootpolya$Y1 - popbootpolya$Y0)
      # bt_polya <- popbootpolya[sample(1:npop, size = nrow(samp), replace=F),]
      
      ## PENCOMP
      t0_pencomp <- Sys.time()
      phat_imputed <- calculate.ptrtspline(num.knot = 3,
                                           propen.model = propen_model,
                                           popdat2 = popdat,
                                           fitdat = popbootpolya,
                                           trtname = trt.varname,
                                           sampdat = samp)
      
      popbootpolya$phat <- phat_imputed$ptrt_synthpop
      phat_pop <- as.numeric(phat_imputed$ptrt_pop)
      phat_samp <- as.numeric(phat_imputed$ptrt_samp)
      imputed_ctrl_pop <- rep(NA,nrow(popdat))
      imputed_trt_pop <- rep(NA,nrow(popdat))
      imputed_ctrl_synthpop <- rep(NA,nrow(popdat))
      imputed_trt_synthpop <- rep(NA,nrow(popdat))
      imputed_ctrl <- rep(NA,nrow(samp))
      imputed_trt <- rep(NA,nrow(samp))
      
      for(trtgrp in 0:1){
        ## separate fitted probabilities by treatment group
        phat_trtgrp <- popbootpolya$phat[popbootpolya[,trt.varname]==trtgrp]
        dfcopy_trtgrp <- popbootpolya[popbootpolya[,trt.varname]==trtgrp,]
        
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
        all=rep(1, nrow(dfcopy_trtgrp))
        ## SINGULARITY ISSUES HERE
        if(!exists("covariateX")){bootflag <<- 1; break}
        # if(sum(table(dfcopy_trtgrp$catecat) < 10) > 0){bootflag <<- 1; break}
        psppM=lme(response ~ covariateX-1, random=list(all=pdIdent(~0+linearB)),
                  control=(maxIter=100))
        
        
        ## COLLECT COEFFICIENTS FROM WFPBB
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

        ## imputations for synthetic population
        impute_potenoutsynthpop <- 
          imputeF(
            newdata= popbootpolya[popbootpolya[,trt.varname]!=trtgrp,],
            model = psppM,
            # fixedcoef=psppM$summary("beta")$estimate,
            # splinecoef = psppM$summary("gamma")$estimate,
            # sigmay = psppM$summary("sigma_y")$estimate,
            knotloc = knots,
            outcome.model =outc_model,
            propen.score.new=popbootpolya$phat[popbootpolya[,trt.varname]!=trtgrp]
          )
        
        if(trtgrp == 0){
          imputed_ctrl[samp[,trt.varname] == 0] <- samp$Ysamp[samp[,trt.varname] == 0]
          imputed_ctrl[samp[,trt.varname] == 1] <- impute_potenout
          # formulaF=function(varList, y.name){
          #   return ( as.formula(paste(y.name, "~ ", paste(c(varList), collapse = "+"))) )
          # }
          popmod <- formulaF(varList = varnames_outcome,y.name = "Y0")
          # impute_potenoutpop <- imputeF(
          #   newdata=popdat,
          #   fixedcoef=psppM$summary("beta")$estimate,
          #   splinecoef = psppM$summary("gamma")$estimate,
          #   sigmay = psppM$summary("sigma_y")$estimate,
          #   knotloc = knots,
          #   outcome.model =popmod,
          #   propen.score.new=phat_pop
          # )
          impute_potenoutpop <- 
            imputeF(
              newdata= popdat,
              model = psppM,
              knotloc = knots,
              outcome.model = popmod,
              propen.score.new=as.numeric(phat_pop)
            )
          
          imputed_ctrl_pop <- impute_potenoutpop
          imputed_ctrl_synthpop[popbootpolya[,trt.varname] == 0] <- 
            popbootpolya[popbootpolya[,trt.varname] == 0, y.varname] %>% 
            unlist %>% as.numeric
          imputed_ctrl_synthpop[popbootpolya[,trt.varname] == 1] <- impute_potenoutsynthpop
          
          knotlocations[[1]] <- knots
        }
        if(trtgrp == 1){
          imputed_trt[samp[,trt.varname] == 1] <- samp$Ysamp[samp[,trt.varname] == 1]
          imputed_trt[samp[,trt.varname] == 0] <- impute_potenout
          popmod <- formulaF(varList = varnames_outcome,y.name = "Y1")
          # impute_potenoutpop <- imputeF(
          #   newdata=popdat,
          #   fixedcoef=psppM$summary("beta")$estimate,
          #   splinecoef = psppM$summary("gamma")$estimate,
          #   sigmay = psppM$summary("sigma_y")$estimate,
          #   knotloc = knots,
          #   outcome.model =popmod,
          #   propen.score.new=phat_pop
          # )
          impute_potenoutpop <- 
            imputeF(
              newdata= popdat,
              model = psppM,
              knotloc = knots,
              outcome.model = popmod,
              propen.score.new=as.numeric(phat_pop)
            )
          
          imputed_trt_pop <- impute_potenoutpop
          imputed_trt_synthpop[popbootpolya[,trt.varname] == 1] <- 
            popbootpolya[popbootpolya[,trt.varname] == 1, y.varname] %>% 
            unlist %>% as.numeric
          imputed_trt_synthpop[popbootpolya[,trt.varname] == 0] <- impute_potenoutsynthpop
          
          knotlocations[[2]] <- knots
        }
      }
      # hist((imputed_ctrl_pop - pop$Y1)/mean(pop$Y1))
      imputed_ctrl_pop[pop$I==1][samp[,trt.varname]==0] <- samp$Ysamp[samp[,trt.varname]==0]
      imputed_trt_pop[pop$I==1][samp[,trt.varname]==1] <- samp$Ysamp[samp[,trt.varname]==1]
      t1_pencomp <- Sys.time()
      time_pencompwfpbb[f] <- t1_pencomp-t0_pencomp
      ## generalized PATE
      pategen_f[f] <- 
        mean(imputed_trt_pop- imputed_ctrl_pop)
      varpategen_f[f] <- var(imputed_trt_pop- imputed_ctrl_pop)/npop
      pategen_synthpop_f[f] <-
        mean(imputed_trt_synthpop - imputed_ctrl_synthpop)
      varpategen_synthpop_f[f] <- 
        var(imputed_trt_synthpop - imputed_ctrl_synthpop)/npop
      
      ## weighted est
      samp$ydiff <- imputed_trt-imputed_ctrl
      des <- svydesign(id=~0, weights=~wts,data=samp)
      patesvymean <- svymean(~ydiff, design=des) 
      wtdest_f[f,1] <-  patesvymean %>% as.numeric()
      wtdest_f[f,2] <- SE(patesvymean)^2
      
      
      # print(f)
    }, error=function(e) return(NA) )
    
  }
  if(is.na(ncol(covariateX))){res <- NA}else{
    

  coeff1 <- rep(0,ncol(covariateX))
  coeff0 <- rep(0, ncol(covariateX))
  # sigmay1 <- 0
  # sigmay0 <- 0
  
  randeff1 <- rep(0, num.knot)
  randeff0 <- rep(0, num.knot)
  
  indskip <- 0
  for(ind in 1:F_draw){
    if(length(coefflist1[[ind]])==0 |length(coefflist0[[ind]])==0 ){
      indskip = indskip + 1; next
    }
    coeff0 <- coeff0 + coefflist0[[ind]]$fixed
    coeff1 <- coeff1 + coefflist1[[ind]]$fixed
    randeff0 <- randeff0 + coefflist0[[ind]]$random$all
    randeff1 <- randeff1 + coefflist1[[ind]]$random$all
    
    # randeff0 <- randeff0 + coefflist0[[ind]]$random
    # randeff1 <- randeff1 + coefflist1[[ind]]$random
    # sigmay0 <- sigmay0 + coefflist0[[ind]]$sigma_y
    # sigmay1 <- sigmay1 + coefflist1[[ind]]$sigma_y
    
  }
  coeff0 <- coeff0/(F_draw-indskip)
  coeff1 <- coeff1/(F_draw-indskip)
  randeff0 <- randeff0/(F_draw-indskip)
  randeff1 <- randeff1/(F_draw-indskip)
  # sigmay0 <- sigmay0/(F_draw-indskip)
  # sigmay1 <- sigmay1/(F_draw-indskip)
  
  coeflistres <- list(
    fixed0 = coeff0,
    fixed1 = coeff1,
    rand0 = randeff0,
    rand1 = randeff1#,
    # sigma_y0 = sigmay0,
    # sigma_y1 = sigmay1
  )

  # (4) return imputed values 
  res <- list(imputed0 = imputed_ctrl, imputed1= imputed_trt,
              pate_gen=mean(pategen_f, na.rm=T),
              wtdest = colMeans(wtdest_f, na.rm=T),
              # pate_aiptw = mean(pategen_aiptw_f, na.rm=T),
              varpategen = mean(varpategen_f,na.rm=T),
              pate_gen_synthpop = mean(pategen_synthpop_f, na.rm=T),
              pate_gen_synthpop_direct = mean(pategen_direct),
              varpategen_synthpop = mean(varpategen_synthpop_f, na.rm=T),
              coeflistres=coeflistres,
              knotlocations = knotlocations,
              phatsamp= phat_samp,
              phatpop = phat_pop,
              # time_aiptw = mean(time_aiptw),
              time_pencompwfpbb = mean(time_pencompwfpbb, na.rm=T)
  )
  }
  return(res)
}






