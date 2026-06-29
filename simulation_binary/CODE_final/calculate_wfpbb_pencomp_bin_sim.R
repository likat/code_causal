## FUNCTION FOR PENCOMP
calculate.wfpbb.pencomp <- function(propen_model, # string for propensity model, e.g. "A~X+Z"
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
  pategen_synthpop_f <- rep(NA, F_draw)
  varpategen_synthpop_f <- rep(NA, F_draw)
  coefflist0 <- vector(mode="list", length=F_draw)
  coefflist1 <- vector(mode="list", length=F_draw)
  knotloc0 <- vector(mode="list", length = F_draw)
  knotloc1 <- vector(mode="list", length = F_draw)
  pategen_direct <- rep(NA, F_draw)
  time_pencompwfpbb <- rep(NA,F_draw)
  # Step 1: Create a temporary directory
  temp_dir <- tempfile("cmdstan_temp_")
  dir.create(temp_dir)
  
  # Step 2: Save current output_dir option and point CmdStanR to temp
  # original_dir <- model$output_dir()
  original_dir <- tempdir()
  
  for(f in 1:F_draw){
    
    #-------------------------------
      # temp <- wtpolyap(ysamp = 1:nrow(bootdf), wts = bootwts, k = npop - nrow(bootdf) )
      temp <- wtpolyap(ysamp = uniqueind, wts = polyawts, k = npop - length(uniqueind) )
      
      # popbootpolya <- bootdf[temp,]
      popbootpolya <- samp[temp,] %>% dplyr::select(c("A", "Ysamp", "X1", "X2","Z"))
      # pategen_direct[f] <- mean(popbootpolya$Y1 - popbootpolya$Y0)
      # bt_polya <- popbootpolya[sample(1:npop, size = nrow(samp), replace=F),]
      
      ## PENCOMP
      t0_pencomp <- Sys.time()
      mod_trt <- glm(propen_model, data = popbootpolya, family="binomial")
      popbootpolya$phat <- predict(mod_trt, newdata=popbootpolya,type="link")
      phat_pop <- predict(mod_trt, newdata= popdat, type = "link")
      phat_samp <- predict(mod_trt, newdata = samp, type = "link")
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
      
        # if(trtgrp == 0){
        #   phat_trtgrp <- logit(1- expit(phat_trtgrp))
        # }
        
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
        # covariateX = cbind(model.matrix(outc_model, dfcopy_trtgrp))
        covariateX = cbind(model.matrix(outc_model, dfcopy_trtgrp), phat_trtgrp)
        
        ## Fit imputation model
        all=rep(1, nrow(dfcopy_trtgrp))
        if(Matrix::rankMatrix(covariateX)[1] < ncol(covariateX)){bootflag <<- 1; break}
        # if(sum(table(dfcopy_trtgrp$catecat) < 10) > 0){bootflag <<- 1; break}
        
        datastan <- list(
          N = nrow(covariateX),
          num_knot= num.knot,
          numpred = ncol(covariateX),
          y = response,
          pred_mat= covariateX,
          spline_mat= linearB
        )
        flag = 0
        numflag = 0
        # psppM <- run_cmdstan_temp(data = datastan, model = mod)
        # datastan$grainsize <- 2500
        datastan$grainsize <- 5375
        
        datastan$onevec <- rep(1, length(datastan$y))
        psppM <-
          mod0$optimize(
            data = datastan,
            jacobian = F,
            iter = 10000,
            show_messages=T,
            save_cmdstan_config=TRUE,
            threads = 2
          )
        microbenchmark(mod0$optimize(
          data = datastan,
          jacobian = F,
          iter = 10000,
          show_messages=T,
          save_cmdstan_config=TRUE,
          threads = 2
        ),times = 10L)
    
          if(psppM$return_codes() != 0){
            flag <- 1
            numflag = numflag + 1
            while(flag == 1){
              # psppM <- run_cmdstan_temp(data = datastan, model = mod)
              psppM <-
                mod$optimize(
                  data = datastan,
                  jacobian = F,
                  iter = 10000,
                  show_messages=F
                )
              if(psppM$return_codes() == 0){flag <- 0}
              numflag = numflag + 1
              if(numflag > 20){break}
            }
          }
          if(numflag > 20){return(NA); break}
          
        # if(psppM$return_codes() != 0){return(NA); break}
        
        ## COLLECT COEFFICIENTS FROM WFPBB
       
        ## imputations for synthetic population
          impute_potenoutsynthpop <-
            imputeF(
              newdata= popbootpolya[popbootpolya[,trt.varname]!=trtgrp,],
              # model = psppM,
              fixedcoef=psppM$summary("beta")$estimate,
              splinecoef = psppM$summary("gamma")$estimate,
              # sigmay = psppM$summary("sigma_y")$estimate,
              knotloc = knots,
              outcome.model =outc_model,
              # propen.score.new=logit(1-expit(popbootpolya$phat[popbootpolya[,trt.varname]!=trtgrp]))
              propen.score.new=popbootpolya$phat[popbootpolya[,trt.varname]!=trtgrp]

            )
        if(trtgrp == 0){
          imputed_ctrl_synthpop[popbootpolya[,trt.varname] == 0] <- 
            popbootpolya[popbootpolya[,trt.varname] == 0, y.varname]%>% unlist() %>% as.numeric()
          imputed_ctrl_synthpop[popbootpolya[,trt.varname] == 1] <- impute_potenoutsynthpop
          knotlocations[[1]] <- knots
          coefflist0[[f]] <- c(psppM$summary("beta")$estimate,
                               psppM$summary("gamma")$estimate)
          
          
        }else if(trtgrp ==1){
          imputed_trt_synthpop[popbootpolya[,trt.varname] == 1] <- 
            popbootpolya[popbootpolya[,trt.varname] == 1, y.varname] %>% 
            unlist %>% as.numeric
          imputed_trt_synthpop[popbootpolya[,trt.varname] == 0] <- impute_potenoutsynthpop
          knotlocations[[2]] <- knots
          coefflist1[[f]] <- c(psppM$summary("beta")$estimate,
                               psppM$summary("gamma")$estimate)
          
        }
      }
      # hist((imputed_ctrl_pop - pop$Y1)/mean(pop$Y1))
      t1_pencomp <- Sys.time()
      time_pencompwfpbb[f] <- t1_pencomp-t0_pencomp
      ## generalized PATE
      
      pategen_synthpop_f[f] <-
        mean(imputed_trt_synthpop - imputed_ctrl_synthpop)
      varpategen_synthpop_f[f] <- 
        var(imputed_trt_synthpop - imputed_ctrl_synthpop)/npop
      
      
    
  }
  
  
  # Get a list of all files in the directory
  # full.names = TRUE ensures the full path to each file is returned
  files_to_delete <- list.files(path = original_dir, full.names = TRUE)
  file.remove(files_to_delete)
  # print(f)
  unlink(temp_dir, recursive = TRUE)
  
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
  res <- list(
              pate_gen_synthpop = mean(pategen_synthpop_f, na.rm=T),
              # pate_gen_synthpop_direct = mean(pategen_direct),
              varpategen_synthpop = mean(varpategen_synthpop_f, na.rm=T),
              coeflistres=coeflistres,
              knotlocations = knotlocations,
              # time_aiptw = mean(time_aiptw),
              time_pencompwfpbb = mean(time_pencompwfpbb, na.rm=T)
  )
  return(res)
}






