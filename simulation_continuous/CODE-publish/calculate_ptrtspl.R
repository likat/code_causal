impute.ptrt=function(newdata, 
                 propen_mod,
                 fixedcoef,
                 splinecoef,
                 knotloc = knots,
                sel.score.new,
                poptrue = F) {
  
  linearB=NULL
  linearB =outer(sel.score.new, knotloc, "-")
  linearB =linearB * (linearB > 0)

  if(poptrue == T){
  designM = cbind(model.matrix(as.formula( paste0(as.character(propen_mod)[1], as.character(propen_mod)[3])), newdata), sel.score.new, linearB)
  }else{
    designM = cbind(model.matrix(propen_mod, newdata), sel.score.new, linearB)
  }
  
  designM=as.matrix(designM)
  fittedcoef <- c(fixedcoef, splinecoef)
  predicted = designM %*% fittedcoef
  
  return(predicted)
  # return(designM)
  
}
calculate.ptrtspline <- function(fitdat,
                                 sampdat,
                                 popdat2,
                                 num.knot,
                                 propen.model = propen_model,
                                 trtname){
  # (3) Fit outcome model with splines
  ## need sensitivity analysis for how many knots
  ## currently borrowing code from zhou 2019
  psel <- 1/fitdat$wts 
  ## sometimes psel is out of bounds when estimating this
  ## should we just use raw psel, no logit?
  # psellogit <- log(psel/(1-psel))
  psellogit <- log(psel)
  # psellogit <- psel
  space=(max(psellogit)-min(psellogit))/(num.knot+1)
  knots=(min(psellogit)+space*(1:num.knot))
  
  ## Construct C2 matrix, truncated linear basis splines
  linearB=NULL
  linearB =outer(psellogit, knots, "-") ## constructing C2 matrix
  linearB =linearB * (linearB > 0)
  
  ## Build C1 matrix, define outcome
  response=fitdat[, trtname] %>% unlist %>% as.numeric
  covariateX = cbind(model.matrix(propen.model, fitdat), psellogit)
  
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
    iter=10000
  )
  ## save imputations for the other treatment group
  imputeptrt_synthpop <- 
    impute.ptrt(
      newdata=fitdat,
      propen_mod = propen.model,
      fixedcoef = psppM$summary("beta")$estimate,
      splinecoef = psppM$summary("gamma")$estimate,
      knotloc = knots,
      sel.score.new=log(1/fitdat$wts)
      )
  imputeptrt_samp <- 
    impute.ptrt(
      newdata=sampdat,
      propen_mod = propen.model,
      fixedcoef = psppM$summary("beta")$estimate,
      splinecoef = psppM$summary("gamma")$estimate,
      knotloc = knots,
      sel.score.new=log(1/sampdat$wts),
      poptrue = F
      )
  
  imputeptrt_pop <- 
    impute.ptrt(
      newdata=popdat2,
      propen_mod = propen.model,
      fixedcoef = psppM$summary("beta")$estimate,
      splinecoef = psppM$summary("gamma")$estimate,
      knotloc = knots,
      sel.score.new=log(1/popdat2$wts),
      poptrue = T
    )
  
  reslist <- list(
    ptrt_synthpop = imputeptrt_synthpop,
    ptrt_samp = imputeptrt_samp,
    ptrt_pop = imputeptrt_pop,
    fixedcoef = psppM$summary("beta")$estimate,
    splinecoef = psppM$summary("gamma")$estimate
  )
  return(reslist)
}
