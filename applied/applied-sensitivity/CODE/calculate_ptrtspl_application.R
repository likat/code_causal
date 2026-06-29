impute.ptrt=function(newdata, 
                 model, 
                 # x.varnames, 
                 knotloc = knots,
                 outcome.model,
                sel.score.new) {
  
  # knots=model$knot.loc
  # space=(max(propen.score.new)-min(propen.score.new))/(num.knot+1)
  # knots=(min(propen.score.new)+space*(1:num.knot))
  
  linearB=NULL
  linearB =outer(sel.score.new, knotloc, "-")
  linearB =linearB * (linearB > 0)
  
  # designM = cbind(model.matrix(outcome.model, newdata), propen.score.new)
  designM = cbind(model.matrix(outcome.model, newdata), linearB)
  designM=as.matrix(designM)
  fittedcoef <- model$coefficients
  # randomcoef <- model$coefficients$random$all
  # predicted = designM %*% fittedcoef + as.matrix(linearB) %*% as.numeric(randomcoef) + rnorm(nrow(newdata), 0, model$sigma)
  predicted = designM %*% fittedcoef
  
  return(predicted)
  
}
calculate.ptrtspline <- function(fitdat,
                                 sampdat,
                                 num.knot = 5,
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
  covariateX = cbind(model.matrix(propen.model, fitdat), linearB)
  # colnames(covariateX) <- c("int", "X1",
  #                     paste0(rep("spl", num.knot), 1:num.knot))
  ## Fit imputation model
  all=rep(1, nrow(fitdat))
  ## SINGULARITY ISSUES HERE
  if(Matrix::rankMatrix(covariateX)[1] < ncol(covariateX)){bootflag <<- 1; break}
  ptrtmod <- glm(response~covariateX-1, family = "binomial")
 temp <-  predict(ptrtmod, fitdat, type = "link")

  ## save imputations for the other treatment group
  imputeptrt_boot <- 
    impute.ptrt(newdata=fitdat,
            model = ptrtmod,
            knotloc = knots,
            outcome.model =propen.model,
            sel.score.new=log(1/fitdat$wts))
  imputeptrt_samp <- 
    impute.ptrt(newdata=sampdat,
                model = ptrtmod,
                knotloc = knots,
                outcome.model =propen.model,
                sel.score.new=log(1/sampdat$wts))
  reslist <- list(
    ptrt_fitdat = imputeptrt_boot,
    ptrt_samp = imputeptrt_samp
  )
  return(reslist)
}
