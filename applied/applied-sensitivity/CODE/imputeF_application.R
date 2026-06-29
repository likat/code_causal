## FUNCTION TO IMPUTE 
imputeF=function(newdata, 
                 model, 
                 knotloc = knots,
                 outcome.model,
                 propen.score.new) {
  
  linearB=NULL
  linearB =outer(propen.score.new, knotloc, "-")
  linearB =linearB * (linearB > 0)
  
  designM = cbind(model.matrix(outcome.model, newdata), linearB)
  designM=as.matrix(designM)
  fittedcoef <- model$coefficients
  predicted = designM %*% fittedcoef
  ## change back to probability scale
  predicted = expit(predicted)
  return(predicted)
  
}


