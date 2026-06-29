## FUNCTION TO IMPUTE 
imputeF=function(newdata, 
                 # model, 
                 fixedcoef,
                 splinecoef,
                 knotloc = knots,
                 outcome.model,
                 propen.score.new) {
  
  linearB=NULL
  linearB =outer(as.numeric(propen.score.new), knotloc, "-")
  linearB =linearB * (linearB > 0)
  
  designM = cbind(model.matrix(outcome.model, newdata), propen.score.new,linearB)
  designM=as.matrix(designM)
  fittedcoef <- c(fixedcoef, splinecoef)
  predictedlogit = designM %*% fittedcoef
  ## change back to probability scale
  predicted = rbinom(length(predictedlogit), size = 1, prob = expit(predictedlogit))
  res <- data.frame(
    predicted = predicted,
    predictedprob = expit(predictedlogit)
  )
  # return(predicted)
return(res)
}


