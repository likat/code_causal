imputeF=function(newdata, 
                 model, 
                 # x.varnames, 
                 knotloc = knots,
                 outcome.model,
                 propen.score.new) {
  
  # knots=model$knot.loc
  # space=(max(propen.score.new)-min(propen.score.new))/(num.knot+1)
  # knots=(min(propen.score.new)+space*(1:num.knot))
  
  linearB=NULL
  linearB =outer(propen.score.new, knotloc, "-")
  linearB =linearB * (linearB > 0)
  
  designM = cbind(model.matrix(outcome.model, newdata), propen.score.new)
  # designM = cbind(model.matrix(outcome.model, newdata))
  designM=as.matrix(designM)
  fittedcoef <- model$coefficients$fixed
  randomcoef <- model$coefficients$random$all 
  # predicted = designM %*% model$coefficients$fixed + as.matrix(linearB) %*% as.vector(unlist(model$random)) + rnorm(nrow(newdata), 0, model$sigma)
  predicted = designM %*% fittedcoef + as.matrix(linearB) %*% as.numeric(randomcoef) + rnorm(nrow(newdata), 0, model$sigma)
  
  return(predicted)
  
}

# imputeF.sel=function(newdata, 
#                  model, 
#                  # x.varnames, 
#                  knotloc = knots,
#                  outcome.model,
#                  propen.score.new,
#                  sel.score.new) {
#   
#   # knots=model$knot.loc
#   # space=(max(propen.score.new)-min(propen.score.new))/(num.knot+1)
#   # knots=(min(propen.score.new)+space*(1:num.knot))
#   
#   linearB=NULL
#   linearB =outer(propen.score.new, knotloc, "-")
#   linearB =linearB * (linearB > 0)
#   
#  
#   designM = cbind(model.matrix(outcome.model, newdata), propen.score.new, sel.score.new)
#   # designM = cbind(model.matrix(outcome.model, newdata))
#   designM=as.matrix(designM)
#   fittedcoef <- model$coefficients$fixed
#   randomcoef <- model$coefficients$random$all 
#   # predicted = designM %*% model$coefficients$fixed + as.matrix(linearB) %*% as.vector(unlist(model$random)) + rnorm(nrow(newdata), 0, model$sigma)
#   predicted = designM %*% fittedcoef + as.matrix(linearB) %*% as.numeric(randomcoef) + rnorm(nrow(newdata), 0, model$sigma)
#   
#   return(predicted)
#   
# }