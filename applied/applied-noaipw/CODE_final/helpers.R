# library loading
# Libraries -------
library(dplyr)
library(stringr)
# require(rstan)
require(polyapost)
# require(LaplacesDemon)

## COMBINING FUNCTION FOR AIPTW
combine.aiptw <- function(boottbl, num_boot = numboot){
  est <-  mean(boottbl, na.rm=T)
  varest <- var(boottbl, na.rm=T)
  cilower <- est - 1.96*sqrt(varest)
  ciupper <- est + 1.96*sqrt(varest)
  res <- list(
    est = est,
    var = varest,
    cilower = cilower,
    ciupper = ciupper
  )
}


## RUBIN COMBINING FUNCTION-----
combine.rubin.applied <- function(boottbl, num_boot=numboot, matrix = FALSE){

  est <- mean(boottbl[,1], na.rm=T)
  withinvar <- mean(boottbl[,2], na.rm=T)
  btwnvar <- var(boottbl[,1], na.rm=T)
  mivar <- withinvar + (1+1/num_boot)*btwnvar
  tdf <- round((num_boot-1)*(1+withinvar/((num_boot+1)*btwnvar)^2))
  cilower_t <- est - qt(0.975,tdf)*sqrt(mivar)
  ciupper_t <- est + qt(0.975,tdf)*sqrt(mivar)
  cilower_emp <- quantile(boottbl[,1],c(0.025, 0.975), na.rm=T)[1]
  ciupper_emp <- quantile(boottbl[,1],c(0.025, 0.975), na.rm=T)[2]
  res <- list(est=est, withinvar=withinvar,btwnvar=btwnvar,mivar=mivar,
              ci_tlength = ciupper_t - cilower_t,
              ci_emplength = as.numeric(ciupper_emp- cilower_emp))
  
  return(res)
}

combine.rubin <- function(boottbl, num_boot=numboot){
  est <- mean(boottbl[,1], na.rm=T)
  withinvar <- mean(boottbl[,2], na.rm=T)
  btwnvar <- var(boottbl[,1], na.rm=T)
  mivar <- withinvar + (1+1/num_boot)*btwnvar
  tdf <- round((num_boot-1)*(1+withinvar/((num_boot+1)*btwnvar)^2))
  cilower_t <- est - qt(0.975,tdf)*sqrt(mivar)
  ciupper_t <- est + qt(0.975,tdf)*sqrt(mivar)
  cilower_emp <- quantile(boottbl[,1],c(0.025, 0.975), na.rm=T)[1]
  ciupper_emp <- quantile(boottbl[,1],c(0.025, 0.975), na.rm=T)[2]
  cicov_t <- between(pate_true, cilower_t, ciupper_t)
  cicov_emp <- between(pate_true, cilower_emp, ciupper_emp)
  res <- list(est=est, withinvar=withinvar,btwnvar=btwnvar,mivar=mivar,
              ci_tlength = ciupper_t - cilower_t,
              ci_emplength = as.numeric(ciupper_emp- cilower_emp),
              cicov_t = cicov_t, cicov_emp = cicov_emp)
  return(res)
}

## HL statistic
get.hl <- function(obs, pred, prob){
  tempdf <- data.frame(obs = obs, pred = pred, prob = prob) %>% 
    mutate(quant = ntile(prob, 10))
  tempdfsum <- tempdf %>% group_by(quant) %>% summarize(
    nq = n(),
    Oq = sum(obs),
    Eq = sum(pred)
  )
  tempdfsum = tempdfsum %>% mutate(hlstat = (Oq - Eq)^2/(Eq*(1-Eq/nq)))
  hlstat = sum(tempdfsum)
  hlstat_pval = 1-pchisq(q=hlstat, df = 8)
  res <- data.frame(stat = hlstat, pval = hlstat_pval)
  return(res)
}

## CONCATENATE LISTS WITH ELEMENTS WITH SAME NAME ----
# from here: https://stackoverflow.com/questions/18538977/combine-merge-lists-by-elements-names
library(purrr)

cat_lists <- function(list1, list2) {  
  
  keys <- unique(c(names(list1), names(list2)))
  map2(list1[keys], list2[keys], c) %>% 
    set_names(keys)  
  
}

# combined_output <- reduce(input_list, cat_lists) ## USAGE


# CONTAINER CREATION FUNCTION------
# purpose: create containers to store results in
# INPUTS:
#   basename <string> = what is the base name of the container?
#   subgroups <real> = number of subgroups we are estimating
#   type <string> = matrix or vector?
#   dim <real> or <vector> = <real> if type is vector, <vector>[2] if type is matrix
# OUTPUT 
#   assign subgroups+1 objects to the global environment
# USAGE
# containerfun(basename="y_bar", subgroups=grps, type="vector",dims = sim)
# containerfun(basename="CI", subgroups=grps, type="matrix",dims = c(sim,2))

containerfun <- function(basename, subgroups=5, type="vector", dims){
  if(type=="vector"){
    assign(basename, rep(NA,dims), envir = .GlobalEnv)
    for(it in 1:subgroups){
      assign(paste(basename,it, sep=""), rep(NA,dims), envir=.GlobalEnv)
    }
  }	
  if(type=="matrix"){
    assign(basename, matrix(nrow = dims[1], ncol=dims[2]), envir = .GlobalEnv)
    for(it in 1:subgroups){
      assign(paste(basename,it, sep=""), matrix(nrow = dims[1], ncol=dims[2]), envir=.GlobalEnv)
    }
  }	
}
formulaF=function(varList, y.name){
  return ( as.formula(paste(y.name, "~ ", paste(c(varList), collapse = "+"))) )
}

# EXPIT ------
# purpose: expit function
expit <- function(x){
  return(exp(x)/(1+exp(x)))
}

# DLL unloading -----
## Workaround for fitting many stan models in the cluster
# to fix the dyn.load() unable to load shared object warning message
# adapted from https://github.com/stan-dev/rstan/issues/448
unload.dll <- function(model=model){
  dso_filename = model@dso@dso_filename
  loaded_dlls = getLoadedDLLs()
  if (dso_filename %in% names(loaded_dlls)) {
    message("Unloading DLL for model dso ", dso_filename)
    model.dll = loaded_dlls[[dso_filename]][['path']]
    dyn.unload(model.dll)
  } else {
    message("No loaded DLL for model dso ", dso_filename)
  }
  
  loaded_dlls = getLoadedDLLs()
  loaded_dlls = loaded_dlls[str_detect(names(loaded_dlls), '^file')]
  if (length(loaded_dlls) > 10) {
    for (dll in head(loaded_dlls, -10)) {
      message("Unloading DLL ", dll[['name']], ": ", dll[['path']])
      dyn.unload(dll[['path']])
    }
  }
  message("DLL Count = ", length(getLoadedDLLs()), ": [", str_c(names(loaded_dlls), collapse = ","), "]")
  
}

# NORMALIZE -------
# purpose: normalize a vector to sum to 1
# input: double vector[n]
# output: normalized vector[n]
normalize <- function(vec){
  return(vec/sum(vec))
}

# NOT IN--------
# purpose: boolean operator; are elements in 'vector' not in 'othervector'?
# use: double vector[n] %notin% another double vector[n]
# output: T or F 
`%notin%` <- Negate(`%in%`)



