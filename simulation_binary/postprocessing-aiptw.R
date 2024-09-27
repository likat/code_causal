iter = 200
N=66000

get.summary <- function(x, truex){
  # sampsize <- readRDS(paste0("../results-clean/sampsize_",type_outcome, type_overlap,".RDS"))
  # fpc <- (N - sampsize) / (N-1)
  
    bias <- mean(x$est - truex)
    rmse <- sqrt(mean((x$est - truex)^2))
    cilower <- x$est - qt(0.975, df = 20-1)*sqrt(x$var)
    ciupper <- x$est + qt(0.975, df = 20-1)*sqrt(x$var)
    cilength <- mean(ciupper - cilower)
    cicov_t <- mean(truex >= cilower & truex <= ciupper)
    
    res <- list(
    bias = bias,
    rmse = rmse,
    cilength = cilength,
    cicov_t = cicov_t
    )
return(res)
}
resnames <- c("est", "var", "cilower", "ciupper")
getres <- function(type_outcome, type_overlap){
  aiptw_est <- vector(mode="list", length=4)
  patetrue <- readRDS(paste0("../results-raw/results_tingting_aiptw", type_outcome, type_overlap,"1.RDS"))
patetrue <- patetrue$patetrue
  
  itersum <- 0
  for(i in 1:iter){
	  res <- readRDS(paste0("../results-raw/results_tingting_aiptw", type_outcome, type_overlap,i,".RDS"))

    aiptw_est <- Map(c, aiptw_est, res$aiptw)

  }
  names(aiptw_est) <- resnames
  aiptwres <- get.summary(aiptw_est, truex = patetrue)
  reswrite <- list(
    patetrue = patetrue,
aiptw_est =aiptwres
)
return(reswrite)
}

outcometypes <- c("nogb", "gb")
overlaptypes <- c("notshared", "shared")

for(outcome in 1:length(outcometypes)){
  for(overlap in 1:length(overlaptypes)){
    outcomecurrent <- outcometypes[outcome]
    overlapcurrent <- overlaptypes[overlap]
    getres(type_outcome = outcomecurrent, type_overlap = overlapcurrent)
saveRDS(getres(type_outcome = outcomecurrent, type_overlap = overlapcurrent), file = paste0("../results-clean/res_aiptw",outcomecurrent,overlapcurrent,".RDS"))
    
  }
}
readRDS("../results-clean/res_aiptwgbnotshared.RDS")
readRDS("../results-clean/res_aiptwnogbnotshared.RDS")
readRDS("../results-clean/res_aiptwnogbshared.RDS")
readRDS("../results-clean/res_aiptwgbshared.RDS")
