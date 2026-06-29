numiter = 200
type_outcome = "gb"
type_overlap = "notshared"

## find bias, rmse, coverage
get.summary <- function(x, truex){
  ## extra entry due to the way we read in the rawdata
  x <- x[-1,]
  bias <- mean(x$est - truex)
  rmse <- sqrt(mean((x$est - truex)^2))
  cilength <- mean(x$ciupper - x$cilower, na.rm=T)
  coverage = mean(truex >= x$cilower & truex <= x$ciupper)
  resdf <- data.frame(
    bias = bias,
    rmse = rmse,
    cilength = cilength,
    coverage = coverage
  )
  return(resdf)
}

## overall function for getting stuff
get.overall.aipw <- function(numiter=200, type_outcome, type_overlap){
  tempres <- 
    readRDS(paste0("../results-raw-2stage/aipw/results_tingting_aipw", type_outcome, type_overlap, "1.RDS"))
ind_estonly <- (2:(length(tempres)-3))
reslist <- vector(mode = "list", length = length(tempres )-2)
reslist <- tempres[ind_estonly]
time_aipw <- 0
# get raw results from all simulations
for(i in 1:numiter){
  count = 1
  tempres <- readRDS(paste0(
    "../results-raw-2stage/aipw/results_tingting_aipw", type_outcome, type_overlap, i, ".RDS"
  ))
  time_aipw <- time_aipw + tempres$time_aipw/numiter
  for(j in ind_estonly){
    reslist[[count]] <- rbind(reslist[[count]], tempres[[j]])
    count = count + 1
  }
}
# get summary
df_summary <- vector(mode = "list", length = length(reslist))
names(df_summary) <- names(reslist)
for(g in 1:length(reslist)){
  df_summary[[g]] <- get.summary(x = reslist[[g]], truex = tempres$patetrue)
}
df_summary$time = time_aipw
resfinal <- list(raw = reslist, overall = df_summary, time = time_aipw)
return(resfinal)
}

# get.overall.aipw(numiter = 200, type_outcome = "gb", type_overlap = "notshared")
outcometypes <- c("nogb", "gb")
overlaptypes <- c("notshared", "shared")
refdf <- data.frame(
  outcometypes = rep(outcometypes, 2),
  overlaptypes = rep(overlaptypes, each = 2),
  numiter = c(200)
)
for(j in 1:nrow(refdf)){
  outcomecurrent <- refdf$outcometypes[j]
  overlapcurrent <- refdf$overlaptypes[j]
  restemp <- get.overall.aipw(type_outcome = outcomecurrent, type_overlap = overlapcurrent)
  saveRDS(restemp$raw, file = paste0("../results-clean/res_",outcomecurrent,overlapcurrent,"_aipw_raw.RDS"))
  saveRDS(restemp$overall, file = paste0("../results-clean/res_",outcomecurrent,overlapcurrent,"_aipw.RDS"))
  
  print(j)
}
readRDS(paste0("../results-clean/res_",outcomecurrent,overlapcurrent,"_aipw.RDS")
)

