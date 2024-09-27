iter = 200
tdf <- 57-25
require(dplyr)

  aiptw_est <- NULL

  itersum <- 0
  for(i in 1:iter){
	  res <- readRDS(paste0("results/step1_res_faps_aiptwbootraw", i,".RDS"))
    aiptw_est <- c(aiptw_est, as.numeric(res$est))

  }
  
  aiptwres <- data.frame(
  est = mean(aiptw_est),
  se = sd(aiptw_est)
  )
  aiptwres <- aiptwres %>% mutate(
    cilower = est - qt(p=0.975, df = tdf)*se,
    ciupper = est + qt(p=0.975, df =tdf)*se
  )


saveRDS(aiptwres, file = paste0("results-clean/res_aiptwres.RDS"))
temp <- readRDS("results-clean/res_aiptwres.RDS")
