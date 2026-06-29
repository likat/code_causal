iter = 200
tdf <- 57-25
require(dplyr)
tempres <- readRDS("results/step1_res_faps_aipw_new.RDS")

resdf <- data.frame(
  est = tempres$est_wtd_pt,
  sd = sd(tempres$est_wtd)
) %>% mutate(
  cilower = est - qt(0.975, tdf) * sd,
  ciupper = est + qt(0.975, tdf) * sd
)
mean(tempres$est_wtd_pt)
sd(tempres$est_wtd)
saveRDS(resdf, file = "results-clean/res_aipw_new.RDS")
