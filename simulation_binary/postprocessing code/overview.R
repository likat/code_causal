require(dplyr)
## Main simulation results
gbshared <- readRDS("res_gbshared.RDS")
gbnotshared <- readRDS("res_gbnotshared.RDS")
nogbnotshared <- readRDS("res_nogbnotshared.RDS")
nogbshared <- readRDS("res_nogbshared.RDS")

## AIPW simulation results
gbshared_aipw <- readRDS("res_gbshared_aipw_bin.RDS")
gbnotshared_aipw <- readRDS("res_gbnotshared_aipw_bin.RDS")
nogbnotshared_aipw <- readRDS("res_nogbnotshared_aipw_bin.RDS")
nogbshared_aipw <- readRDS("res_nogbshared_aipw_bin.RDS")


arrange.res <- function(resdfname, resdfname.aipw){
  resdf <- get(resdfname)
  resdf.aipw <- get(resdfname.aipw)
  restab <- data.frame(
    namesim = resdfname,
    estimator = c("ATE Oracle Model Impute",
                  "ATE Incorrect Model Impute",
                  "ATE PENCOMP",
                  "WATE Oracle Model Impute",
                  "WATE Incorrect Model Impute",
                  "WATE PENCOMP",
                  "SYNGATE WFPBB-PENCOMP",
                  "selwt-AIPW"),
    absbias = abs(c(resdf$summary_ate_nospline$bias,
                 resdf$summary_ate_nospline_wrongmod$bias,
                 resdf$summary_ate_pencomp$bias,
                 resdf$summary_wtdest_nospline$bias,
                 resdf$summary_wtdest_nospline_wrongmod$bias,
                 resdf$summary_wtdest_pencomp$bias,
                 resdf$summary_synthpop_pencompwfpbb$bias,
                 resdf.aipw$aipw_unwtd_correcttrt$bias
    )) %>% round(3),
    rmse = c(resdf$summary_ate_nospline$rmse,
             resdf$summary_ate_nospline_wrongmod$rmse,
             resdf$summary_ate_pencomp$rmse,
             resdf$summary_wtdest_nospline$rmse,
             resdf$summary_wtdest_nospline_wrongmod$rmse,
             resdf$summary_wtdest_pencomp$rmse,
             resdf$summary_synthpop_pencompwfpbb$rmse,
             resdf.aipw$aipw_unwtd_correcttrt$rmse
    )%>% round(3),
    cilength = c(resdf$summary_ate_nospline$cilength,
                 resdf$summary_ate_nospline_wrongmod$cilength,
                 resdf$summary_ate_pencomp$cilength,
                 resdf$summary_wtdest_nospline$cilength,
                 resdf$summary_wtdest_nospline_wrongmod$cilength,
                 resdf$summary_wtdest_pencomp$cilength,
                 resdf$summary_synthpop_pencompwfpbb$cilength,
                 resdf.aipw$aipw_unwtd_correcttrt$cilength
    )%>% round(3),
    cicov = c(resdf$summary_ate_nospline$cicov,
              resdf$summary_ate_nospline_wrongmod$cicov,
              resdf$summary_ate_pencomp$cicov,
              resdf$summary_wtdest_nospline$cicov,
              resdf$summary_wtdest_nospline_wrongmod$cicov,
              resdf$summary_wtdest_pencomp$cicov,
              resdf$summary_synthpop_pencompwfpbb$cicov,
              resdf.aipw$aipw_unwtd_correcttrt$coverage
    ),
    patetrue = resdf$patetrue,
    time_per_synthpop_wfpbb_seconds = resdf$time_pencompwfpbb %>% as.numeric,
    time_per_bootstrap_aipw = resdf.aipw$time %>% as.numeric
  )
}
temp1 <- arrange.res("gbshared", resdfname.aipw = "gbshared_aipw")
temp2 <- arrange.res("nogbshared", resdfname.aipw = "nogbshared_aipw")
temp3 <- arrange.res("gbnotshared", resdfname.aipw = "gbnotshared_aipw")
temp4 <- arrange.res("nogbnotshared", resdfname.aipw = "nogbnotshared_aipw")

results_all <- dplyr::bind_rows(temp1, temp2, temp3, temp4)
saveRDS(results_all, file = "results_all_bin.RDS")
write.csv(results_all, file = "results_all_bin.csv")
