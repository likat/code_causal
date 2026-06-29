## read in cleaned results
noaiptw <- readRDS(paste0("results-clean/res_noaiptw_application_numknot15.RDS"))
aiptw <- readRDS("../causal-application/results-clean/res_aipw_new.RDS")
naive <- readRDS("results-clean/step1_res_faps_naive.RDS")
## Assemble table

resdf <- data.frame(
  name = c(
    "Naive model-based estimator",
    "Naive HT estimator",
    "ATE Model imputation",
    "Weighted PENCOMP",
    "WFPBB-PENCOMP",
    "AIPW"
  ),
  est = c(
    naive$est_naivemod,
    naive$est_ht,
    noaiptw$ate_linimipute$est,
          # noaiptw$ate_pencomp$est,
          # noaiptw$ate_ptrtzonly$est,
          # noaiptw$wtd_linimpute$est,
          noaiptw$wtd_pencomp$est,
          # noaiptw$wtd_ptrtzonly$est,
          noaiptw$synthpop_ptrtzonly$est,
          aiptw$est),
  sd = c(naive$se_mod,
         naive$se_ht,
          noaiptw$ate_linimipute$mivar%>% sqrt,
          # noaiptw$wtd_linimpute$mivar%>% sqrt,
          noaiptw$wtd_pencomp$mivar%>% sqrt,
          # noaiptw$wtd_ptrtzonly$mivar%>% sqrt,
          noaiptw$synthpop_ptrtzonly$mivar %>% sqrt,
          aiptw$sd),
  cilower =c(naive$cilower_naivemod,
             naive$cilower_ht,
             noaiptw$ate_linimipute$cilower,
             # noaiptw$wtd_linimpute$mivar%>% sqrt,
             noaiptw$wtd_pencomp$cilower,
             # noaiptw$wtd_ptrtzonly$mivar%>% sqrt,
             noaiptw$synthpop_ptrtzonly$cilower,
             aiptw$cilower),
  ciupper =c(naive$ciupper_naivemod,
             naive$ciupper_ht,
             noaiptw$ate_linimipute$ciupper,
             # noaiptw$wtd_linimpute$mivar%>% sqrt,
             noaiptw$wtd_pencomp$ciupper,
             # noaiptw$wtd_ptrtzonly$mivar%>% sqrt,
             noaiptw$synthpop_ptrtzonly$ciupper,
             aiptw$ciupper)
)

saveRDS(resdf, file = "results-clean/res_all_applied.RDS")
