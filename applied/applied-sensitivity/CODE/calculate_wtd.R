calculate.wtd <- function(imputedlist, bootdf = boot_all){
  bootdf$ydiff <- imputedlist[[2]]-imputedlist[[1]]
  des <- svydesign(id=~0, weights=~wts,data=bootdf)
  patesvymean <- svymean(~ydiff, design=des) 
  wtdest <- patesvymean %>% as.numeric()
  ydiff_var <- SE(patesvymean)^2
  # catesvymean <- svyby(~ydiff, ~catecat, des, svymean)
  # cate <- catesvymean$ydiff
  # cate_var <- catesvymean$se^2
  # ydiff <- imputedlist[[2]]-imputedlist[[1]]
  # wtdest <- sum(ydiff*bootdf$wts/sum(boot_all$wts))
  # ydiff_var <- (1/sum(bootdf$wts))^2*sum(bootdf$wts^2*var(bootdf$ydiff))
  # return(list(wtdest,ydiff_var, cate=cate, catevar=cate_var))
  return(list(wtdest,ydiff_var))
  
}