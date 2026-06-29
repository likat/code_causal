
calculate.sate <- function(dat){
  temp_est <- mean(dat$Ysamp[dat$A==1], na.rm=T) - mean(dat$Ysamp[dat$A==0], na.rm=T)
  return(temp_est)
}
