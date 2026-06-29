calculate.sate.x <- function(dat, xname){
  temptbl1 <- dat %>% filter(A==1) %>% group_by(get(xname)) %>% summarise(y1=mean(Ysamp))
  temptbl0 <- dat %>% filter(A==0) %>% group_by(get(xname)) %>% summarise(y0=mean(Ysamp))
  return(temptbl1$y1-temptbl0$y0)
}
