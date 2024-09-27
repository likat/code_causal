# library(tidyverse)
# library(tidycensus)
# library(foreign)
set.seed(42)
expit <- function(x){
  return(exp(x)/(1+exp(x)))
}

N <- 66000 ## approx size of the ACS data
pop <- data.frame(
  X1 = rnorm(N, mean=0, sd=1),
  X2 = rnorm(N, mean=0, sd=1),
  logZ = rnorm(N,mean=0, sd=0.7)
) %>% mutate(Z = exp(logZ))


## TREATMENT ASSIGNMENT COEFFICIENTS
a0 <- 0
ax1 <- 1
ax2 <- 1
ax1x2 <- 0.5
# ax <- c(0, -0.5, 0.35,0.75)
# az <- -0.75
# azx1 <- -0.5
az <- -0.75
azx1 <- 0

pop$ptrt <-  expit(a0 + ax1*pop$X1 + ax2*pop$X2 + ax1x2*pop$X1*pop$X2 +
                     az*pop$Z + azx1*pop$Z*pop$X1)

## OUTCOME MODEL COEFFICIENTS
# We need azx term to necessitate the PENCOMP.
set.seed(50)
b0 <- 0
# b1 <- 2 #A
# b2 <- c(0,-1.5,1,0.5)#c(0, runif(n=3, min = -1, max = 2)) # X
bx1 <- 1 # X1
bx2 <- 1 # X2
bA <- 5
bAx1 <- 2
bAx2 <- 0
bz <- 4
bAz <- 2
# bz <- 0
# bAz <- 0
# bAz <- 0
# bax1z <- 0.8
bax1z <- 0
bax2z <- 0.8
# bax2z <- 0


# b3 <- 1 #log(hhwt)Z
# b4 <- c(0, -0.5,0.4, -0.8) # xz
# b5 <- c(0,-0.75, 0.5, -0.25)#c(0,runif(n=3, min = -1, max = 1)) # AX
# b6 <- 1 #AZ
# b7 <- c(0,-0.75, 0.5, -0.25)#c(0,runif(n=3, min = -1, max = 1)) #azx

## SELECTION PROBABILITY
# Sample size
nsamp <- round(0.1*N)
pop$pselect <- sampling::inclusionprobabilities(pop$Z, nsamp)
quant_psel <- quantile(pop$pselect)
pop <- pop %>% mutate(catecat = case_when(
  pselect < quant_psel[2] ~ 1,
  pselect < quant_psel[3] ~ 2,
  pselect < quant_psel[4] ~ 3,
  TRUE~ 4
))
pop$wts <- 1/pop$pselect



#### 
##======= BENCHMARK JOINT #####
####
## How is treatment assigned? First work with the simplest case of fully randomized treatment
