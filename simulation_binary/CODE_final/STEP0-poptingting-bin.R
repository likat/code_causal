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

# Assign cluster (100 clusters, sample 20)
numclus <- 400
numclussamp <- 200
# numclussamp <- 10
nstrat = 1
numclussamp_perstrat <- numclussamp / nstrat
pop$cluster <- rep(1:numclus, each= N/numclus)
## cluster specific random effect
cluseff <- rnorm(numclus)
pop$cluseff <- cluseff[pop$cluster]
pop$ptrt <-  expit(a0 + ax1*pop$X1 + ax2*pop$X2 + ax1x2*pop$X1*pop$X2 +
                     az*pop$Z + azx1*pop$Z*pop$X1)
pop <- pop %>% mutate(
  # stratum = case_when(
  #   cluster <= 100 ~ 1,
  #   cluster <= 200 ~ 2
  # )
  stratum = 1
  
)
## OUTCOME MODEL COEFFICIENTS
# We need azx term to necessitate the PENCOMP.
set.seed(50)
b0 <- 0
bx1 <- 0.2 # X1
bx2 <- 0.1 # X2
bA <- 0.5
# bAx1 <- 0.8
# bAx2 <- 0.5
bAx1 <- 0.8
bAx2 <- 0.5

bz <- -0.2
bAz <- 0.4


# b3 <- 1 #log(hhwt)Z
# b4 <- c(0, -0.5,0.4, -0.8) # xz
# b5 <- c(0,-0.75, 0.5, -0.25)#c(0,runif(n=3, min = -1, max = 1)) # AX
# b6 <- 1 #AZ
# b7 <- c(0,-0.75, 0.5, -0.25)#c(0,runif(n=3, min = -1, max = 1)) #azx

## SELECTION PROBABILITY
# Sample size
nsamp <- round(0.2*N)
pop$pselect1stage <- sampling::inclusionprobabilities(pop$Z, round(0.1*N))
temp <- sampling::inclusionprobabilities(pop$Z, nsamp)
pop$pselectf2 <- temp
pop$pselectf1f2 <- numclussamp_perstrat / (numclus/nstrat) * pop$pselectf2

# quant_psel <- quantile(pop$pselect)
# pop <- pop %>% mutate(catecat = case_when(
#   pselect < quant_psel[2] ~ 1,
#   pselect < quant_psel[3] ~ 2,
#   pselect < quant_psel[4] ~ 3,
#   TRUE~ 4
# ))
pop$wts1stage <- 1/pop$pselect1stage
pop$wts2stage <- 1/pop$pselectf1f2

saveRDS(pop, file = "simpop-bin.RDS")
