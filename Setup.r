
library(mvtnorm)
library(Matrix)
set.seed(3809)

source("SdeBM.r")
source("SolveSdeBM.r")
source('RunKF.r')


source('Data_sim.r')
#source('Data_real.r')

mxsamp = max(all_t)
alen = max(dim(ax))
all_state <- all_s
all_time <- as.matrix( all_t )
nai = 5 # 5 animals


# ## initial value for simulated dataset
# ## diffusion parameters
 alpha <- 1.2 #0.8
 beta = 0
 rho   <- 5  #5.2
 sigma <- 0.7   # 2.7
 Bsigma <- 2

# ## proposal SDs
 prop.alpha = 0.2 #0.5 #0.2
 prop.beta  = 0.1 #0.5 #0.1
 prop.rho   = 1 #0.5 #0.2
 prop.sigma = 1 #0.3 #0.2
 prop.theta = 0.5
 prop.Bsigma = 1



## initial value for real dataset
## diffusion parameters
# alpha <- 10
# beta  <- 0 
# rho   <- 5.2
# sigma <- 2.7
# Bsigma<- 2

# ## proposal SDs
# prop.alpha = 0.2
# prop.beta  = 0.1 #0.5 #0.1
# prop.rho   = 1
# prop.sigma = 1
# prop.theta = 0.5
# prop.Bsigma = 1



## we fix Kappa, and also if the swiching rate from OU to BM is too small (close to zero), then it would be 
## hard to find a actual switching point. Because for many cases, animal would do OU and the actual switchig 
## probability become nai*lambda_ou_to_bm/kappa, which is very small. Therefore we need set low bound and high bound for lambdas
Kappa = 3.5
sw = 0
osw = -10^4


################## initialise accept number ################################################## 
alaccept  = 0
#beaccept  = 0
sigaccept = 0
rhaccept  = 0
#theaccept = 0
Bsigaccept = 0
staccept = 0

dat <- Sys.Date()
#filetheta3 <- paste("blacBGtheta",dat,".txt", sep = "")
filetheta4 <- paste("blacBGalpha",dat,".txt", sep = "")
#filetheta5 <- paste("blacBGbeta",dat,".txt", sep = "")
filetheta6 <- paste("blacBGsigma",dat,".txt", sep = "")
filetheta7 <- paste("blacBGrho",dat,".txt", sep = "")
filetheta8 <- paste("blacBGBsig",dat,".txt", sep = "")

filetheta9 <- paste("KappaState",dat,".txt", sep = "")
filetheta10 <- paste("KappaLambda",dat,".txt", sep = "")

###initial H for missing or no missing value  for  kalman filter
  Id = diag(nai)
  ha = matrix(c(0),nrow=nai)
  H = matrix(c(ha,ha,Id),nrow=nai)
  oH = H

H_t <- list()                                         #
#
for(i in 1:alen){                                     #
  #
  NAs <- which(is.na(ax[i, ]))                        #   
  #
  if (length(NAs) == 0){H_t[[i]] = H}                 #            
  #
  else{                                               #
    H_t[[i]] = H[-NAs, , drop=FALSE]                  #
  }                                                   #
}                                                     #


lik = 0
olik = -Inf #-10^2


#######################################################
# H_t <- list()                                         #
# #
# for(i in 1:alen){                                     #
#   #
#   N <- which(is.na(raw[i, ]))/2                       #   
#   #
#   if (length(N) == 0){H_t[[i]] = H}                   #            
#   #
#   else{                                               #
#     NAs <- N[seq(1, length(N), 2)]                    #
#     H_t[[i]] = H[-NAs, , drop=FALSE]                  #
#   }                                                   #
# }                                                     #


#######################################################