####################################
### Simulate Data

library(shinystan)
library(rstan)
n = 500
k = 3  # not including intercept
nMissing = n * k * .15
set.seed(3)
propCols = round(gtools::rdirichlet(1, rep(3,k)) *nMissing)
missingPosN = unlist(lapply(1:k, function(i) sample.int(500, size = propCols[i])))
nMissing <- length(missingPosN)
orderN = order(missingPosN)
missingPosN = sort(missingPosN)
missingPosK = rep(1:k, times = propCols)+1
missingPosK = missingPosK[orderN]
missingRows = unique(missingPosN)
nMissingRows = length(missingRows)
missingPerRow = as.vector(table(missingPosN))
wholeRows = setdiff(1:n, missingRows)

beta = c(3, 8, -4, -7)
sigma = 1
x = cbind(1, matrix(rbinom(n*k, 1, .5),n,k))
y = rnorm(n, x %*% beta, sigma) 
k = k+1 # update to include intercept
library(rstan)
LKJParam = 2

xTmp = x
for(i in 1:nMissing)
  xTmp[missingPosN[i],missingPosK[i]]<- 100000

nZero = sum(xTmp==0);  # number of zeros
nOne = sum(xTmp[,-1]==1);  # number of ones.
onePosN = unlist(lapply(2:k, function(i) which(xTmp[,i]==1)))   # row position of missing variable; these should be sorted from lowest to highest.
onePosK = unlist(lapply(2:k, function(i) rep(i, sum(xTmp[,i]==1))))    # column position of missing variable, corresponding to missingPosN
zeroPosN = unlist(lapply(2:k, function(i) which(xTmp[,i]==0)))   # row position of missing variable; these should be sorted from lowest to highest.
zeroPosK = unlist(lapply(2:k, function(i) rep(i, sum(xTmp[,i]==0))))   # column position of missing variable, corresponding to missingPosN

#setwd( "/media/peterson/B965-8339/R/discreteMissingData")

fit = stan("discreteMissingDataTest.stan", chains = 1, iter = 10)
fit = stan(fit = fit, chains = 5, cores=5, iter = 1000, control=list(adapt_delta=.9, max_treedepth = 15))

print(fit, "beta")
 launch_shinystan(fit)
 
library(mvtnorm)
library(magrittr)
 
####################################################
####### Now, do the same as above but with correlation
 
rcor <- function(n, d, eta = 1, chol = FALSE, forceList=FALSE){
  k = choose(d,2)
  out = lapply(1:n, function(i) {
    z = matrix(0, d,d)
    z[upper.tri(z)] = rbeta(k, eta, eta) * 2 - 1  ## Constrain to (-1,1)
    # Many have to convert to lower.tri if cholesky factors are to be done properly
    
    ## Calculate determinant 
    z = z + diag(d)
    w = matrix(0, d,d)
    #Convert z into cholesky factor
    w[1,] = z[1,]
    w[1,1] = 1
    for(i in 2:d){ ## See section 56.11 of Stan manual v2.9 (chokesly Correlation Transforms) for this.
      for(j in i:d)
        w[i,j] = z[i,j] * sqrt(1-sum(w[1:(i-1),j]^2))
    }
    if(chol==FALSE) w = t(w) %*% w
    w
  })
  if(n==1 && forceList==FALSE) out = out[[1]]
  out
}   # random correlation; eta is LKJ parameter

set.seed(3)
n = 500
k = 3  # not including intercept
nMissing = n * k * .40
propCols = round(gtools::rdirichlet(1, rep(3,k)) *nMissing)
missingPosN = unlist(lapply(1:k, function(i) sample.int(500, size = propCols[i])))
nMissing <- length(missingPosN)
orderN = order(missingPosN)
missingPosN = sort(missingPosN)
missingPosK = rep(1:k, times = propCols)+1
missingPosK = missingPosK[orderN]
missingRows = unique(missingPosN)
nMissingRows = length(missingRows)
missingPerRow = as.vector(table(missingPosN))
wholeRows = setdiff(1:n, missingRows)
LKJParam = 2
beta = c(3, 8, -4, -7)
sigma = 1

Omega = rcor(1,k,1)
x = cbind(1, rmvnorm(n, sigma = Omega) > 0)
y = rnorm(n, x %*% beta, sigma) 
k = k+1 # update to include intercept
xTmp = x
for(i in 1:nMissing)
  xTmp[missingPosN[i],missingPosK[i]]<- 100000
nZero = sum(xTmp==0);  # number of zeros
nOne = sum(xTmp[,-1]==1);  # number of ones.
onePosN = unlist(lapply(2:k, function(i) which(xTmp[,i]==1)))   # row position of missing variable; these should be sorted from lowest to highest.
onePosK = unlist(lapply(2:k, function(i) rep(i, sum(xTmp[,i]==1))))    # column position of missing variable, corresponding to missingPosN
zeroPosN = unlist(lapply(2:k, function(i) which(xTmp[,i]==0)))   # row position of missing variable; these should be sorted from lowest to highest.
zeroPosK = unlist(lapply(2:k, function(i) rep(i, sum(xTmp[,i]==0))))   # column position of missing variable, corresponding to missingPosN
rm(xTmp)


fit = stan("discreteMissingDataTest.stan", chains = 1, iter = 10)
fit = stan(fit = fit, chains = 5, cores=5, iter = 1000, control=list(adapt_delta=.9, max_treedepth = 15))

print(fit, "beta")
print(fit2, "beta")
