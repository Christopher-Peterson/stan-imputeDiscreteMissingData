####################################
### Simulate Data, hierarchically

library(shinystan)
library(rstan)
library(mvtnorm)
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


n = 300
nClust = 3
nPerClust = n/nClust
k = 3  # not including intercept
nMissing = n * k * 0.15
clustID = rep(1:nClust, each = nPerClust)

set.seed(2)
if(nMissing>0){
  propCols = round(gtools::rdirichlet(1, rep(3,k)) *nMissing)
  missingPosN = unlist(lapply(1:k, function(i) sample.int(n, size = propCols[i])))
  nMissing <- length(missingPosN)
  orderN = order(missingPosN)
  missingPosN = sort(missingPosN)
  missingPosK = rep(1:k, times = propCols)+1
  missingPosK = missingPosK[orderN]
  missingRows = unique(missingPosN)
  nMissingRows = length(missingRows)
  missingPerRow = as.vector(table(missingPosN))
  wholeRows = setdiff(1:n, missingRows)
  whichClustMissing <- unique(clustID[missingRows])
  nClustMissing <- length(whichClustMissing)
  missingIDTranslate = rep(0, nClust); missingIDTranslate[whichClustMissing] <- 1:nClustMissing
  missingClustID <-   missingIDTranslate[clustID[missingPosN]] 
}
set.seed(3)
betaHier = c(3, 8, -4, -7)
L = rcor(1,k+1, 3, chol=TRUE)
betaSD = c(2, .3, 1, 1.8)
Omega = tcrossprod(diag(betaSD) %*% t(L))
beta = rmvnorm(nClust,betaHier, sigma=Omega)
sigma = 1

# Make x vary by level.
# This will require x to be generated as a normal distribution then "stepped"
baseXMeans = rnorm(k)
clustXMeans = rmvnorm(nClust, baseXMeans)
x = cbind(1, do.call(rbind, lapply(1:n, function(i) rnorm(k, clustXMeans[clustID[i],]))) > 0)
y = sapply(1:n, function(i) rnorm(1, x[i,] %*% beta[clustID[i],], sigma))
k = k+1 # update to include intercept
LKJParam = 2


xTmp = x
if(nMissing > 0){ for(i in 1:nMissing)
  xTmp[missingPosN[i],missingPosK[i]]<- 100000
  nZero = sum(xTmp==0);  # number of zeros
  nOne = sum(xTmp[,-1]==1);  # number of ones.
  onePosN = unlist(lapply(2:k, function(i) which(xTmp[,i]==1)))   # row position of missing variable; these should be sorted from lowest to highest.
  onePosK = unlist(lapply(2:k, function(i) rep(i, sum(xTmp[,i]==1))))    # column position of missing variable, corresponding to missingPosN
  zeroPosN = unlist(lapply(2:k, function(i) which(xTmp[,i]==0)))   # row position of missing variable; these should be sorted from lowest to highest.
  zeroPosK = unlist(lapply(2:k, function(i) rep(i, sum(xTmp[,i]==0))))   # column position of missing variable, corresponding to missingPosN
  rm(xTmp)
}

if(nMissing == 0) {
  fit = stan("hierNoMissingTest.stan", chains = 1, iter = 10)
  fit = stan(fit = fit,seed=2, chains = 4, cores=4, iter = 100, control=list(adapt_delta=.9, max_treedepth = 13))
}
setwd("/home/peterson/Documents/Files/R/discreteMissingData")
fit = stan("discreteMissingDataTestHier.stan", chains = 1, iter = 10)
fit = stan(fit = fit,seed=2, chains = 5, cores=5, iter = 1000, control=list(adapt_delta=.9, max_treedepth = 13))
print(fit,"beta")




tmpY = y[clustID==1]
tmpx1 = x[clustID==1,2]
tmpx2 = x[clustID==1,3]
tmpx3 = x[clustID==1,4]
lm(tmpY~tmpx1+tmpx2+tmpx3)
x11()
boxplot(tmpY, tmpx1)




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
 
set.seed(3)
n = 500
k = 3  # not including intercept
nMissing = n * k * .30
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
