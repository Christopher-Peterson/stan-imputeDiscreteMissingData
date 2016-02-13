n = 500
k = 3  # not including intercept
nMissing = n * k * .15
set.seed(3)
propCols = round(gtools::rdirichlet(1, rep(3,k)) *nMissing)
missingPosN = unlist(lapply(1:k, function(i) sample.int(500, size = propCols[i])))
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


fit = stan("~/Documents/Files/RProjects/discreteMissingDataTest.stan", chains = 1, iter = 10)
fit = stan(fit = fit, chains = 2, iter = 100)









