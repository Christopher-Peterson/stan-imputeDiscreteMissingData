library(shinystan)
library(rstan)
n = 500
k = 3  # not including intercept
nMissing = n * k * .05
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

setwd( "/media/peterson/B965-8339/R/discreteMissingData")

fit = stan("discreteMissingDataTest.stan", chains = 1, iter = 10)
fit = stan(fit = fit, chains = 5, cores=5, iter = 1000, control=list(adapt_delta=.9999, max_treedepth = 15))
launch_shinystan(fit)

#fit2 = stan("discreteMissingDataTest.stan", chains = 1, iter = 10)
#fit2 = stan(fit = fit2, chains = 5, cores=5, iter = 1000)
#launch_shinystan(fit2)



### Make a test here to evaluate dmd_probs
if(FALSE) {
testProbs = c(.1, .8, .3, .7)
probCombs = function(probs){
  qs = 1-probs
  lp = log(probs)
  lq = log(qs)
  k = length(probs)
  combinations= do.call(c, lapply(1:(k-1), function(j) combn(k, j, simplify = FALSE)))
  c(sum(lq), sapply(1:length(combinations), function(i)
      sum(c(lp[combinations[[i]]], lq[-combinations[[i]]]))), sum(lp))
}


compProbs = probCombs(testProbs)
probs = testProbs
k = length(testProbs)

dmdpTest = stan("dmd_probs.stan", chains = 1, iter = 1, algorithm = "Fixed_param")





bitMask = function(i, d, revTwos = 2^((d-1):0)){
    iReal <- i;
    out = numeric(d)
#    browser()
  for(j in 1:d){
      outStep <- ((iReal-revTwos[j]) > 0)*1 ;
      iReal <- iReal - revTwos[j] * outStep;
      out[j] <- outStep;
    }
    out;
  }

bitMask(3,4)

lapply(1:16, function(i) bitMask(i, 4))





}