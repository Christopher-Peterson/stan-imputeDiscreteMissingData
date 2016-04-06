## Hier Impute .R
library(mvtnorm)
library(rstan)
library(magrittr)

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

 
set.seed(2);
n = 1000;
d = 5;
nMissing = round(n * d * 1/10);
nClust = 10;

whichMissingD = sort(sample(1:d, nMissing, replace = TRUE));
whichMissingN = unlist(lapply(1:d, function(i) sample(1:n, sum(whichMissingD==i))));
while(sum(as.vector(table(whichMissingN))>=(d-1))>0) # This makes sure there's no rows where all but 1 covariate are missing
  whichMissingN = unlist(lapply(1:d, function(i) sample(1:n, sum(whichMissingD==i))));

alphaMean = .3;
alphaSD = .1;
alpha = rnorm(nClust, alphaMean, alphaSD);
beta = c(-.3, .4, .1, -.7, .5)
sigma = 1
clustID = rep(1:nClust, each = n/nClust)
corMatMean = rcor(1,d, 4)
corMatShared = rcor(1,d, 2, chol = TRUE)
corMatIndiv = rcor(nClust, d, 100, chol = TRUE)
corMat = lapply(1:nClust, function(i)
  cov2cor(crossprod(corMatShared + corMatIndiv[[i]])))
xMean = rmvnorm(nClust, sigma = corMatMean)
xFull = do.call(rbind, lapply(1:nClust, function(i) rmvnorm(n/nClust, mean = xMean[i,], sigma = corMat[[i]])))

y = rnorm(n,alpha[clustID] + xFull %*% beta)
x = xFull
missingReal = numeric(nMissing)

for(i in 1:nMissing){
  x[whichMissingN[i], whichMissingD[i]] = -500
  missingReal[i] = xFull[whichMissingN[i], whichMissingD[i]]
}
LKJParamMean = 4
LKJParamShared = 2
LKJEtaScale = .4


stanImputeHier = stan("imputeHier2.stan", chains = 2, iter = 600, cores = 1)
 r eadSum(stanCorGen, alpha, beta, sigma)
  readSum(stanCor, x2Missing)
  readSum(stanCor, x2Mean, x2SD)

  