## Multiple Imputation Test
#### Objective: examine multiple imputation for missing covariates in Stan
library(magrittr)
library(rstan);
library(mvtnorm)
library(shinystan)
source("~/R/quotr.R")
readSum = function(fit, ..., cols = -c(2,3,9,10))
  summary(fit, quotr(...))$summary[,cols] %>% round(3)


## First: simulate some data; start with a non-hierarchical regression model with two covariates, one of which is partially Missing.
#####################################
## No correlation between x's.
set.seed(2);
n = 2000;
x1 = rnorm(n);
x2Complete = rnorm(n);
nMissing = round(n/4);
nComplete = n - nMissing;
whichMissing = sample(1:n, nMissing);
whichComplete = setdiff(1:n, whichMissing)
x2 = x2Complete;
x2[whichMissing] = -500;
xFull = cbind(x1,x2Complete);

# coefs:

alpha = .3;
beta = c(-.3, .4)
sigma = 1
y = rnorm(n,alpha + xFull %*% beta)


## write stan model, then set stan data

stanTest1 = stan("impute1.stan", chains = 4, iter = 500)
  readSum(stanTest1, alpha, beta1, beta2, sigma)
  readSum(stanTest1, x2Missing)
  readSum(stanTest1, x2Mean, x2SD)

# does correlation improve it any?
stanCor = stan("impute2.stan", chains = 4, iter = 500)
  readSum(stanCor, alpha, beta1, beta2, sigma)
  readSum(stanCor, x2Missing)
  readSum(stanCor, x2Mean, x2SD, x2Cor)

#######################
## New assumption: correlation between x1 and x2.
set.seed(2);
n = 2000;
nMissing = round(n/4);
nComplete = n - nMissing;
whichMissing = sample(1:n, nMissing);
whichComplete = setdiff(1:n, whichMissing)

alpha = .3;
beta = c(-.3, .4)
sigma = 1
rho = -.6 # correlation between x1 and x2
corMat = matrix(c(1, rho, rho, 1), 2,2)

xFull = rmvnorm(n, sigma = corMat)
x1 = xFull[,1]
x2 = xFull[,2]
x2[whichMissing] = -500;
y = rnorm(n,alpha + xFull %*% beta)

# No correlation
stanNC = stan("impute1.stan", chains = 4, iter = 500)
  readSum(stanNC, alpha, beta1, beta2, sigma)
  readSum(stanNC, x2Missing)
  readSum(stanNC, x2Mean, x2SD)

# with correlation
stanCor = stan("impute2.stan", chains = 4, iter = 500)
  readSum(stanCor, alpha, beta1, beta2, sigma)
  readSum(stanCor, x2Missing)
  readSum(stanCor, x2Mean, x2SD)

whichMissingD = rep(2,nMissing)
whichMissingN = whichMissing
d = 2
x = xFull
x[whichMissingN, whichMissingD] = -500
LKJParam = 2

stanCorGen = stan("impute3.stan", chains = 1, iter =200)
  readSum(stanCorGen, alpha, beta, sigma)
  readSum(stanCor, x2Missing)
  readSum(stanCor, x2Mean, x2SD)


##################################
## Expansion: multiple predictors, all of which can have missing data


rcor <- function(n, d, eta = 1, chol = FALSE, list1=FALSE){
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
  if(n==1 && list1==FALSE) out = out[[1]]
  out
}   # random correlation; eta is LKJ parameter

 
set.seed(2);
n = 2000;
d = 5;
nMissing = round(n * d * 1/10);

whichMissingD = sort(sample(1:d, nMissing, replace = TRUE));
whichMissingN = unlist(lapply(1:d, function(i) sample(1:n, sum(whichMissingD==i))));
while(sum(as.vector(table(whichMissingN))>=(d-1))>0) # This makes sure there's no rows where all but 1 covariate are missing
  whichMissingN = unlist(lapply(1:d, function(i) sample(1:n, sum(whichMissingD==i))));

alpha = .3;
beta = c(-.3, .4, .1, -.7, .5)
sigma = 1
corMat = rcor(1,d, 3)

xFull = rmvnorm(n, sigma = corMat)
y = rnorm(n,alpha + xFull %*% beta)
x = xFull
x[whichMissingN, whichMissingD] = -500
LKJParam = 2


stanCorGen = stan("impute3.stan", chains = 4, iter = 600, cores = 4)
  readSum(stanCorGen, alpha, beta, sigma)
  readSum(stanCor, x2Missing)
  readSum(stanCor, x2Mean, x2SD)
