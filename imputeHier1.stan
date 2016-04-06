/* Random Intercepts multiple imputation
 * This has a couple of assumptions:
 * There are correlations among the group-level x means for each covariate AND among individual-level x's.
 */
data {
  int n;
  int nClust;
  int<lower = 1, upper = nClust> clustID[n];
  int d;
  int nMissing;
  vector[n] y;
  row_vector[d] x[n];
  int whichMissingD[nMissing]; // Column position of missing index
  int whichMissingN[nMissing]; // Row position of missing index
 real LKJParam;
}
parameters {
  vector[nClust] alphaRaw;
  real<lower=0> alphaSD;
  real alphaMean;
  
  vector[d] beta;
  real<lower=0> sigma;
  
  matrix[d,nClust] xMeanRaw; // The mean offset of each cluster's x-values from the grand mean
  row_vector[d] xHyperMean; // The grand mean of the x-values
  vector<lower=0>[d] xMeanSD; // The standard deviation for each cluster's x-means around the grand mean
#  vector<lower=0>[d] xSDRaw[nClust]; // The standard deviation of a the x-values around the cluster mean, per-cluster
#  vector<lower=0>[d] xHyperSD; // scaling parameter for the xSD;
  vector<lower=0>[d] xSD; // scaling parameter for the xSD;
  cholesky_factor_corr[d] xCor; // Correlation between different covariates
  cholesky_factor_corr[d] xMeanCor; // How is this different from the above?  
  
  vector[nMissing] xMissing;
}
transformed parameters{
  matrix[nClust,d] xMean;
  vector[nClust] alpha;
#  vector<lower=0>[d] xSD[nClust]
  alpha <- alphaRaw * alphaSD + alphaMean;
  xMean <- rep_vector(1.0,nClust) * xHyperMean + (diag_pre_multiply(xMeanSD, xCor) * xMeanRaw)'; 
#  for(i in 1:nClust)
#    xSD[i] <- xSDRaw[i] .* xHyperSD;
}
model{
  row_vector[d] xFull[n];
  row_vector[d] xMeanFull[n];
  vector[n] eta;
  row_vector[d] xEta[n];
  xFull <- x;
  for(i in 1:nMissing)
    xFull[whichMissingN[i],whichMissingD[i]] <- xMissing[i];
  alphaRaw ~ normal(0,1);
  alphaMean ~ normal(0,4);
  alphaSD ~ cauchy(0,1);
  
  beta ~ normal(0,4);
  sigma ~ cauchy(0, 1);
  to_vector(xMeanRaw) ~ normal(0,1);
  xHyperMean ~ normal(0,1);
  xMeanCor ~ lkj_corr_cholesky(LKJParam);
  xCor ~ lkj_corr_cholesky(LKJParam);
#  for(i in 1:nClust)
#    xSDRaw[i] ~ cauchy(0, 1);
  xMeanSD ~ cauchy(0, 1);
  xSD ~ cauchy(0, 1);
  for(i in 1:n){
    eta[i] <- alpha[clustID[i]] + xFull[i] * beta;
    xMeanFull[i] <- xMean[clustID[i]];
  }
  xFull ~ multi_normal_cholesky(xMeanFull, diag_pre_multiply(xSD, xCor));
  y ~ normal(eta, sigma); 
}
generated quantities{
  matrix[d,d] Sigma;
  matrix[d,d] SigmaMean;
  Sigma <- multiply_lower_tri_self_transpose(xCor);
  SigmaMean <- multiply_lower_tri_self_transpose(xMeanCor);
  
}