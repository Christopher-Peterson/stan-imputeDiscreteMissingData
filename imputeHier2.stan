/* Random Intercepts multiple imputation
 * This has a couple of assumptions:
 * There are correlations among the group-level x means for each covariate AND among individual-level x's.
 */
functions{
  matrix cov2Cor(matrix L_Sigma){//L1 and L2 should be cholesky factors with the same dimensionality
    matrix[cols(L_Sigma), cols(L_Sigma)] out;
    int d;
    d <- cols(L_Sigma);
    {
      matrix[d, d] Sigma;
      matrix[d, d] DiagMat;
      Sigma <- multiply_lower_tri_self_transpose(L_Sigma);
      DiagMat <- diag_matrix(diagonal(Sigma));
      for(i in 1:d) DiagMat[i,i] <- sqrt(1/DiagMat[i,i]);
    
      out <- DiagMat * Sigma * DiagMat;
    }
    return cholesky_decompose(out);
  }
}
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
  vector[nMissing] missingReal; // Actual value of missing data (for checking)
 real LKJParamMean;
 real LKJParamShared;
 real LKJEtaScale;
 
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
  cholesky_factor_corr[d] xSharedCor; // Correlation between different covariates
  cholesky_factor_corr[d] xMeanCor; // How is this different from the above?  
  cholesky_factor_corr[d] xClustCor[nClust]; // Correlation between different covariates
  real<lower=0> etaLKJ;
  vector[nMissing] xMissing;
}
transformed parameters{
  matrix[nClust,d] xMean;
  vector[nClust] alpha;
  matrix[d,d] xCor[nClust];

  alpha <- alphaRaw * alphaSD + alphaMean;
  xMean <- rep_vector(1.0,nClust) * xHyperMean + (diag_pre_multiply(xMeanSD, xMeanCor) * xMeanRaw)'; 
  
  for(i in 1:nClust)
    xCor[i] <- cov2Cor(xClustCor[i] + xSharedCor);
}
model{
  row_vector[d] xFull[n];
  vector[n] eta;
  xFull <- x;
  for(i in 1:nMissing)
    xFull[whichMissingN[i],whichMissingD[i]] <- xMissing[i];
  alphaRaw ~ normal(0,1);
  alphaMean ~ normal(0,4);
  alphaSD ~ cauchy(0,1);
  etaLKJ ~ normal(0, LKJEtaScale);
  beta ~ normal(0,4);
  sigma ~ cauchy(0, 1);
  to_vector(xMeanRaw) ~ normal(0,1);
  xHyperMean ~ normal(0,1);
  xMeanCor ~ lkj_corr_cholesky(LKJParamMean);
  xSharedCor ~ lkj_corr_cholesky(LKJParamShared);
  for(i in 1:nClust)
    xSharedCor ~ lkj_corr_cholesky(1/etaLKJ);
  xMeanSD ~ cauchy(0, 1);
  xSD ~ cauchy(0, 1);
  for(i in 1:n){
    eta[i] <- alpha[clustID[i]] + xFull[i] * beta;
  }
  for(i in 1:n)
    xFull[i] ~ multi_normal_cholesky(xMean[clustID[i]], diag_pre_multiply(xSD, xCor[clustID[i]]));
  y ~ normal(eta, sigma); 
}
generated quantities{
  matrix[d,d] Sigma[nClust];
  matrix[d,d] SigmaMean;
  vector[nMissing] missDiff;
  real xRMSE;
  
  for(i in 1:nClust)
    Sigma[i] <- multiply_lower_tri_self_transpose(xCor[i]);
  SigmaMean <- multiply_lower_tri_self_transpose(xMeanCor);
  missDiff <- missingReal - xMissing;
  xRMSE <- sqrt(missDiff' * missDiff);
}