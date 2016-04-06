// Objective: generalize this to a matrix with many possible missing data points
data {
  int n;
  int d;
  int nMissing;
  vector[n] y;
  row_vector[d] x[n];
  int whichMissingD[nMissing]; // Column position of missing index
  int whichMissingN[nMissing]; // Row position of missing index
# int whichCompleteD[n*d - nMissing]; // Column
# int whichCompleteN[n*d - nMissing]; // Row
 real LKJParam;
}
parameters {
  real alpha;
  vector[d] beta;
  real<lower=0> sigma;
  vector[d] xMean;
  vector<lower=0>[d] xSD;
  cholesky_factor_corr[d] xCor;
  vector[nMissing] xMissing;
}
model{
  row_vector[d] xFull[n];
  vector[n] eta;
  xFull <- x;
  for(i in 1:nMissing)
    xFull[whichMissingN[i],whichMissingD[i]] <- xMissing[i];
  alpha ~ normal(0,4);
  beta ~ normal(0,4);
  sigma ~ cauchy(0, 1);
  xMean ~ normal(0,1);
  xCor ~ lkj_corr_cholesky(LKJParam);
  xSD ~ cauchy(0, 1);
  xFull ~ multi_normal_cholesky(xMean, diag_pre_multiply(xSD, xCor));
  for(i in 1:n)
    eta[i] <- alpha + xFull[i] * beta;
  y ~ normal(eta, sigma); 
}
generated quantities{
  matrix[d,d] Sigma;
  Sigma <- multiply_lower_tri_self_transpose(xCor);
}