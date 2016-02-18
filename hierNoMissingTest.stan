functions{
  /// Functions to marginalize out missing categorical covariates for imputation
  /// Written by Christopher R. Peterson
  int dmd_pow2(int k){ // returns 2^k
    int out;
    out <- 1;
    for(i in 1:k)
      out <- out * 2;
    return out;
  }
  // convert i into a binary bit-mask of d length; revTwos is the powers of 2 from 1 to d-1, in reverse order.
  row_vector dmd_bit_mask(int i, int d, int[] revTwos){
    row_vector[d] out;
    int iReal;
    int outStep;
    iReal <- i;
    for(j in 1:d){
      outStep <- int_step(iReal-revTwos[j]);
      iReal <- iReal - revTwos[j] * outStep;
      out[j] <- outStep;
    }
    return out;
  }

  // take a row with d missing data points and return a 
  // vector of all possible combinations of linear predictors
  vector dmd_etas(row_vector x, vector beta, int[] missingPos){
    int d;
    int k;
    int outSize;
    vector[dmd_pow2(size(missingPos))] out;
    d <- size(missingPos);
    k <- cols(x);
    outSize <- dmd_pow2(d);
    { 
      row_vector[d] posState; //<lower=0, upper=1>
      int twos[d]; 
      row_vector[k] tmpX;
      //Initialize
        tmpX <- x;
        for(i in 1:d)// Twos need to be in reverse order.
          twos[i] <- dmd_pow2(d-i);
      // create each output
      for(i in 1:outSize) {
        for(j in 1:d) // increment posState; this is conceptually similar to using bitwise operators for binary masking
        posState <- dmd_bit_mask(i, d, twos); 
        tmpX[missingPos] <- posState;
        out[i] <- tmpX * beta;
      }
    }
    return out;
}
   
  // Take the p parameter for each missing value in a row and 
  // return the joint distribution, expressed as a vector.  
  // These should match with the linear predictors generated by dmd_etas()
  // Return values are on log scale
  vector dmd_probs(row_vector p){
    int d;
    int outSize;
    vector[dmd_pow2(cols(p))] out;
    d <- cols(p);
      outSize <- dmd_pow2(d);
    { row_vector[d] posState; //<lower=0, upper=1>
      row_vector[d] posSign; //  -1 if posState is 1, 1 if posState is 0; //<lower=0, upper=1>
      int twos[d]; 
      row_vector[d] tmpP;
      // Twos need to be in reverse order.
        for(i in 1:d)
          twos[i] <- dmd_pow2(d-i);
      // create each output
      for(i in 1:outSize) {
        posState <- dmd_bit_mask(i, d, twos); 
        posSign <- -2*posState + 1;
        tmpP <- posState + posSign .* p;  // equivalent to p or 1-p's, combined.
        out[i] <- sum(log(tmpP));
      }
    }
    return out;
  }  //returns normal_log for a single row with some missing data
  real dmd_normal_single_log(real y, row_vector x, real sigma, vector beta, row_vector p,  int[] missingPos){
    int d;
    real out;
    int vecSize;
    d <- size(missingPos);
    vecSize <- dmd_pow2(d);
    {
      vector[vecSize] lprobs;
      vector[vecSize] etas;
      vector[vecSize] lpVec;
      lprobs <- dmd_probs(p);
      etas <- dmd_etas(x, beta, missingPos);
      for(i in 1:vecSize)
        lpVec[i] <- normal_log(y, etas[i], sigma) + lprobs[i];
      out <- log_sum_exp(lpVec);
    }
    return out;
}
  
  // This function assumes everything is fully hierarchical, and all lengths are n.
  real dmd_normal_log(vector y, matrix x, vector sigma, vector[] beta, row_vector p, 
    int[] missingRows, /*Which rows have missing data; must be sorted*/
    int[] missingPerRow, /*How many on that row are missing*/ 
    int[] wholeRows, /*which rows have all data*/
    int[] missingPos ){
    int n; // rows
    int k; // columns
    int d; // missing values
    real lp_out;
    int nRowsMissing;
    int nRowsWhole;
    
    n <- rows(y);
    k <- cols(x);
    d <- size(missingPos);
    nRowsMissing <- size(missingRows);
    nRowsWhole <- n - nRowsMissing;
    {// In this block: take the logLikelihood of all rows that are NOT missing
      vector[nRowsWhole] eta;
      
      for(i in 1:nRowsWhole)
        eta[i] <- x[wholeRows[i]] * beta[wholeRows[i]];
      
      lp_out <- normal_log(y[wholeRows], eta, sigma[wholeRows]);
    }
    {// In this block: use a for loop to calculate the lp of each other row, then add it to lp_out
      int start;
      int stop;
      start <- 1;
      stop <- 0;
      for(i in 1:nRowsMissing){
        stop <- stop + missingPerRow[i];
        lp_out <- lp_out + dmd_normal_single_log(y[missingRows[i]], x[missingRows[i]], 
          sigma[missingRows[i]], beta[missingRows[i]], p[start:stop], missingPos[start:stop]);
        start <- start + missingPerRow[i];
        
      }
    }
    return lp_out;
}
  row_vector Phi_rvec(row_vector x){
    row_vector[cols(x)] out;
    for(i in 1:cols(x))
      out[i] <- Phi_approx(x[i]);
    return out;
  }
  row_vector inv_logit_rvec(row_vector x){
    row_vector[cols(x)] out;
    for(i in 1:cols(x))
      out[i] <- inv_logit(x[i]);
    return out;
  }
}
data{
  int n;
  int<lower=2> k; // num Covariates, including intercept
  int nClust; // Number of clusters/hierarchical groups
  int clustID[n]; // cluster membership for each row;
  matrix[n, k] x; // column 1 should be all 1's.; this is dropped for the fully latent variable version
  vector[n] y;
  real LKJParam; 
  /*
  int nZero;  // number of zero's
  int nOne;  // number of ones.
  int onePosN[nOne]; // row position of missing variable, sorted from lowest to highest.
  int onePosK[nOne]; // column position of missing variable, corresponding to missingPosN
  int zeroPosN[nZero]; // row position of missing variable, sorted from lowest to highest.
  int zeroPosK[nZero]; // column position of missing variable, corresponding to missingPosN
*/}
/*transformed data{
  matrix[n,k] x;
  for(i in 1:n)
    x[i,1]<- 1;
  for(i in 1:nOne)
    x[onePosN[i], onePosK[i]] <- 1;
  for(i in 1:nZero)
    x[zeroPosN[i], zeroPosK[i]] <- 0;
}*/
parameters{
  row_vector[k] betaHier; 
  vector<lower=0>[k] betaHierSD; 
  matrix[k,nClust] betaRaw;
  real<lower=0> sigma;
  cholesky_factor_corr[k] L;
}

transformed parameters{
  matrix[nClust,k] beta;
  beta <- rep_vector(1, nClust) * betaHier + 
  (diag_pre_multiply(betaHierSD, L) * betaRaw)';
}
model{
  vector[n] eta;
  betaHier ~ normal(0, 3);
  betaHierSD ~ cauchy(0, 2);
  to_vector(betaRaw) ~ normal(0,1);
  sigma ~ cauchy(0,2.5);
  L ~ lkj_corr_cholesky(LKJParam);
  for(i in 1:n)
    eta[i] <- dot_product(x[i], beta[clustID[i]]);
  y ~ normal(eta, sigma);
}