// Stan Test
data {
 int n;
 int nMissing;
 vector[n] y;
 vector[n] x1;
 vector[n] x2; // This is the one with missing data
 int whichMissing[nMissing];
 int whichComplete[n - nMissing];
}
parameters {
  real alpha;
  real beta1;
  real beta2;
  real<lower=0> sigma;
  real x2Mean;
  real x2Cor;
  real<lower=0> x2SD;
  vector[nMissing] x2Missing;
}
model{
  alpha ~ normal(0,4);
  beta1 ~ normal(0,4);
  beta2 ~ normal(0,4);
  sigma ~ cauchy(0, 1);
  x2Mean ~ normal(0,1);
  x2Cor ~ normal(0, 1);
  x2SD ~ cauchy(0, 1);
  // Assume no relationship between x1 and x2
  x2[whichComplete] ~ normal(x2Mean + x2Cor*x1[whichComplete], x2SD);
  x2Missing ~ normal(x2Mean+ x2Cor*x1[whichMissing], x2SD);
  
  y[whichComplete] ~ normal(alpha + beta1*x1[whichComplete] + beta2*x2[whichComplete], sigma);
  y[whichMissing] ~ normal(alpha + beta1*x1[whichMissing] + beta2*x2Missing, sigma);
}