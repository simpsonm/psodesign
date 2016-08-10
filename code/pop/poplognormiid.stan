data {
  int<lower=1> nobs;
  int<lower=1> ndelta;
  int<lower=1> nbeta;
  vector[nobs] lzobs;
  matrix[nobs, ndelta] S;
  matrix[nobs, nbeta] X;
  real<lower=0> a_sig;
  real<lower=0> b_sig;
  vector[nbeta] mu_beta;
  real<lower=0> sig_beta;
  real<lower=0> a_phi;
  real<lower=0> b_phi;
}
transformed data {
  vector[ndelta] zero_delta;
  vector[nobs] Xmu_beta;
  for(n in 1:ndelta){
    zero_delta[n] <- 0.0;
  }
  Xmu_beta <- X*mu_beta;
}
parameters {
  vector[nbeta] beta_raw;
  vector[ndelta] delta;
  real<lower=0> sig2;
  real<lower=0> phi2;
}
transformed parameters{
  vector[nobs] muz;
  vector[nbeta] beta;
  real<lower=0> phi;
  real<lower=0> sig;
  muz <- X*beta_raw + Xmu_beta + S*delta;
  beta <- beta_raw + mu_beta;
  phi <- sqrt(phi2);
  sig <- sqrt(sig2);
  
}
model {
  for(n in 1:nobs)
    lzobs[n] ~ normal(muz[n], phi);
  delta ~ normal(0.0, sig);
  beta_raw ~ normal(0.0, sig_beta);
  sig2 ~ inv_gamma(a_sig, b_sig);
  phi2 ~ inv_gamma(a_phi, b_phi);
}
