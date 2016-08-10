data {
  int<lower=1> nobs;
  int<lower=1> ndelta;
  int<lower=1> nbeta;
  int<lower=1> nell;
  vector[nobs] lzobs;
  matrix[nobs, ndelta] S;
  matrix[nobs, nbeta] X;
  real<lower=ndelta - 1> d_omega;
  cov_matrix[ndelta] Einv_omega;
  vector[nbeta] mu_beta;
  real<lower=0> sig_beta;
  real<lower=0> a_phi;
  real<lower=0> b_phi;
}
transformed data {
  vector[ndelta] zero_delta;
  vector[nobs] Xmu_beta;
  matrix[ndelta, ndelta] Einv_chol;
  Einv_chol <- cholesky_decompose(Einv_omega);
  for(n in 1:ndelta){
    zero_delta[n] <- 0.0;
  }
  Xmu_beta <- X*mu_beta;
}
parameters {
  vector[nbeta] beta_raw;
  vector[ndelta] delta;
  vector[nell-ndelta] zwishoff;
  vector<lower=0>[ndelta] zwishdiag;
  real<lower=0> phi2;
}
transformed parameters{
  vector[nobs] muz;
  vector[nbeta] beta;
  real<lower=0> phi;
  matrix[ndelta, ndelta] Zmat;
  matrix[ndelta, ndelta] L;
  beta <- beta_raw + mu_beta;
  for(i in 1:ndelta){
    for(j in 1:ndelta){
      if(i == j){
	Zmat[j,i] <- sqrt(zwishdiag[i]);
      } else if(j > i){
	## compiler will complain about integer division here. it's fine.
	Zmat[j,i] <- zwishoff[j + (i-1)*ndelta - i*(i-1)/2 - i];  
      } else {
	Zmat[j,i] <- 0.0;
      }
    }
  }
  L <- Einv_chol*Zmat;
  muz <- X*beta_raw + Xmu_beta + S*delta;
  phi <- sqrt(phi2);
  beta <- beta_raw + mu_beta;
}
model {
  for(n in 1:nobs)
    lzobs[n] ~ normal(muz[n], phi);
    delta ~ multi_normal_prec(zero_delta, tcrossprod(L));
  beta_raw ~ normal(0.0, sig_beta);
  for(i in 1:ndelta){
    for(j in i:ndelta){
      if(i == j){
	zwishdiag[i] ~ chi_square(d_omega - i + 1);
      } else {
	## compiler will complain about integer division here. it's fine.
	zwishoff[j + (i-1)*ndelta - i*(i-1)/2 - i] ~ normal(0.0, 1);
      }
    }
  }
  phi2 ~ inv_gamma(a_phi, b_phi);
}
generated quantities {
  vector[ndelta*(ndelta + 1)/2] ell;
  for(i in 1:ndelta){
    for(j in i:ndelta){
      ell[j + (i-1)*ndelta - i*(i-1)/2] <- L[j,i];
    }
  }
}

