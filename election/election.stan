data {
  //comment
  int<lower=1> nobs;
  int<lower=1> nstate;
  int<lower=1> nregion;
  int<lower=1> nedu;
  int<lower=1> nage;
  int<lower=1> nbeta;
  int<lower=0, upper=1> y[nobs];
  matrix[nobs, nbeta] xmat;
  matrix[nobs, nstate] statemat;
  matrix[nobs, nage*nedu] ageedumat;
  matrix[nstate, nregion] regionmat;
  vector[nstate] prev;
  real<lower=0> sig2a;
  real<lower=0> sig2b;
  real betamn;
  real<lower=0> betavar;
}
transformed data {
  real<lower=0> betasd;
  betasd <- sqrt(betavar);
}
parameters {
  vector[nbeta] betay;
  vector[nstate] alphastate;
  vector[nage*nedu] alphaageedu;
  vector[1] betaprev;
  vector[nregion] alpharegion;
  real<lower=0> sig2state;
  real<lower=0> sig2ageedu;
  real<lower=0> sig2region;
}
transformed parameters{
  vector[nobs] mu;
  real<lower=0, upper=1> p[nobs];
  real<lower=0> sigstate;
  real<lower=0> sigregion;
  real<lower=0> sigageedu;
  sigstate <- sqrt(sig2state);
  sigregion <- sqrt(sig2region);
  sigageedu <- sqrt(sig2ageedu);
  mu <- xmat*betay + statemat*alphastate + ageedumat*alphaageedu;
  for(n in 1:nobs){
    p[n] <- exp(mu[n])/(1 + exp(mu[n]));
  }
}
model {
  y ~ bernoulli(p);
  for(n in 1:nstate){
    alphastate[n] ~ normal(regionmat[n,] * alpharegion + prev[n] * betaprev, sigstate);
  }
  alpharegion ~ normal(0.0, sigregion);
  alphaageedu ~ normal(0.0, sigageedu);
  sig2ageedu ~ inv_gamma(sig2a, sig2b);
  sig2state ~ inv_gamma(sig2a, sig2b);
  sig2region ~ inv_gamma(sig2a, sig2b);
  betay ~ normal(betamn, betasd);
  betaprev ~ normal(betamn, betasd);
}
