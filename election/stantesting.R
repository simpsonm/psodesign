library(rstan)
polldat <- read.csv("polldat.csv")
statedat <- read.csv("statedat.csv")
xzmat <- as.matrix(read.csv("xzmat.csv"))
wmat <- as.matrix(read.csv("wmat.csv"))

datlist <- list(polldat = polldat, statedat = statedat, x.mat = x.mat, state.mat = state.mat,
                age.edu.mat = age.edu.mat, region.mat = region.mat,
                age.mat = age.mat, edu.mat = edu.mat,
                betamn = 0, betavar = 1000, sig2a = 1, sig2b = 1)

y <- polldat$bush
female <- polldat$female
black <- polldat$black
prev <- statedat$prev
age <- polldat$age
edu <- polldat$edu
region <- statedat$region
state <- polldat$state

datliststan <- list(y = y, female = female, black = black, state = state, age = age, edu = edu,
                    region = region, prev = prev, nobs = length(y), nstate = max(state),
                    nregion = max(region), nedu = max(edu), nage = max(age),
                    sig2a = 1, sig2b = 1, betamn = 0, betavar = 1000)

stantest <- stan(file = "election.stan", data = datliststan, chains = 1, iter = 2000)

stanpars <- c(summary(stantest, pars = c("beta0", "betaf", "betab", "betafb", "alphastate", "alphaageedu", "betaprev", "alpharegion"))$summary[,1],
              log(summary(stantest, pars = c("sig2state", "sig2ageedu", "sig2region"))$summary[,1]))

stanpars2 <- c(summary(stantest, pars = c("beta0", "betaf", "betab", "betafb", "betaprev", "alphastate", "alphaageedu", "alpharegion"))$summary[,1],
              log(summary(stantest, pars = c("sig2state", "sig2ageedu", "sig2region"))$summary[,1]))


gelmanlpost(stanpars, datlist)
gelmanlpost(mu, datlist)

log_prob(stantest, stanpars2, FALSE)
log_prob(stantest, mu[c(1:4, 4 + 51 + 16 + 1, 5 + 1:(51 + 16 + 5 + 3))], FALSE)

gelmanlpost(stanpars, datlist) - gelmanlpost(mu, datlist)
log_prob(stantest, stanpars2, FALSE) - log_prob(stantest, mu[c(1:4, 4 + 51 + 16 + 1, 5 + 1:(51 + 16 + 5 + 3))], FALSE)


gelmanlpost(stanpars, datlist) - log_prob(stantest, stanpars2, FALSE)
gelmanlpost(mu, datlist) - log_prob(stantest, mu[c(1:4, 4 + 51 + 16 + 1, 5 + 1:(51 + 16 + 5 + 3))], FALSE)



stanex <- extract(stantest)

plot(ts(c(stanex$beta0)))

source("electionfun.R")
source("mcmcfun.R")


rwgibbstest <- gelmanrwgibbs(1000, rep(0, 80), datlist)

rwgibbstest <- gelmanrwgibbs(10000, rwgibbstest$draws[1000,], datlist,
                             logrwsds = rwgibbstest$logrwsds)

par(mfrow=c(3,3))
for(i in 1:9){
  plot(ts(rwgibbstest$draws[,i]))
}

library(MCMCpack)

summary(mcmc(rwgibbstest$draws[-c(1:500),]))[[1]][1:5,]
