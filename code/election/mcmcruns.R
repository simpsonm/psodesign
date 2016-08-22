source("../psofun.R")
source("electionfun.R")
source("mcmcfun.R")
load("datlistsmall.RData")
load("datlistplus.RData")
library(rstan)

niter <- 20000
nswarm <- 100
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496
## nbhd <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)}) ## ring-1
## nbhd <- sapply(1:nswarm, function(x){return( (x + -3:3 - 1)%%nswarm + 1)}) ## ring-3
## nbhd <- matrix(1:nswarm, ncol=nswarm, nrow=nswarm) ## global
models <- c("small", "poll")

stansmall <- stan(file = "election.stan", data = datlistsmall, chains = 1, iter = 2000)
stanplus <- stan(file = "electionplus.stan", data = datlistplus, chains = 1, iter = 2000)

imwgsmall <- gelmanindwithingibbs(mcmciter, musmall, musmall, df, datlistsmall)

imwgplus <- gelmanplusindwithingibbs(mcmciter, muplus, muplus, df, datlistplus)


H <- gelmanpluslposthess(mu, datlist)

Sigmatest <- chol2inv(chol(-H))

indmetropplustest <- indmetrop(mcmciter, gelmanpluslpost, gelmanpluslposthess,
                               mu, mu, df, lpbest, datlist,
                               tune = FALSE, Sigma = Sigmatest)



rwgibbstest <- gelmanrwgibbs(1000, rep(0, 80), datlist)

rwgibbstest <- gelmanrwgibbs(10000, rwgibbstest$draws[1000,], datlist,
                             logrwsds = rwgibbstest$logrwsds)
