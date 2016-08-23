## it's weakly identified, but I don't care!!! Stan is slow too (really expensive to compute)

library(rstan)
source("electionfun.R")
source("mcmcfun.R")
source("../psofun.R")

load("datlistsmall.RData")
load("datlistplus.RData")

nbeta <- datlistplus$nbeta
nageedu <- datlistplus$nageedu
npoll <- datlistplus$npoll
nregion <- datlistplus$nregion
nstate <- datlistplus$nstate

stansmall <- stan(file = "election.stan", data = datlistsmall, chains = 1, iter = 2000)
stanplus <- stan(file = "electionplus.stan", data = datlistplus, chains = 1, iter = 2000)

rwgibbstest <- gelmanrwgibbs(1000, rep(0, 80), datlistsmall)

rwgibbstest <- gelmanrwgibbs(10000, rwgibbstest$draws[1000,], datlistsmall,
                             logrwsds = rwgibbstest$logrwsds)

rwgibbstest <- gelmanrwgibbs(20000, rwgibbstest$draws[10000,], datlistsmall,
                             logrwsds = rwgibbstest$logrwsds)

sighat <- cov(rwgibbstest$draws[-c(1:5000), 1:(4 + 51 + 16)])

rwblockgibbstest <- gelmanblockrwgibbs(20000, rep(0, 80), datlistsmall, sighat)

rwblockgibbstest <- gelmanblockrwgibbs(50000, rwblockgibbstest$draws[20000,], datlistsmall,
                                       sighat, logrwsd = rwblockgibbstest$logrwsd,
                                       rwtarget = 0.3)

rwplusgibbstest <- gelmanplusrwgibbs(1000, rep(0, 89), datlistplus)

rwplusgibbstest <- gelmanplusrwgibbs(10000, rwplusgibbstest$draws[1000,], datlistplus,
                                     logrwsds = rwplusgibbstest$logrwsds)

sighat <- cov(rwplusgibbstest$draws[-c(1:5000), 1:(4 + 51 + 16 + npoll)])

rwplusblockgibbstest <- gelmanplusblockrwgibbs(20000, rep(0, 89), datlistplus, sighat)

rwplusblockgibbstest <- gelmanplusblockrwgibbs(50000, rwplusblockgibbstest$draws[20000,],
                                               datlistplus, sighat,
                                               logrwsd = rwplusblockgibbstest$logrwsd)

par(mfrow=c(3,3))
for(i in 1:9){
  plot(ts(rwplusgibbstest$draws[,i]))
}

apply(rwplusgibbstest$acc, 2, mean)
apply(rwplusgibbstest$accs, 2, mean)




xzmat <- datlistsmall$xzmat
alphabetas <- rwblockgibbstest$draws[,1:(4 + 51 + 16)]
mus <- t(tcrossprod(xzmat, alphabetas))

library(MCMCpack)
library(coda)

effsizestest <- effectiveSize(mcmc(mus))

par(mfrow=c(5,4))
for(i in which(effsizestest<10)){
  plot(ts(mus[,i]))
}


par(mfrow=c(3,3))
for(i in 1:9){
  plot(ts(rwgibbstest$draws[,i]))
}

par(mfrow=c(3,3))
for(i in 1:9){
  plot(ts(rwblockgibbstest$draws[,i]))
}



summary(stansmall, pars = c("betay"))$summary
library(MCMCpack)

summary(mcmc(rwgibbstest$draws[-c(1:5000),]))[[1]][1:4,]
summary(mcmc(rwblockgibbstest$draws[-c(1:10000),]))[[1]][1:4,]











stanplus2 <- stan(file = "electionplus2.stan", data = datlistplus, chains = 1, iter = 2000)

rstan::traceplot(stansmall, pars = c("alpha"))

stansmallmod <- stan_model("election.stan")
stanplusmod <- stan_model("electionplus.stan")

stansmalloptim <- optimizing(stansmallmod, datlistsmall)
stanplusoptim <- optimizing(stanplusmod, datlistplus)

musmall <- c(stansmalloptim$par[1:77], log(stansmalloptim$par[78:80]))
muplus <- c(stanplusoptim$par[1:85], log(stanplusoptim$par[86:89]))

lpbestsmall <- gelmanlpost(musmall, datlistsmall)
lpbestplus <- gelmanpluslpost(muplus, datlistplus)

Hsmall <- gelmanlposthess(musmall, datlistsmall)
Hplus <- gelmanpluslposthess(muplus, datlistplus)

nswarm <- 20
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496
nbhd <- sapply(1:nswarm, function(x){return( (x + -3:3 - 1)%%nswarm + 1)}) ## ring-3

plusinit <- matrix(runif((nbeta + nstate + nageedu + npoll + 1 + nregion + 4)*nswarm, -10, 10),
                   ncol = nswarm) + c(muplus)
plusinit[,1] <- muplus

smallinit <- matrix(runif((nbeta + nstate + nageedu + 1 + nregion + 3)*nswarm, -1, 1),
                   ncol = nswarm) + c(musmall)
smallinit[,1] <- musmall

psoiter <- 1000

psoplus <- pso(psoiter, nswarm, 0.9, cognitive, social, plusinit, nbhd, gelmanpluslpost,
               datlistplus, tune = TRUE, rate = 0.5, ccc = 0.01)

psosmall <- pso(psoiter, nswarm, 0.9, cognitive, social, smallinit, nbhd, gelmanlpost,
                datlistsmall, tune = TRUE, rate = 0.5, ccc = 0.01)

psoiter <- 10000
  
psoplus <- pso(psoiter, nswarm, 0.9, cognitive, social, psoplus$pos, nbhd, gelmanpluslpost,
               datlistplus, tune = TRUE, rate = 0.5, ccc = 0.01)

psosmall <- pso(psoiter, nswarm, 0.9, cognitive, social, psosmall$pos, nbhd, gelmanlpost,
                datlistsmall, tune = TRUE, rate = 0.5, ccc = 0.01)


muplus <- psoplus$argmax
lpbestplus <- psoplus$max

musmall <- psosmall$argmax
lpbestsmall <- psosmall$max

plot(ts(psoplus$maxes))

par(mfrow=c(2,1))
plot(ts(psoplus$maxes[-c(1:900)]))
plot(ts(psosmall$maxes[-c(1:900)]))


mcmciter <- 10000
df <- 100

imwgsmall <- gelmanindwithingibbs(mcmciter, musmall, musmall, df, datlistsmall)

H <- gelmanpluslposthess(muplus, datlistplus)

imwgplus <- gelmanplusindwithingibbs(mcmciter, muplus, muplus, df, datlistplus)

mean(imwgplus$acc)

plot(ts(imwgplus$draws[,1:10]))

effectiveSize(mcmc(imwgplus$draws))

plot(ts(imwgplus$draws[,1:10 + 51 + 16 + 8 + 4]))


psoplustest <- pso(psoiter, nswarm, 0.9, cognitive, social, psoplusinit, nbhd, gelmanpluslpost,
                   datlist, tune = TRUE, rate = 0.5, ccc = 0.01)

psoplustest <- pso(psoiter*10, nswarm, 0.9, cognitive, social, psoplustest$pos, nbhd,
                   gelmanpluslpost, datlist, tune = TRUE, rate = 0.5, ccc = 0.01)

save(psoplustest, file = "psoplustest.RData")
load("psoplustest.RData")

mu <- psoplustest$argmax
lpbest <- psoplustest$max

plot(ts(psoplustest$maxes[-c(1:90000)]))


H <- gelmanpluslposthess(mu, datlist)

Sigmatest <- chol2inv(chol(-H))



df <- 100
mcmciter <- 1000

indmetropplustest <- indmetrop(mcmciter, gelmanpluslpost, gelmanpluslposthess,
                               mu, mu, df, lpbest, datlist,
                               tune = FALSE, Sigma = Sigmatest)

mean(indmetropplustest$acc)


imwgplustest <- gelmanplusindwithingibbs(10000, rep(0, 89), mu, df, datlist)

imwgplusstarttest <- gelmanplusindwithingibbs(10000, rep(0, 88), mustart, df, datlist)

c(mean(imwgplustest$acc), mean(imwgplusstarttest$acc))
c(min(effectiveSize(mcmc(imwgplustest$draws))), min(effectiveSize(mcmc(imwgplusstarttest$draws))))

effectiveSize(mcmc(imwgplustest$draws))
effectiveSize(mcmc(imwgplusstarttest$draws))




par(mfrow=c(3,3))
for(i in 4 + 51 + 1:9){
  plot(ts(imwgplusstarttest$draws[,i]))
}






rwgibbstest <- gelmanrwgibbs(1000, rep(0, 80), datlist)

rwgibbstest <- gelmanrwgibbs(10000, rwgibbstest$draws[1000,], datlist,
                             logrwsds = rwgibbstest$logrwsds)

par(mfrow=c(3,3))
for(i in 4 + 51 + 16 + 1:9){
  plot(ts(rwgibbstest$draws[,i]))
}

par(mfrow=c(3,3))
for(i in 1:9){
  plot(ts(rwgibbstest$draws[,i]))
}

library(MCMCpack)

summary(mcmc(rwgibbstest$draws[-c(1:5000),]))[[1]][1:5,]





all.equal(apply(xzmat[,4 + 51 + 1:16], 1, sum), xzmat[,1])

all.equal(apply(xzmat[, 4 + 1:51], 1, sum), xzmat[,1])

all.equal(apply(poll.mat, 1, sum), x.mat[,1])

all.equal(apply(age.edu.mat, 1, sum), x.mat[,1])

all.equal(apply(state.mat, 1, sum), x.mat[,1])
