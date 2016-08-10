library(rstan)
source("popfun.R")
source("mcmcfun.R")
source("psofun.R")

load("popdat/popdat.RData")
ndelta <- 5
z <- popdat[,1]
xmat <- as.matrix(popdat[,2])
smat <- popdat[,-c(1:2)]
aphi <- 1
bphi <- 1
asig <- 1
bsig <- 1
d <- ndelta + 1
Einv <- diag(ndelta)
mubeta <- array(rep(0, ncol(xmat)), dim = ncol(xmat))
sigbeta <- 10

fullstandat <- list(nobs = length(z), ndelta = ndelta, nell = ndelta*(ndelta + 1)/2,
                    nbeta = ncol(xmat), zobs = z,
                    lzobs = log(z), S = smat[,1:ndelta], X = xmat,
                    a_phi = aphi, b_phi = bphi, a_sig = asig, b_sig = bsig,
                    mu_beta = mubeta, sig_beta = sigbeta,
                    d_omega = d, Einv_omega = Einv)

stanfullmod <- stan_model("poppoisfull.stan")

stanfullfit <- stan(file = "poppoisfull.stan", data = fullstandat, chains = 1, iter = 1)

system.time(stanfulloptim <- optimizing(stanfullmod, fullstandat))

nswarm <- 100
psoiter <- 10000
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496


nell <- ndelta*(ndelta + 1)/2
nbeta <- 1
poispsofullinit <- matrix(runif((ndelta + nbeta + nell)*nswarm, -1, 1), ncol = nswarm)

ominvscale <- diag(ndelta)
K <- commutation.matrix(ndelta)
M <- elimination.matrix(ndelta)
N2 <- 2*N.matrix(ndelta)
R <- matrix(0, nell, nell)
MKprime <- tcrossprod(M, K)
diag(R)[(ndelta+1)*(1:ndelta) - (1:ndelta)*(1 + 1:ndelta)/2 + 1:ndelta - ndelta] <-
  (ndelta + 2 - 1:ndelta)
dvelldvell <- - M%*%tcrossprod(kronecker(diag(ndelta), ominvscale), M)
datlistfull <- list(dat = popdat, nbeta = nbeta, ndelta = ndelta,
                    betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                    aphi = 1, bphi = 1, omdf = ndelta + 1, ominvscale = diag(ndelta),
                    K = K, M = M, N2 = N2, R = R, MKprime = MKprime, dvelldvell = dvelldvell)
nbhd <- sapply(1:nswarm, function(x){return( (x + -3:3 - 1)%%nswarm + 1)})

psoiter <- 10000
psopoisfull <- pso(psoiter, nswarm, inertia, cognitive, social, poispsofullinit, nbhd,
                   poislpostfull, datlistfull)

psomu <- psopoisfull$argmax
psolpbest <- psopoisfull$max

stanmuprime <- stanfulloptim$par[c(22, 1 + 1:(ndelta + nell))]
standiag <- stanmuprime[nbeta + ndelta + nell - ndelta + 1:ndelta]
stanoff <- stanmuprime[nbeta + ndelta + 1:(nell - ndelta)]
stanL <- matrix(stanfulloptim$par[grep("L", names(stanfulloptim$par))], ndelta, ndelta)
stanell <- stanL[lower.tri(stanL, TRUE)]
stanell2 <- stanfulloptim$par[grep("ell", names(stanfulloptim$par))]
stanmu <- c(stanmuprime[c(1, 1 + 1:ndelta)], stanell2)

L1 <- matrix(stanfulloptim$par[grep("L", names(stanfulloptim$par))], ndelta, ndelta)
L2 <- matrix(0, ndelta, ndelta)
L2[lower.tri(L2, TRUE)] <- stanell

stanlpbest <- poislpostfull(stanmu, datlistfull)

log_prob(stanfullfit,
         c(stanfulloptim$par[1:(nbeta + ndelta + nell - ndelta)],
           log(stanfulloptim$par[1:ndelta + nbeta + ndelta + nell - ndelta])),
         FALSE)

rbfgs <- optim(rep(1, length(psomu)), poislpostfull, datlist = datlistfull, method = "BFGS",
               control = list(fnscale = -1))

rbfgs <- optim(rep(1, nbeta + ndelta + nell), poislpostfull, poisgradfull,
               datlist = datlistfull, method = "BFGS",
               control = list(fnscale = -1, maxit = 10000, reltol = .Machine$double.eps))

mu <- rbfgs$par
lpbest <- rbfgs$value

poislpostfull(stanmu, datlistfull)

metroptest <- indmetrop(10000, poislpostfull, poislposthessfull, stanmu, stanmu, 100, stanlpbest,
                        datlistfull, tune = FALSE)

metroptest2 <- indmetrop(10000, poislpostfull, poislposthessfull, psomu, psomu, 100,
                        psolpbest, datlistfull, tune = FALSE)

mean(metroptest$acc)

mean(metroptest2$acc)


poisgradfull(rbfgs$par, datlistfull)

poisgradfull(psomu, datlistfull)

H <- poislposthessfull(rbfgs$par, datlistfull)

Sigma <- chol2inv(chol(-H))




stanfulloptim$value






poisfulltest <- stan(file = "poppoisfull.stan", data = fullstandat, chains = 1, iter = 1)
poisfullfit <- stan(fit = poisfulltest, data = fullstandat, chains = 1, iter = 20000,
                    control = list(adapt_delta = 0.9, max_treedepth = 11))

poisiidtest <- stan(file = "poppoisiid.stan", data = fullstandat, chains = 1, iter = 1000,
                    warmup = 1000,
                    control = list(adapt_delta = 0.9, max_treedepth = 11))
poisiidfit <- stan(fit = poisiidtest, data = fullstandat, chains = 1, iter = 1000,
                   warmup = 0,
                   control = list(adapt_delta = 0.9, max_treedepth = 11))


lognormfulltest <- stan(file = "poplognormfull.stan", data = fullstandat, chains = 1, iter = 1)
lognormfullfit <- stan(fit = lognormfulltest, data = fullstandat, chains = 1, iter = 1000,
                       control = list(adapt_delta = 0.9, max_treedepth = 11))

lognormfullfit <- stan(fit = lognormfullfit, data = fullstandat, chains = 1, iter = 10000,
                       warmup = 2000,
                       control = list(adapt_delta = 0.9, max_treedepth = 11))

lognormiidtest <- stan(file = "poplognormiid.stan", data = fullstandat, chains = 1, iter = 1)
lognormiidfit <- stan(fit = lognormiidtest, data = fullstandat, chains = 1, iter = 10000)

### poisson will likely need higher adapt_delta & max_treedepth:
### default: control = list(adapt_delta = 0.8, max_treedepth = 10)
### try:     control = list(adapt_delta = 0.9, max_treedepth = 11)








library(matrixcalc)
nbeta <- 1
nell <- ndelta*(ndelta + 1)/2
ominvscale <- diag(ndelta)
K <- commutation.matrix(ndelta)
M <- elimination.matrix(ndelta)
N2 <- 2*N.matrix(ndelta)
R <- matrix(0, nell, nell)
MKprime <- tcrossprod(M, K)
diag(R)[(ndelta+1)*(1:ndelta) - (1:ndelta)*(1 + 1:ndelta)/2 + 1:ndelta - ndelta] <-
  (ndelta + 2 - 1:ndelta)
dvelldvell <- - M%*%tcrossprod(kronecker(diag(ndelta), ominvscale), M)
datlistfull <- list(dat = popdat, nbeta = nbeta, ndelta = ndelta,
                    betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                    aphi = 1, bphi = 1, omdf = ndelta + 1, ominvscale = diag(ndelta),
                    K = K, M = M, N2 = N2, R = R, MKprime = MKprime, dvelldvell = dvelldvell)
datlistiid <- datlistfull

source("mcmcfun.R")

mcmciter <- 20000
poisfullrwgibbsout <- poisfullrwgibbs(mcmciter, rep(0, nbeta + ndelta + nell), datlistfull)
sighatfull <- cov(poisfullrwgibbsout$draws[-c(1:10000),1:(nbeta + ndelta)])
cholhatfull <- t(chol(sighatfull))
blockinit <- poisfullrwgibbsout$draws[mcmciter + 1,]
poisfullblockrwgibbsout <- poisfullblockrwgibbs(50000, blockinit, datlistfull, cholhatfull)

mcmciter <- 20000
lnormfullrwgibbsout <- lnormfullrwgibbs(mcmciter, rep(0, nbeta + ndelta + nell+1), datlistfull)
sighatfull <- cov(lnormfullrwgibbsout$draws[-c(1:10000),1:(nbeta + ndelta)])
cholhatfull <- t(chol(sighatfull))
blockinit <- lnormfullrwgibbsout$draws[mcmciter + 1,]
lnormfullblockrwgibbsout <- lnormfullblockrwgibbs(50000, blockinit, datlistfull, cholhatfull)


mcmciter <- 20000
poisiidrwgibbsout <- poisiidrwgibbs(mcmciter, rep(0, nbeta + ndelta + 1), datlistiid)
sighatiid <- cov(poisiidrwgibbsout$draws[-c(1:10000),1:(nbeta + ndelta)])
cholhatiid <- t(chol(sighatiid))
blockinit <- poisiidrwgibbsout$draws[mcmciter + 1,]
poisiidblockrwgibbsout <- poisiidblockrwgibbs(50000, blockinit, datlistiid, cholhatiid)

mcmciter <- 20000
lnormiidrwgibbsout <- lnormiidrwgibbs(mcmciter, rep(0, nbeta + ndelta + 1 + 1), datlistiid)
sighatiid <- cov(lnormiidrwgibbsout$draws[-c(1:10000),1:(nbeta + ndelta)])
cholhatiid <- t(chol(sighatiid))
blockinit <- lnormiidrwgibbsout$draws[mcmciter + 1,]
lnormiidblockrwgibbsout <- lnormiidblockrwgibbs(50000, blockinit, datlistiid, cholhatiid)



cbind(summary(mcmc(poisfullblockrwgibbsout$draws[-c(1:10000),]))[[1]][,1],
summary(poisfullfit, pars = c("beta", "delta", "ell"))$summary[,1])

cbind(summary(mcmc(lnormfullblockrwgibbsout$draws[-c(1:10000),]))[[1]][,1],
summary(lognormfullfit, pars = c("beta", "delta", "ell", "phi2"))$summary[,1])

cbind(summary(mcmc(poisiidblockrwgibbsout$draws[-c(1:10000),]))[[1]][,1],
summary(poisiidfit, pars = c("beta", "delta", "sig2"))$summary[,1])

cbind(summary(mcmc(lnormiidblockrwgibbsout$draws[-c(1:10000),]))[[1]][,1],
summary(lognormiidfit, pars = c("beta", "delta", "sig2", "phi2"))$summary[,1])


