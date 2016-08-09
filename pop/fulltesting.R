source("popfun.R")
source("mcmcfun.R")
source("psofun.R")
load("popdat/popdat.RData")
library(rstan)
library(MCMCpack)



nbeta <- 1
nswarms <- c(100, 1000)
ndeltas <- c(8)
psoiter <- 1000
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496
mcmciter <- 50000
df <- 100

out <- NULL
psoarglist <- list()
bbpso1arglist <- list()
bbpso2arglist <- list()

for(nswarm in nswarms){
  i <- which(nswarm == nswarms)
  psoarglist[[i]] <- list()
  bbpso1arglist[[i]] <- list()
  bbpso2arglist[[i]] <- list()
  for(ndelta in ndeltas){
    j <- which(ndelta == ndeltas)
    print(c(i,j))
    nbhd <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)})
    nell <- ndelta*(ndelta + 1)/2
    poispsofullinit <- matrix(runif((ndelta + nbeta + nell)*nswarm, -1, 1), ncol = nswarm)
    lnormpsofullinit <- matrix(runif((ndelta + nbeta + nell + 1)*nswarm, -1, 1), ncol = nswarm)
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
    psopoisfull <- pso(psoiter, nswarm, inertia, cognitive, social, poispsofullinit, nbhd,
                      poislpostfull, datlistfull)
    poisfullmu <- psopoisfull$argmax
    poisfullinit <- psopoisfull$argmax
    poisfulllpbest <- psopoisfull$max
    psopoisfull2 <- bbpso(psoiter, nswarm, 1, 0.1, poispsofullinit, nbhd,
                   poislpostfull, 1, TRUE, c(0.5, 0.5), datlistfull)
    poisfullmu2 <- psopoisfull2$argmax
    poisfullinit2 <- psopoisfull2$argmax
    poisfulllpbest2 <- psopoisfull2$max
    psopoisfull5 <- bbpso(psoiter, nswarm, 1, 0.3, poispsofullinit, nbhd,
                         poislpostfull, 5, TRUE, c(0.5, 0.5), datlistfull)
    poisfullmu5 <- psopoisfull5$argmax
    poisfullinit5 <- psopoisfull5$argmax
    poisfulllpbest5 <- psopoisfull5$max
    poisfullind <- indmetrop(mcmciter, poislpostfull, poislposthessfull,
                         poisfullmu, poisfullmu, df,
                         poisfulllpbest, datlistfull)
    poisfullind2 <- indmetrop(mcmciter, poislpostfull, poislposthessfull,
                             poisfullmu2, poisfullmu2, df,
                             poisfulllpbest2, datlistfull)
    poisfullind5 <- indmetrop(mcmciter, poislpostfull, poislposthessfull,
                             poisfullmu5, poisfullmu5, df,
                             poisfulllpbest5, datlistfull)
    out <- rbind(out, data.frame(ndelta = ndelta, pso = "pso", nswarm = nswarm,
                                 max = poisfulllpbest,
                                 acc = mean(poisfullind$acc[-c(1:10000)])))
    out <- rbind(out, data.frame(ndelta = ndelta, pso = "bbpso1", nswarm = nswarm,
                                 max = poisfulllpbest2,
                                 acc = mean(poisfullind2$acc[-c(1:10000)])))
    out <- rbind(out, data.frame(ndelta = ndelta, pso = "bbpso2", nswarm = nswarm,
                                 max = poisfulllpbest5,
                                 acc = mean(poisfullind5$acc[-c(1:10000)])))
    psoarglist[[i]][[j]] <- poisfullmu
    bbpso1arglist[[i]][[j]] <- poisfullmu2
    bbpso2arglist[[i]][[j]] <- poisfullmu5
  }
}


out



source("mcmcfun.R")
load("popdat/popdat.RData")
library(MCMCpack)

z <- popdat[,1]
xmat <- as.matrix(popdat[,2])
smat <- popdat[,-c(1:2)]
nbeta <- 1
ndelta <- 5



source("popfun.R")
source("psofun.R")




nell <- ndelta*(ndelta + 1)/2
poispsofullinit <- matrix(runif((ndelta + nbeta + nell)*nswarm, -1, 1), ncol = nswarm)
lnormpsofullinit <- matrix(runif((ndelta + nbeta + nell + 1)*nswarm, -1, 1), ncol = nswarm)
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

mcmciter <- 50000

lnormfullrwgibbsout <- lnormfullrwgibbs(mcmciter,
                                        rep(0, nbeta + ndelta + ndelta*(ndelta + 1)/2 + 1),
                                        datlistfull)

sighatfull <- cov(lnormfullrwgibbsout$draws[-c(1:10000),1:(nbeta + ndelta)])
cholhatfull <- t(chol(sighatfull))
blockinit <- lnormfullrwgibbsout$draws[mcmciter + 1,]
lnormfullblockrwgibbsout <- lnormfullblockrwgibbs(mcmciter, blockinit, datlistfull, cholhatfull,
                                                  rwtarget = c(0.4, 0.5))

poisfullrwgibbsout <- poisfullrwgibbs(mcmciter,
                                        rep(0, nbeta + ndelta + ndelta*(ndelta + 1)/2),
                                        datlistfull)

sighatfull <- cov(poisfullrwgibbsout$draws[-c(1:10000),1:(nbeta + ndelta)])
cholhatfull <- t(chol(sighatfull))
blockinit <- poisfullrwgibbsout$draws[mcmciter + 1,]
poisfullblockrwgibbsout <- poisfullblockrwgibbs(mcmciter, blockinit, datlistfull, cholhatfull,
                                                rwtarget = c(0.4, 0.5))




c(min(effectiveSize(mcmc(poisfullrwgibbsout$draws[-c(1:20000),]))), min(effectiveSize(mcmc(poisfullblockrwgibbsout$draws[-c(1:20000),]))))

c(min(effectiveSize(mcmc(lnormfullrwgibbsout$draws[-c(1:20000),]))), min(effectiveSize(mcmc(lnormfullblockrwgibbsout$draws[-c(1:20000),]))))




niter <- 1000
nswarm <- 100
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496
nbhd <- matrix(1:nswarm, ncol=nswarm, nrow=nswarm)
lnormfullinit <- matrix(runif((ndelta + nbeta + ndelta*(ndelta + 1)/2 + 1)*nswarm, -1, 1), ncol = nswarm)

psolnormfull <- pso(niter, nswarm, inertia, cognitive, social, lnormfullinit, nbhd,
                    lnormlpostfull, )

lfullmu <- psolnormfull$argmax
lfullinit <- psolnormfull$argmax
lfullpbest <- psolnormfull$max

df <- 8

dat <- popdat
dat[,1] <- exp(rnorm(nrow(dat), xmat*10 + smat[,1:2]%*%c(-2,13)))

psolnormfull <- pso(niter, nswarm, inertia, cognitive, social, lnormfullinit, nbhd,
                    lnormlpostfull, dat, nbeta, ndelta)

lfullmu <- psolnormfull$argmax
lfullinit <- psolnormfull$argmax
lfullpbest <- psolnormfull$max


testH <- lnormlposthessfull(lfullmu, dat, nbeta, ndelta)

lnormfullindtest <- indmetrop(mcmciter, lnormlpostfull, lnormlposthessfull, lfullmu, lfullmu, df,
                             lfullpbest, dat, nbeta, ndelta)

mean(lnormfullindtest$acc)




dat <- popdat
dat[,1] <- rpois(nrow(dat), exp(xmat*4 + smat[,1:2]%*%c(-5,8)))


poisfullinit <- matrix(runif((ndelta + nbeta + ndelta*(ndelta + 1)/2)*nswarm, -1, 1), ncol = nswarm)
psopoisfull <- pso(niter, nswarm, inertia, cognitive, social, poisfullinit, nbhd,
                   poislpostfull, dat, nbeta, ndelta)


poisfullmu <- psopoisfull$argmax
poisfullinit <- psopoisfull$argmax
poisfullpbest <- psopoisfull$max

df <- 2



L <- matrix(0, ndelta, ndelta)
L[lower.tri(L, TRUE)] <- poisfullmu[-c(1:(nbeta + ndelta))]
omega <- tcrossprod(L)
omega

poisfullindtest <- indmetrop(250*mcmciter, poislpostfull, poislposthessfull,
                             poisfullmu, poisfullmu, df,
                             poisfullpbest, dat, nbeta, ndelta)


mean(poisfullindtest$acc)


poisfullindtest$draws2 <- poisfullindtest$draws

elldraws <- poisfullindtest$draws2[,-c(1:(nbeta + ndelta))]
L <- matrix(0, ndelta, ndelta)
for(i in 1:nrow(elldraws)){
  L[lower.tri(L, TRUE)] <- elldraws[i,]
  omega <- tcrossprod(L)
  elldraws[i,] <- omega[lower.tri(omega, TRUE)]
}
poisfullindtest$draws2[,-c(1:(nbeta + ndelta))] <- elldraws

fullinit <- poisfullmu

poisfullrwgibbstest <- poisfullrwgibbs(4*mcmciter, fullinit, dat, nbeta, ndelta)
sighatfull <- cov(poisfullrwgibbstest$draws[-c(1:10000),1:(nbeta + ndelta)])
cholhatfull <- t(chol(sighatfull))
blockinitfull <- poisfullrwgibbstest$draws[mcmciter + 1,]
poisfullblockrwgibbstest <- poisfullblockrwgibbs(4*mcmciter, fullinit, dat,
                                                nbeta, ndelta, cholhatfull)



poisfullindtest2 <- poisfullindmetrop(mcmciter*250, poisfullinit, poisfullmu, poisfullpbest,
                                      df, dat, nbeta, ndelta)

poisfullindtest2$draws2 <- poisfullindtest2$draws
elldraws <- poisfullindtest2$draws2[,-c(1:(nbeta + ndelta))]
L <- matrix(0, ndelta, ndelta)
for(i in 1:nrow(elldraws)){
  L[lower.tri(L, TRUE)] <- elldraws[i,]
  omega <- tcrossprod(L)
  elldraws[i,] <- omega[lower.tri(omega, TRUE)]
}
poisfullindtest2$draws2[,-c(1:(nbeta + ndelta))] <- elldraws


mean(poisfullindtest$acc)
mean(poisfullindtest2$acc)


rbind(apply(poisfullrwgibbstest$draws[-c(1:10000),], 2, mean),
      apply(poisfullblockrwgibbstest$draws[-c(1:10000),], 2, mean),
      apply(poisfullindtest$draws2[-c(1:10000),], 2, mean),
      apply(poisfullindtest2$draws2[-c(1:10000),], 2, mean))

c(effectiveSize(poisfullrwgibbstest$draws[,4]),
  effectiveSize(poisfullblockrwgibbstest$draws[,4]),
  effectiveSize(poisfullindtest$draws2[,4]),
  effectiveSize(poisfullindtest2$draws2[,4]))

summary(mcmc(poisfullrwgibbstest$draws[-c(1:10000),]))[[2]]
summary(mcmc(poisfullblockrwgibbstest$draws[-c(1:10000),]))[[2]]
summary(mcmc(poisfullindtest$draws2[-c(1:10000),]))[[2]]
summary(mcmc(poisfullindtest2$draws2[-c(1:10000),]))[[2]]



L <- matrix(0, ndelta, ndelta)
L[lower.tri(L, TRUE)] <- poisfullmu[-c(1:(nbeta + ndelta))]
omega <- tcrossprod(L)
omega



rbind(summary(poisfullfit)$summary[c(8,2,3,4,6,7),1],
      apply(poisfullrwgibbstest$draws[-c(1:10000),], 2, mean),
      apply(poisfullblockrwgibbstest$draws[-c(1:10000),], 2, mean),
      apply(poisfullindtest$draws2[-c(1:10000),], 2, mean),
      apply(poisfullindtest2$draws2[-c(1:10000),], 2, mean))


summary(mcmc(poisfullrwgibbstest$draws[-c(1:10000),]))[[2]]
summary(mcmc(poisfullblockrwgibbstest$draws[-c(1:10000),]))[[2]]
summary(mcmc(poisfullindtest$draws2[-c(1:10000),]))[[2]]
summary(mcmc(poisfullindtest2$draws2[-c(1:10000),]))[[2]]





library(rstan)
aphi <- 1
bphi <- 1
asig <- 1
bsig <- 1
d <- ndelta + 1
Einv <- diag(ndelta)
mubeta <- array(rep(0, ncol(xmat)), dim = ncol(xmat))
sigbeta <- 10
fullstandat <- list(nobs = length(z), ndelta = ndelta, nbeta = ncol(xmat), zobs = dat[,1],
                    lzobs = log(z), S = smat[,1:ndelta], X = xmat,
                    a_phi = aphi, b_phi = bphi, a_sig = asig, b_sig = bsig,
                    mu_beta = mubeta, sig_beta = sigbeta,
                    d_omega = d, Einv_omega = Einv)


poisfulltest <- stan(file = "poppoisfull.stan", data = fullstandat, chains = 1, iter = 1)
poisfullfit <-  stan(fit = poisfulltest,        data = fullstandat, chains = 1, iter = 50000)



lnormlpostfull <- function(pars, dat, nbeta, ndelta, betamn = 0, betavar = 100,
                           omdf = ndelta + 1, ominvscale = diag(ndelta),
                           aphi = 1, bphi = 1){
  z <- dat[,1]
  lz <- z ##log(z)
  xmat <- dat[,1 + 1:nbeta]
  beta <- pars[1:nbeta]
  nell <- ndelta*(ndelta + 1)/2
  lphi2 <- pars[nbeta + ndelta + nell + 1]
  phi2 <- exp(lphi2)
  if(ndelta > 0){
    delta <- pars[nbeta + 1:ndelta]
    smat <- dat[, 1 + nbeta + 1:ndelta]
    ell <- pars[nbeta + ndelta + 1:nell]
    L <- matrix(0, ndelta, ndelta)
    L[lower.tri(L, diag = TRUE)] <- ell
    elldiag <- diag(L)
    sigma <- crossprod(L)
    omega <- chol2inv(t(L))
    y <- as.matrix(xmat)%*%beta + as.matrix(smat)%*%delta
  } else {
    y <- as.matrix(xmat)%*%beta
  }
  datpart <- sum(dnorm(lz, y, sqrt(phi2), log = TRUE))
  betapart <- -crossprod(beta - betamn)/(2*betavar)
  phipart <- -lphi2*aphi - bphi/phi2
  if(ndelta > 0){
    ellpart <- - sum((omdf + 1 + 1:ndelta)/2*log(elldiag^2)) -
      (crossprod(delta, omega)%*%delta + sum(diag(ominvscale%*%omega)))/2
    out <- datpart + betapart + phipart + ellpart
  } else {
    out <- datpart + betapart + phipart
  }
  if(is.null(out) | is.na(out)){
    return(-Inf)
  }
  return(out)
}

lnormlposthessfull <- function(pars, dat, nbeta, ndelta, betamn = 0, betavar = 100,
                               omdf = ndelta + 1, ominvscale = diag(ndelta),
                               aphi = 1, bphi = 1){
  npars <- length(pars)
  ndat <- nrow(dat)
  z <- dat[, 1]
  lz <- z ##log(z)
  xmat <- dat[, 1 + 1:nbeta]
  beta <- pars[1:nbeta]
  nell <- ndelta*(ndelta+1)/2
  lphi2 <- pars[nbeta + ndelta + nell + 1]
  phi2 <- exp(lphi2)
  if(ndelta > 0){
    smat <- dat[, 1 + nbeta + 1:ndelta]
    delta <- pars[nbeta + 1:ndelta]
    ell <- pars[nbeta + ndelta + 1:nell]
    L <- matrix(0, ndelta, ndelta)
    L[lower.tri(L, diag = TRUE)] <- ell
    elldiag <- diag(L)
    omega <- crossprod(L)
    y <- as.matrix(xmat)%*%beta + as.matrix(smat)%*%delta
  } else {
    y <- as.matrix(xmat)%*%beta
  }
  xx <- crossprod(xmat)
  dbetadbeta <- - xx/phi2 - diag(1/betavar, nbeta)
  dlphi2dlphi2 <- - (bphi + crossprod(lz - y)/2)
  if(nbeta > 1){
    dlphi2dbeta <- - apply(drop(lz - y)*xmat, 2, sum)/phi2
  } else {
    dlphi2dbeta <- - sum(drop(lz - y)*xmat)/phi2
  }
  out <- matrix(0, npars, npars)
  out[1:nbeta, 1:nbeta] <- dbetadbeta
  out[nbeta + ndelta + nell + 1, nbeta + ndelta + nell + 1] <- dlphi2dlphi2
  out[nbeta + ndelta + nell + 1, 1:nbeta] <- dlphi2dbeta
  out[1:nbeta, nbeta + ndelta + nell + 1] <- dlphi2dbeta
  if(ndelta > 0){
    dlphi2ddelta <- - apply(drop(lz - y)*smat, 2, sum)/phi2
    ss <- crossprod(smat)
    xs <- crossprod(xmat, smat)
    ddeltaddelta <- - ss/phi2 - omega
    dbetaddelta <- - xs/phi2
    K <- commutation.matrix(ndelta)
    M <- elimination.matrix(ndelta)
    N2 <- 2*N.matrix(ndelta)
    Indelta <- diag(ndelta)
    ddeltadvell <- -kronecker(t(delta), Indelta)%*%N2%*%tcrossprod(kronecker(L, Indelta), M)
    R <- matrix(0, nell, nell)
    diag(R)[(ndelta+1)*(1:ndelta) - (1:ndelta)*(1 + 1:ndelta)/2 + 1:ndelta - ndelta] <-
      (ndelta + 1 + 1:ndelta)
    MKprime <- tcrossprod(M, K)
    dvelldvell <- - (MKprime%*%tcrossprod(kronecker(tcrossprod(delta), Indelta), MKprime) +
                     diag(nell)) +
      - (M%*%tcrossprod(kronecker(Indelta, ominvscale), M)) - R%*%diag(1/ell^2)
    out[nbeta + 1:ndelta, nbeta + 1:ndelta] <- ddeltaddelta
    out[1:nbeta, nbeta + 1:ndelta] <- dbetaddelta
    out[nbeta + 1:ndelta, 1:nbeta] <- t(dbetaddelta)
    out[nbeta + ndelta + 1:nell, nbeta + ndelta + 1:nell] <- dvelldvell
    out[nbeta + 1:ndelta, nbeta + ndelta + 1:nell] <- ddeltadvell
    out[nbeta + ndelta + 1:nell, nbeta + 1:ndelta] <- t(ddeltadvell)
    out[nbeta + ndelta + nell + 1, 1:ndelta + nbeta] <- dlphi2ddelta
    out[1:ndelta + nbeta, nbeta + ndelta + nell + 1] <- dlphi2ddelta
  }
  return(out)
}



poislpostfull <- function(pars, dat, nbeta, ndelta, betamn = 0, betavar = 100,
                          omdf = ndelta + 1, ominvscale = diag(ndelta)){
  z <- dat[,1]
  xmat <- dat[,1 + 1:nbeta]
  beta <- pars[1:nbeta]
  if(ndelta > 0){
    nell <- ndelta*(ndelta + 1)/2
    delta <- pars[nbeta + 1:ndelta]
    smat <- dat[, 1 + nbeta + 1:ndelta]
    ell <- pars[nbeta + ndelta + 1:nell]
    L <- matrix(0, ndelta, ndelta)
    L[lower.tri(L, diag = TRUE)] <- ell
    elldiag <- diag(L)
    omega <- crossprod(L)
    y <- as.matrix(xmat)%*%beta + as.matrix(smat)%*%delta
  } else {
    y <- as.matrix(xmat)%*%beta
  }
  lambda <- exp(y)
  lampart <- sum(z*y - lambda)
  betapart <- -crossprod(beta - betamn)/(2*betavar)
  if(ndelta > 0){
    ellpart <- sum((omdf + 1 - 1:ndelta)/2*log(elldiag^2)) -
      (crossprod(crossprod(L, delta)) + sum(diag(ominvscale%*%omega)))/2
    out <- lampart + ellpart + betapart
  } else {
    out <- lampart + betapart
  }
  if(is.null(out) | is.na(out)){
    return(-Inf)
  }
  return(out)
}

poislposthessfull <- function(pars, dat, nbeta, ndelta, betamn = 0, betavar = 100,
                              omdf = ndelta + 1, ominvscale = diag(ndelta)){
  npars <- length(pars)
  ndat <- nrow(dat)
  z <- dat[, 1]
  xmat <- dat[, 1 + 1:nbeta]
  beta <- pars[1:nbeta]
  if(ndelta > 0){
    nell <- ndelta*(ndelta+1)/2
    smat <- dat[, 1 + nbeta + 1:ndelta]
    delta <- pars[nbeta + 1:ndelta]
    ell <- pars[nbeta + ndelta + 1:nell]
    L <- matrix(0, ndelta, ndelta)
    L[lower.tri(L, diag = TRUE)] <- ell
    elldiag <- diag(L)
    omega <- crossprod(L)
    y <- as.matrix(xmat)%*%beta + as.matrix(smat)%*%delta
  } else {
    y <- as.matrix(xmat)%*%beta
  }
  lambda <- exp(y)
  xlam <- xmat*c(lambda)
  xxlam <- crossprod(xmat, xlam)
  dbetadbeta <- - xxlam - diag(1/betavar, nbeta)
  if(ndelta > 0){
    slam <- smat*c(lambda)
    sslam <- crossprod(smat, slam)
    xslam <- crossprod(xmat, slam)
    ddeltaddelta <- - sslam - omega
    dbetaddelta <- - xslam
    K <- commutation.matrix(ndelta)
    M <- elimination.matrix(ndelta)
    N2 <- 2*N.matrix(ndelta)
    Indelta <- diag(ndelta)
    ddeltadvell <- -kronecker(t(delta), Indelta)%*%N2%*%tcrossprod(kronecker(L, Indelta), M)
    R <- matrix(0, nell, nell)
    diag(R)[(ndelta+1)*(1:ndelta) - (1:ndelta)*(1 + 1:ndelta)/2 + 1:ndelta - ndelta] <- (ndelta + 2 - 1:ndelta)
    MKprime <- tcrossprod(M, K)
    dvelldvell <- - (MKprime%*%tcrossprod(kronecker(tcrossprod(delta), Indelta), MKprime) +
                     diag(nell)) +
      - (M%*%tcrossprod(kronecker(Indelta, ominvscale), M)) - R%*%diag(1/ell^2)
    out <- matrix(0, npars, npars)
    out[1:nbeta, 1:nbeta] <- dbetadbeta
    out[nbeta + 1:ndelta, nbeta + 1:ndelta] <- ddeltaddelta
    out[1:nbeta, nbeta + 1:ndelta] <- dbetaddelta
    out[nbeta + 1:ndelta, 1:nbeta] <- t(dbetaddelta)
    out[nbeta + ndelta + 1:nell, nbeta + ndelta + 1:nell] <- dvelldvell
    out[nbeta + 1:ndelta, nbeta + ndelta + 1:nell] <- ddeltadvell
    out[nbeta + ndelta + 1:nell, nbeta + 1:ndelta] <- t(ddeltadvell)
  } else {
    out <- dbetadbeta
  }
  return(out)
}

poisfullindmetrop <- function(niter, init, mu, lpbest, df, dat, nbeta, ndelta,
                              betamn = 0, betavar = 100,
                              omdf = ndelta + 1, ominvscale = diag(1, ndelta)){
  z <- dat[, 1]
  xmat <- dat[, 1 + 1:nbeta]
  smat <- dat[, 1 + nbeta + 1:ndelta]
  npar <- length(init)
  par <- init
  draws <- matrix(0, ncol = npar, nrow = niter + 1)
  nell <- ndelta*(ndelta + 1)/2
  omeganums <- apply(which(lower.tri(diag(ndelta), TRUE), TRUE), 1, paste, collapse = ".")
  colnames(draws) <- c(paste("beta", 1:nbeta, sep = "."),
                       paste("delta", 1:ndelta, sep = "."),
                       paste("omega", omeganums, sep = "."))
  acc <- rep(0, niter)
  lpbests <- rep(0, niter)
  draws[1,] <- init
  hess <- poislposthessfull(mu, dat, nbeta, ndelta, betamn, betavar, omdf, ominvscale)
  Sigma <- chol2inv(chol(-hess))
  Rsig <- chol(Sigma)
  lpold <- lpprop <- poislpostfull(par, dat, nbeta, ndelta, betamn, betavar, omdf, ominvscale)
  ljold <- dmt(par, mu, Sigma, df, TRUE)
  for(iter in 1:niter){
    if(iter %% 50 == 0){
          print(iter)
    }
    prop <- rmtfixed(1, mu, Sigma, df, Rsig)
    lpprop <- poislpostfull(prop, dat, nbeta, ndelta, betamn, betavar, omdf, ominvscale)
    ljprop <- dmt(prop, mu, Sigma, df, TRUE)
    la <- lpprop - lpold + ljold - ljprop
    u <- runif(1)
    if(log(u) < la){
      par <- prop
      acc[iter] <- 1
      lpold <- lpprop
      ljold <- ljprop
    }
    if(lpbest < lpprop){
      lpbest <- lpprop
      mu <- par
      hess <- poislposthessfull(mu, dat, nbeta, ndelta, betamn, betavar, omdf, ominvscale)
      Sigma <- chol2inv(chol(-hess))
      Rsig <- chol(Sigma)
    }
    lpbests[iter] <- lpbest
    draws[iter + 1,] <- par
  }
  return(list(draws = draws, acc = acc, state = par, lpbests = lpbests, lpbest = lpbest,
              mu = mu, Rsig = Rsig, Sigma = Sigma))
}
