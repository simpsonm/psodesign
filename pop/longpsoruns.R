source("popfun.R")
source("mcmcfun.R")
source("psofun.R")
load("popdat/popdat.RData")
library(rstan)
library(MCMCpack)



nbeta <- 1
nswarms <- c(100)
ndeltasiid <- c(10, 30, 50)
ndeltasfull <- c(5, 10, 15)
psoiter <- 10000
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496
alpha <- .2*niter
beta <- 1

psolist <- list()
bbpsoxpmclist <- list()
atbbpsoxpmc1list <- list()
atbbpsoxpmc2list <- list()
psolist[["pois"]] <- list()
bbpsoxpmclist[["pois"]] <- list()
atbbpsoxpmc1list[["pois"]] <- list()
atbbpsoxpmc2list[["pois"]] <- list()
psolist[["lnorm"]] <- list()
bbpsoxpmclist[["lnorm"]] <- list()
atbbpsoxpmc1list[["lnorm"]] <- list()
atbbpsoxpmc2list[["lnorm"]] <- list()

psolist[["pois"]][["iid"]] <- list()
bbpsoxpmclist[["pois"]][["iid"]] <- list()
atbbpsoxpmc1list[["pois"]][["iid"]] <- list()
atbbpsoxpmc2list[["pois"]][["iid"]] <- list()
psolist[["lnorm"]][["iid"]] <- list()
bbpsoxpmclist[["lnorm"]][["iid"]] <- list()
atbbpsoxpmc1list[["lnorm"]][["iid"]] <- list()
atbbpsoxpmc2list[["lnorm"]][["iid"]] <- list()
psolist[["pois"]][["full"]] <- list()
bbpsoxpmclist[["pois"]][["full"]] <- list()
atbbpsoxpmc1list[["pois"]][["full"]] <- list()
atbbpsoxpmc2list[["pois"]][["full"]] <- list()
psolist[["lnorm"]][["full"]] <- list()
bbpsoxpmclist[["lnorm"]][["full"]] <- list()
atbbpsoxpmc1list[["lnorm"]][["full"]] <- list()
atbbpsoxpmc2list[["lnorm"]][["full"]] <- list()

{
for(nswarm in nswarms){
  i <- which(nswarm == nswarms)
  nbhd <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)})
  psolist[["pois"]][["iid"]][[i]] <- list()
  bbpsoxpmclist[["pois"]][["iid"]][[i]] <- list()
  atbbpsoxpmc1list[["pois"]][["iid"]][[i]] <- list()
  atbbpsoxpmc2list[["pois"]][["iid"]][[i]] <- list()
  psolist[["lnorm"]][["iid"]][[i]] <- list()
  bbpsoxpmclist[["lnorm"]][["iid"]][[i]] <- list()
  atbbpsoxpmc1list[["lnorm"]][["iid"]][[i]] <- list()
  atbbpsoxpmc2list[["lnorm"]][["iid"]][[i]] <- list()
  psolist[["pois"]][["full"]][[i]] <- list()
  bbpsoxpmclist[["pois"]][["full"]][[i]] <- list()
  atbbpsoxpmc1list[["pois"]][["full"]][[i]] <- list()
  atbbpsoxpmc2list[["pois"]][["full"]][[i]] <- list()
  psolist[["lnorm"]][["full"]][[i]] <- list()
  bbpsoxpmclist[["lnorm"]][["full"]][[i]] <- list()
  atbbpsoxpmc1list[["lnorm"]][["full"]][[i]] <- list()
  atbbpsoxpmc2list[["lnorm"]][["full"]][[i]] <- list()
  for(j in 1:3){
    print(c(i,j))
    ndeltaiid <- ndeltasiid[j]
    ndeltafull <- ndeltasfull[j]
    nell <- ndeltafull*(ndeltafull + 1)/2
    poispsofullinit <- matrix(runif((ndeltafull + nbeta + nell)*nswarm, -1, 1), ncol = nswarm)
    lnormpsofullinit <- matrix(runif((ndeltafull + nbeta + nell + 1)*nswarm, -1, 1),
                               ncol = nswarm)
    ominvscale <- diag(ndeltafull)
    K <- commutation.matrix(ndeltafull)
    M <- elimination.matrix(ndeltafull)
    N2 <- 2*N.matrix(ndeltafull)
    R <- matrix(0, nell, nell)
    MKprime <- tcrossprod(M, K)
    diag(R)[(ndeltafull+1)*(1:ndeltafull) - (1:ndeltafull)*(1 + 1:ndeltafull)/2 +
            1:ndeltafull - ndeltafull] <- (ndeltafull + 2 - 1:ndeltafull)
    dvelldvell <- - M%*%tcrossprod(kronecker(diag(ndeltafull), ominvscale), M)
    datlistfull <- list(dat = popdat, nbeta = nbeta, ndelta = ndeltafull,
                    betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                    aphi = 1, bphi = 1, omdf = ndeltafull + 1, ominvscale = diag(ndeltafull),
                    K = K, M = M, N2 = N2, R = R, MKprime = MKprime, dvelldvell = dvelldvell)
    poispsoiidinit <- matrix(runif((ndeltaiid + nbeta + 1)*nswarm, -1, 1), ncol = nswarm)
    lnormpsoiidinit <- matrix(runif((ndeltaiid + nbeta + 1 + 1)*nswarm, -1, 1), ncol = nswarm)
    datlistiid <- list(dat = popdat, nbeta = nbeta, ndelta = ndeltaiid,
                    betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                    aphi = 1, bphi = 1)
    print("Poisson IID")
    psolist[["pois"]][["iid"]][[i]][[j]] <-
      pso(psoiter, nswarm, inertia, cognitive, social, poispsoiidinit, nbhd,
          poislpostiid, datlistiid)
    bbpsoxpmclist[["pois"]][["iid"]][[i]][[j]] <-
      bbpso(psoiter, nswarm, 1, 0, poispsoiidinit, nbhd, poislpostiid, Inf, FALSE,
            c(0.5, 0.5), datlistiid)
    atbbpsoxpmc1list[["pois"]][["iid"]][[i]][[j]] <-
      bbpso(psoiter, nswarm, 1, 0.5, poispsoiidinit, nbhd, poislpostiid, 1, TRUE,
            c(0.5, 0.5), datlistiid)
    atbbpsoxpmc2list[["pois"]][["iid"]][[i]][[j]] <-
      bbpso(psoiter, nswarm, 1, 0.3, poispsoiidinit, nbhd, poislpostiid, 1, TRUE,
            c(0.5, 0.5), datlistiid)
    print("Poisson Full")
    psolist[["pois"]][["full"]][[i]][[j]] <-
      pso(psoiter, nswarm, inertia, cognitive, social, poispsofullinit, nbhd,
          poislpostfull, datlistfull)
    bbpsoxpmclist[["pois"]][["full"]][[i]][[j]] <-
      bbpso(psoiter, nswarm, 1, 0, poispsofullinit, nbhd, poislpostfull, Inf, FALSE,
            c(0.5, 0.5), datlistfull)
    atbbpsoxpmc1list[["pois"]][["full"]][[i]][[j]] <-
      bbpso(psoiter, nswarm, 1, 0.5, poispsofullinit, nbhd, poislpostfull, 1, TRUE,
            c(0.5, 0.5), datlistfull)
    atbbpsoxpmc2list[["pois"]][["full"]][[i]][[j]] <-
      bbpso(psoiter, nswarm, 1, 0.3, poispsofullinit, nbhd, poislpostfull, 1, TRUE,
            c(0.5, 0.5), datlistfull)
    print("lognormal IID")
    psolist[["lnorm"]][["iid"]][[i]][[j]] <-
      pso(psoiter, nswarm, inertia, cognitive, social, lnormpsoiidinit, nbhd,
          lnormlpostiid, datlistiid)
    bbpsoxpmclist[["lnorm"]][["iid"]][[i]][[j]] <-
      bbpso(psoiter, nswarm, 1, 0, lnormpsoiidinit, nbhd, lnormlpostiid, Inf, FALSE,
            c(0.5, 0.5), datlistiid)
    atbbpsoxpmc1list[["lnorm"]][["iid"]][[i]][[j]] <-
      bbpso(psoiter, nswarm, 1, 0.5, lnormpsoiidinit, nbhd, lnormlpostiid, 1, TRUE,
            c(0.5, 0.5), datlistiid)
    atbbpsoxpmc2list[["lnorm"]][["iid"]][[i]][[j]] <-
      bbpso(psoiter, nswarm, 1, 0.3, lnormpsoiidinit, nbhd, lnormlpostiid, 1, TRUE,
            c(0.5, 0.5), datlistiid)
    print("lognormal Full")
    psolist[["lnorm"]][["full"]][[i]][[j]] <-
      pso(psoiter, nswarm, inertia, cognitive, social, lnormpsofullinit, nbhd,
          lnormlpostfull, datlistfull)
    bbpsoxpmclist[["lnorm"]][["full"]][[i]][[j]] <-
      bbpso(psoiter, nswarm, 1, 0, lnormpsofullinit, nbhd, lnormlpostfull, Inf, FALSE,
            c(0.5, 0.5), datlistfull)
    atbbpsoxpmc1list[["lnorm"]][["full"]][[i]][[j]] <-
      bbpso(psoiter, nswarm, 1, 0.5, lnormpsofullinit, nbhd, lnormlpostfull, 1, TRUE,
            c(0.5, 0.5), datlistfull)
    atbbpsoxpmc2list[["lnorm"]][["full"]][[i]][[j]] <-
      bbpso(psoiter, nswarm, 1, 0.3, lnormpsofullinit, nbhd, lnormlpostfull, 1, TRUE,
            c(0.5, 0.5), datlistfull)
    save(psolist, file = "psolist.RData")
    save(bbpsoxpmclist, file = "bbpsoxpmclist.RData")
    save(atbbpsoxpmc1list, file = "atbbpsoxpmc1list.RData")
    save(atbbpsoxpmc2list, file = "atbbpsoxpmc2list.RData")
  }
}
}


dipsolist <- list()
atpso1list <- list()
atpso2list <- list()
atpso1list[["pois"]] <- list()
atpso2list[["pois"]] <- list()
atpso1list[["lnorm"]] <- list()
atpso2list[["lnorm"]] <- list()
dipsolist[["pois"]] <- list()
dipsolist[["lnorm"]] <- list()

atpso1list[["pois"]][["iid"]] <- list()
atpso1list[["lnorm"]][["iid"]] <- list()
atpso1list[["pois"]][["full"]] <- list()
atpso1list[["lnorm"]][["full"]] <- list()

atpso2list[["pois"]][["iid"]] <- list()
atpso2list[["lnorm"]][["iid"]] <- list()
atpso2list[["pois"]][["full"]] <- list()
atpso2list[["lnorm"]][["full"]] <- list()

dipsolist[["pois"]][["iid"]] <- list()
dipsolist[["lnorm"]][["iid"]] <- list()
dipsolist[["pois"]][["full"]] <- list()
dipsolist[["lnorm"]][["full"]] <- list()


{
for(nswarm in nswarms){
  i <- which(nswarm == nswarms)
  nbhd <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)})
  atpso1list[["pois"]][["iid"]][[i]] <- list()
  atpso1list[["lnorm"]][["iid"]][[i]] <- list()
  atpso1list[["pois"]][["full"]][[i]] <- list()
  atpso1list[["lnorm"]][["full"]][[i]] <- list()
  atpso2list[["pois"]][["iid"]][[i]] <- list()
  atpso2list[["lnorm"]][["iid"]][[i]] <- list()
  atpso2list[["pois"]][["full"]][[i]] <- list()
  atpso2list[["lnorm"]][["full"]][[i]] <- list()
  dipsolist[["pois"]][["iid"]][[i]] <- list()
  dipsolist[["lnorm"]][["iid"]][[i]] <- list()
  dipsolist[["pois"]][["full"]][[i]] <- list()
  dipsolist[["lnorm"]][["full"]][[i]] <- list()
  for(j in 1:3){
    print(c(i,j))
    ndeltaiid <- ndeltasiid[j]
    ndeltafull <- ndeltasfull[j]
    nell <- ndeltafull*(ndeltafull + 1)/2
    poispsofullinit <- matrix(runif((ndeltafull + nbeta + nell)*nswarm, -1, 1), ncol = nswarm)
    lnormpsofullinit <- matrix(runif((ndeltafull + nbeta + nell + 1)*nswarm, -1, 1),
                               ncol = nswarm)
    ominvscale <- diag(ndeltafull)
    K <- commutation.matrix(ndeltafull)
    M <- elimination.matrix(ndeltafull)
    N2 <- 2*N.matrix(ndeltafull)
    R <- matrix(0, nell, nell)
    MKprime <- tcrossprod(M, K)
    diag(R)[(ndeltafull+1)*(1:ndeltafull) - (1:ndeltafull)*(1 + 1:ndeltafull)/2 +
            1:ndeltafull - ndeltafull] <- (ndeltafull + 2 - 1:ndeltafull)
    dvelldvell <- - M%*%tcrossprod(kronecker(diag(ndeltafull), ominvscale), M)
    datlistfull <- list(dat = popdat, nbeta = nbeta, ndelta = ndeltafull,
                    betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                    aphi = 1, bphi = 1, omdf = ndeltafull + 1, ominvscale = diag(ndeltafull),
                    K = K, M = M, N2 = N2, R = R, MKprime = MKprime, dvelldvell = dvelldvell)
    poispsoiidinit <- matrix(runif((ndeltaiid + nbeta + 1)*nswarm, -1, 1), ncol = nswarm)
    lnormpsoiidinit <- matrix(runif((ndeltaiid + nbeta + 1 + 1)*nswarm, -1, 1), ncol = nswarm)
    datlistiid <- list(dat = popdat, nbeta = nbeta, ndelta = ndeltaiid,
                    betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                    aphi = 1, bphi = 1)
    print("Poisson IID")
    atpso1list[["pois"]][["iid"]][[i]][[j]] <-
      pso(psoiter, nswarm, inertia, cognitive, social, poispsoiidinit, nbhd,
          poislpostiid, datlistiid, tune = TRUE, rate = 0.3, ccc = 0.1)
    atpso2list[["pois"]][["iid"]][[i]][[j]] <-
      pso(psoiter, nswarm, inertia, cognitive, social, poispsoiidinit, nbhd,
          poislpostiid, datlistiid, tune = TRUE, rate = 0.5, ccc = 0.1)
    dipsolist[["pois"]][["iid"]][[i]][[j]] <-
      pso(psoiter, nswarm, inertia, cognitive, social, poispsoiidinit, nbhd,
          poislpostiid, datlistiid, tune = TRUE, style = "deterministic", rate = 0.5, ccc = 0.1)
    print("Poisson Full")
    atpso1list[["pois"]][["full"]][[i]][[j]] <-
      pso(psoiter, nswarm, inertia, cognitive, social, poispsofullinit, nbhd,
          poislpostfull, datlistfull, tune = TRUE, rate = 0.3, ccc = 0.1)
    atpso2list[["pois"]][["full"]][[i]][[j]] <-
      pso(psoiter, nswarm, inertia, cognitive, social, poispsofullinit, nbhd,
          poislpostfull, datlistfull, tune = TRUE, rate = 0.5, ccc = 0.1)
    dipsolist[["pois"]][["full"]][[i]][[j]] <-
      pso(psoiter, nswarm, inertia, cognitive, social, poispsofullinit, nbhd,
          poislpostfull, datlistfull, tune = TRUE, style = "deterministic", rate = 0.5, ccc = 0.1)
    print("lognormal IID")
    atpso1list[["lnorm"]][["iid"]][[i]][[j]] <-
      pso(psoiter, nswarm, inertia, cognitive, social, lnormpsoiidinit, nbhd,
          lnormlpostiid, datlistiid, tune = TRUE, rate = 0.3, ccc = 0.1)
    atpso2list[["lnorm"]][["iid"]][[i]][[j]] <-
      pso(psoiter, nswarm, inertia, cognitive, social, lnormpsoiidinit, nbhd,
          lnormlpostiid, datlistiid, tune = TRUE, rate = 0.5, ccc = 0.1)
    print("lognormal Full")
    atpso1list[["lnorm"]][["full"]][[i]][[j]] <-
      pso(psoiter, nswarm, inertia, cognitive, social, lnormpsofullinit, nbhd,
          lnormlpostfull, datlistfull, tune = TRUE, rate = 0.3, ccc = 0.1)
    atpso2list[["lnorm"]][["full"]][[i]][[j]] <-
      pso(psoiter, nswarm, inertia, cognitive, social, lnormpsofullinit, nbhd,
          lnormlpostfull, datlistfull, tune = TRUE, rate = 0.5, ccc = 0.1)
    save(atpso1list, file = "atpso1list.RData")
    save(atpso2list, file = "atpso2list.RData")
  }
}
}


## i indexes swarm size, j indexes number of random effects

load("psolist.RData")
load("bbpsoxpmclist.RData")
load("atbbpsoxpmc1list.RData")
load("atbbpsoxpmc2list.RData")
load("atpso1list.RData")
load("atpso2list.RData")

## pois/lnorm full/iid nswarm ndelta
mcmciter <- 10000
niters <- c(100, 1000, 5000, 10000)
ndeltasiid <- c(10, 30, 50)
ndeltasfull <- c(5, 10, 15)
nswarms <- c(100)
ranefs <- c("iid", "full")
models <- c("pois", "lnorm")
nbeta <- 1
longaccout <- NULL
df <- 100

for(ndeltak in 1:3){
  ndeltaiid <- ndeltasiid[ndeltak]
  ndeltafull <- ndeltasfull[ndeltak]
  nell <- ndeltafull*(ndeltafull + 1)/2
  ominvscale <- diag(ndeltafull)
  K <- commutation.matrix(ndeltafull)
  M <- elimination.matrix(ndeltafull)
  N2 <- 2*N.matrix(ndeltafull)
  R <- matrix(0, nell, nell)
  MKprime <- tcrossprod(M, K)
  diag(R)[(ndeltafull+1)*(1:ndeltafull) - (1:ndeltafull)*(1 + 1:ndeltafull)/2 +
          1:ndeltafull - ndeltafull] <- (ndeltafull + 2 - 1:ndeltafull)
  dvelldvell <- - M%*%tcrossprod(kronecker(diag(ndeltafull), ominvscale), M)
  datlistfull <- list(dat = popdat, nbeta = nbeta, ndelta = ndeltafull,
                      betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                      aphi = 1, bphi = 1, omdf = ndeltafull + 1, ominvscale = diag(ndeltafull),
                      K = K, M = M, N2 = N2, R = R, MKprime = MKprime, dvelldvell = dvelldvell)
  datlistiid <- list(dat = popdat, nbeta = nbeta, ndelta = ndeltaiid,
                     betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                     aphi = 1, bphi = 1)
  for(nswarmk in c(1)){
    nswarm <- nswarms[nswarmk]
    for(niter in niters){
      for(ranef in ranefs){
        for(model in models){
          print(c(ndeltak, nswarm, niter, ranef, model))
          datlist <- switch(ranef, iid = datlistiid, full = datlistfull)
          lpost <- switch(paste(model, ranef, sep=""),
                          poisiid = poislpostiid,
                          lnormiid = lnormlpostiid,
                          poisfull = poislpostfull,
                          lnormfull = lnormlpostfull)
          lposthess <- switch(paste(model, ranef, sep=""),
                              poisiid = poislposthessiid,
                              lnormiid = lnormlposthessiid,
                              poisfull = poislposthessfull,
                              lnormfull = lnormlposthessfull)
          ndelta <- switch(ranef, iid = ndeltaiid, full = ndeltafull)
          indmetropwithingibbs <- switch(paste(model, ranef, sep=""),
                                         poisiid = poisiidindgibbs,
                                         poisfull = poisfullindgibbs,
                                         lnormiid = lnormiidindgibbs,
                                         lnormfull = lnormfullindgibbs)
          mu <- psolist[[model]][[ranef]][[nswarmk]][[ndeltak]]$argmaxes[,niter+1]
          lpbest <- psolist[[model]][[ranef]][[nswarmk]][[ndeltak]]$maxes[niter+1]
          out <- indmetrop(mcmciter, lpost, lposthess, mu, mu, 100, lpbest, datlist)
          accout <- mean(out$acc[-c(1:10000)])
          accout <- mean(out$acc)
          longaccout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                         nswarm = nswarm, niter = niter,
                                         pso = "PSO", mcmc = "IMH", acc = accout),
                              longaccout)
          out <- indmetropwithingibbs(mcmciter, mu, mu, df, lpbest, datlist)
          accout <- mean(out$acc[-c(1:10000)])
          accout <- mean(out$acc)
          longaccout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                         nswarm = nswarm, niter = niter,
                                         pso = "PSO", mcmc = "IMHwG", acc = accout),
                              longaccout)
          mu <- bbpsoxpmclist[[model]][[ranef]][[nswarmk]][[ndeltak]]$argmaxes[,niter+1]
          lpbest <- bbpsoxpmclist[[model]][[ranef]][[nswarmk]][[ndeltak]]$maxes[niter+1]
          out <- indmetrop(mcmciter, lpost, lposthess, mu, mu, 100, lpbest, datlist)
          accout <- mean(out$acc[-c(1:10000)])
          accout <- mean(out$acc)
          longaccout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                         nswarm = nswarm, niter = niter,
                                         pso = "BBPSOxp-MC", mcmc = "IMH", acc = accout),
                              longaccout)
          out <- indmetropwithingibbs(mcmciter, mu, mu, df, lpbest, datlist)
          accout <- mean(out$acc[-c(1:10000)])
          accout <- mean(out$acc)
          longaccout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                         nswarm = nswarm, niter = niter,
                                         pso = "BBPSOxp-MC", mcmc = "IMHwG", acc = accout),
                              longaccout)
          mu <- atbbpsoxpmc1list[[model]][[ranef]][[nswarmk]][[ndeltak]]$argmaxes[,niter+1]
          lpbest <- atbbpsoxpmc1list[[model]][[ranef]][[nswarmk]][[ndeltak]]$maxes[niter+1]
          out <- indmetrop(mcmciter, lpost, lposthess, mu, mu, 100, lpbest, datlist)
          accout <- mean(out$acc[-c(1:10000)])
          accout <- mean(out$acc)
          longaccout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                         nswarm = nswarm, niter = niter,
                                         pso = "AT-BBPSOxp-MC-1", mcmc = "IMH", acc = accout),
                              longaccout)
          out <- indmetropwithingibbs(mcmciter, mu, mu, df, lpbest, datlist)
          accout <- mean(out$acc[-c(1:10000)])
          accout <- mean(out$acc)
          longaccout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                         nswarm = nswarm, niter = niter,
                                         pso = "AT-BBPSOxp-MC-1", mcmc = "IMHwG", acc = accout),
                              longaccout)
          mu <- atbbpsoxpmc2list[[model]][[ranef]][[nswarmk]][[ndeltak]]$argmaxes[,niter+1]
          lpbest <- atbbpsoxpmc2list[[model]][[ranef]][[nswarmk]][[ndeltak]]$maxes[niter+1]
          out <- indmetrop(mcmciter, lpost, lposthess, mu, mu, 100, lpbest, datlist)
          accout <- mean(out$acc[-c(1:10000)])
          accout <- mean(out$acc)
          longaccout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                         nswarm = nswarm, niter = niter,
                                         pso = "AT-BBPSOxp-MC-2", mcmc = "IMH", acc = accout),
                              longaccout)
          out <- indmetropwithingibbs(mcmciter, mu, mu, df, lpbest, datlist)
          accout <- mean(out$acc[-c(1:10000)])
          accout <- mean(out$acc)
          longaccout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                         nswarm = nswarm, niter = niter,
                                         pso = "AT-BBPSOxp-MC-2", mcmc = "IMHwG", acc = accout),
                              longaccout)
        }
      }
    }
  }
}

save(longaccout, file = "longaccout.RData")

## AT-PSO MCMC

load("longaccout.RData")

for(ndeltak in 1:3){
  ndeltaiid <- ndeltasiid[ndeltak]
  ndeltafull <- ndeltasfull[ndeltak]
  nell <- ndeltafull*(ndeltafull + 1)/2
  ominvscale <- diag(ndeltafull)
  K <- commutation.matrix(ndeltafull)
  M <- elimination.matrix(ndeltafull)
  N2 <- 2*N.matrix(ndeltafull)
  R <- matrix(0, nell, nell)
  MKprime <- tcrossprod(M, K)
  diag(R)[(ndeltafull+1)*(1:ndeltafull) - (1:ndeltafull)*(1 + 1:ndeltafull)/2 +
          1:ndeltafull - ndeltafull] <- (ndeltafull + 2 - 1:ndeltafull)
  dvelldvell <- - M%*%tcrossprod(kronecker(diag(ndeltafull), ominvscale), M)
  datlistfull <- list(dat = popdat, nbeta = nbeta, ndelta = ndeltafull,
                      betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                      aphi = 1, bphi = 1, omdf = ndeltafull + 1, ominvscale = diag(ndeltafull),
                      K = K, M = M, N2 = N2, R = R, MKprime = MKprime, dvelldvell = dvelldvell)
  datlistiid <- list(dat = popdat, nbeta = nbeta, ndelta = ndeltaiid,
                     betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                     aphi = 1, bphi = 1)
  for(nswarmk in c(1)){
    nswarm <- nswarms[nswarmk]
    for(niter in niters){
      for(ranef in ranefs){
        for(model in models){
          print(c(ndeltak, nswarm, niter, ranef, model))
          datlist <- switch(ranef, iid = datlistiid, full = datlistfull)
          lpost <- switch(paste(model, ranef, sep=""),
                          poisiid = poislpostiid,
                          lnormiid = lnormlpostiid,
                          poisfull = poislpostfull,
                          lnormfull = lnormlpostfull)
          lposthess <- switch(paste(model, ranef, sep=""),
                              poisiid = poislposthessiid,
                              lnormiid = lnormlposthessiid,
                              poisfull = poislposthessfull,
                              lnormfull = lnormlposthessfull)
          ndelta <- switch(ranef, iid = ndeltaiid, full = ndeltafull)
          indmetropwithingibbs <- switch(paste(model, ranef, sep=""),
                                         poisiid = poisiidindgibbs,
                                         poisfull = poisfullindgibbs,
                                         lnormiid = lnormiidindgibbs,
                                         lnormfull = lnormfullindgibbs)
          mu <- atpso1list[[model]][[ranef]][[nswarmk]][[ndeltak]]$argmaxes[,niter+1]
          lpbest <- atpso1list[[model]][[ranef]][[nswarmk]][[ndeltak]]$maxes[niter+1]
          try({out <- indmetrop(mcmciter, lpost, lposthess, mu, mu, 100, lpbest, datlist)
          accout <- mean(out$acc[-c(1:10000)])
          accout <- mean(out$acc)
          longaccout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                         nswarm = nswarm, niter = niter,
                                         pso = "AT-PSO-0.3", mcmc = "IMH", acc = accout),
                              longaccout)})
          try({out <- indmetropwithingibbs(mcmciter, mu, mu, df, lpbest, datlist)
          accout <- mean(out$acc[-c(1:10000)])
          accout <- mean(out$acc)
          longaccout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                         nswarm = nswarm, niter = niter,
                                         pso = "AT-PSO-0.3", mcmc = "IMHwG", acc = accout),
                              longaccout)})
          mu <- atpso2list[[model]][[ranef]][[nswarmk]][[ndeltak]]$argmaxes[,niter+1]
          lpbest <- atpso2list[[model]][[ranef]][[nswarmk]][[ndeltak]]$maxes[niter+1]
          try({out <- indmetrop(mcmciter, lpost, lposthess, mu, mu, 100, lpbest, datlist)
          accout <- mean(out$acc[-c(1:10000)])
          accout <- mean(out$acc)
          longaccout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                         nswarm = nswarm, niter = niter,
                                         pso = "AT-PSO-0.5", mcmc = "IMH", acc = accout),
                              longaccout)})
          try({out <- indmetropwithingibbs(mcmciter, mu, mu, df, lpbest, datlist)
          accout <- mean(out$acc[-c(1:10000)])
          accout <- mean(out$acc)
          longaccout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                         nswarm = nswarm, niter = niter,
                                         pso = "AT-PSO-0.5", mcmc = "IMHwG", acc = accout),
                              longaccout)})
        }
      }
    }
  }
}

save(longaccout, file = "longaccout.RData")






#### some experiments
load("psolist.RData")

par(mfrow=c(2,2))
models <- c("pois", "lnorm")
ranefs <- c("iid", "full")
k <- 3
t1s <- c(1000, 9000)
t2s <- c(1300, 10000)
for(i in 1:2){
  for(j in 1:2){
    model <- models[i]
    ranef <- ranefs[j]
    t1 <- t1s[j]
    t2 <- t2s[j]
    x <- psolist[[model]][[ranef]][[1]][[k]]$maxes[t1:t2]
    y <- atpso1list[[model]][[ranef]][[1]][[k]]$maxes[t1:t2]
    z <- atpso2list[[model]][[ranef]][[1]][[k]]$maxes[t1:t2]
    plot(ts(x), ylim=c(min(c(x,y,z)), max(c(x,y,z))))
    lines(1:(t2 - t1 + 1), y, col = "red")
    lines(1:(t2 - t1 + 1), z, col = "blue")
  }
}

## AT-PSO looks good by the above!!!!!!!!

nswarm <- 100
nbhd <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)})
ndeltafull <- 2
nell <- ndeltafull*(ndeltafull + 1)/2
poispsofullinit <- matrix(runif((ndeltafull + nbeta + nell)*nswarm, -1, 1), ncol = nswarm)
lnormpsofullinit <- matrix(runif((ndeltafull + nbeta + nell + 1)*nswarm, -1, 1),
                               ncol = nswarm)
ominvscale <- diag(ndeltafull)
K <- commutation.matrix(ndeltafull)
M <- elimination.matrix(ndeltafull)
N2 <- 2*N.matrix(ndeltafull)
R <- matrix(0, nell, nell)
MKprime <- tcrossprod(M, K)
diag(R)[(ndeltafull+1)*(1:ndeltafull) - (1:ndeltafull)*(1 + 1:ndeltafull)/2 +
        1:ndeltafull - ndeltafull] <- (ndeltafull + 2 - 1:ndeltafull)
dvelldvell <- - M%*%tcrossprod(kronecker(diag(ndeltafull), ominvscale), M)
datlistfull <- list(dat = popdat, nbeta = nbeta, ndelta = ndeltafull,
                    betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                    aphi = 1, bphi = 1, omdf = ndeltafull + 1, ominvscale = diag(ndeltafull),
                    K = K, M = M, N2 = N2, R = R, MKprime = MKprime, dvelldvell = dvelldvell)

psotest <- pso(psoiter, nswarm, inertia, cognitive, social, lnormpsofullinit, nbhd,
               lnormlpostfull, datlistfull)

mu <- psotest$argmax
lpbest <- psotest$max
ndelta <- ndeltafull

xsmat <- popdat[, 1 + 1:(nbeta + ndelta)]
xsxs <- crossprod(xsmat)
Ebeta <- diag(1/100, nbeta)
cbetadelta <- c(0/100, rep(0, ndelta))
lz <- log(popdat[,1])
xslz <- crossprod(xsmat, lz)

mubd <- mu[1:(nbeta + ndelta)]
muell <- mu[nbeta + ndelta + 1:nell]
mulphi2 <- mu[nbeta + ndelta + nell + 1]
phi2 <- exp(mulphi2)
L <- matrix(0, ndelta, ndelta)
L[lower.tri(L, TRUE)] <- muell
omega <- tcrossprod(L)

hess <- lnormlposthessfull(mu, datlistfull)
SigmaBig <- chol2inv(chol(-hess))
Sigma1 <- SigmaBig[1:(nbeta + ndelta), 1:(nbeta + ndelta)]
Sigma12 <- SigmaBig[1:(nbeta + ndelta), nbeta + ndelta + 0:nell + 1]
Sigma22 <- SigmaBig[nbeta + ndelta + 0:nell + 1, nbeta + ndelta + 0:nell + 1]
Sigma12Sigma22 <- Sigma12%*%chol2inv(chol(Sigma22))
Sigma <- Sigma1 - tcrossprod(Sigma12Sigma22, Sigma12)
Sigma <- (Sigma + t(Sigma))/2

omegabd <- xsxs/phi2 + bdiag(Ebeta, omega)
sigmagibbs <- chol2inv(chol(omegabd))
mugibbs <- drop(as.matrix(sigmagibbs%*%(xslz/phi2 + cbetadelta)))
mubd


mu2 <- mu
mu2[1:(nbeta + ndelta)] <- mugibbs
hess <- lnormlposthessfull(mu2, datlistfull)
SigmaBig <- chol2inv(chol(-hess))
Sigma1 <- SigmaBig[1:(nbeta + ndelta), 1:(nbeta + ndelta)]
Sigma12 <- SigmaBig[1:(nbeta + ndelta), nbeta + ndelta + 0:nell + 1]
Sigma22 <- SigmaBig[nbeta + ndelta + 0:nell + 1, nbeta + ndelta + 0:nell + 1]
Sigma12Sigma22 <- Sigma12%*%chol2inv(chol(Sigma22))
Sigma <- Sigma1 - tcrossprod(Sigma12Sigma22, Sigma12)
Sigma <- (Sigma + t(Sigma))/2


Sigma

(as.matrix(sigmagibbs) - Sigma)^2



out <- indmetrop(10000, lnormlpostfull, lnormlposthessfull, mu, mu, 100, lpbest, datlistfull,
                 tune = FALSE)

out2 <- lnormfullindgibbs(10000, mu, mu, 100, lpbest, datlistfull, FALSE)

mean(out$acc)
mean(out2$acc)
