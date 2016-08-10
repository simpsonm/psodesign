source("../psofun.R")
source("popfun.R")
source("mcmcfun.R")
load("popdat/popdat.RData")
library(MCMCpack)

nbeta <- 1
nswarm <- 50
ndeltasiid <- c(10, 30, 50)
ndeltasfull <- c(5, 10, 15)
psoiter <- 2000
niter <- psoiter
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496
alpha <- .2*psoiter
beta <- 1
rates <- c(0.3, 0.5)
ranefs <- c("iid", "full")
models <- c("pois", "lnorm")
nbhd <- sapply(1:nswarm, function(x){return( (x + -3:3 - 1)%%nswarm + 1)}) ## ring-3
ccc <- 0.1

psolist <- list()
bbpsolist <- list()
dipsolist <- list()
atpso1list <- list()
atpso2list <- list()
atbbpso1list <- list()
atbbpso2list <- list()

psolist[["pois"]] <- list()
bbpsolist[["pois"]] <- list()
dipsolist[["pois"]] <- list()
atpso1list[["pois"]] <- list()
atpso2list[["pois"]] <- list()
atbbpso1list[["pois"]] <- list()
atbbpso2list[["pois"]] <- list()

psolist[["lnorm"]] <- list()
bbpsolist[["lnorm"]] <- list()
dipsolist[["lnorm"]] <- list()
atpso1list[["lnorm"]] <- list()
atpso2list[["lnorm"]] <- list()
atbbpso1list[["lnorm"]] <- list()
atbbpso2list[["lnorm"]] <- list()

psolist[["pois"]][["iid"]] <- list()
bbpsolist[["pois"]][["iid"]] <- list()
dipsolist[["pois"]][["iid"]] <- list()
atpso1list[["pois"]][["iid"]] <- list()
atpso2list[["pois"]][["iid"]] <- list()
atbbpso1list[["pois"]][["iid"]] <- list()
atbbpso2list[["pois"]][["iid"]] <- list()

psolist[["pois"]][["full"]] <- list()
bbpsolist[["pois"]][["full"]] <- list()
dipsolist[["pois"]][["full"]] <- list()
atpso1list[["pois"]][["full"]] <- list()
atpso2list[["pois"]][["full"]] <- list()
atbbpso1list[["pois"]][["full"]] <- list()
atbbpso2list[["pois"]][["full"]] <- list()

psolist[["lnorm"]][["iid"]] <- list()
bbpsolist[["lnorm"]][["iid"]] <- list()
dipsolist[["lnorm"]][["iid"]] <- list()
atpso1list[["lnorm"]][["iid"]] <- list()
atpso2list[["lnorm"]][["iid"]] <- list()
atbbpso1list[["lnorm"]][["iid"]] <- list()
atbbpso2list[["lnorm"]][["iid"]] <- list()

psolist[["lnorm"]][["full"]] <- list()
bbpsolist[["lnorm"]][["full"]] <- list()
dipsolist[["lnorm"]][["full"]] <- list()
atpso1list[["lnorm"]][["full"]] <- list()
atpso2list[["lnorm"]][["full"]] <- list()
atbbpso1list[["lnorm"]][["full"]] <- list()
atbbpso2list[["lnorm"]][["full"]] <- list()


datlistiids <- list()
datlistfulls <- list()
for(idelta in 1:length(ndeltasiid)){
  ndeltafull <- ndeltasfull[idelta]
  ndeltaiid <- ndeltasiid[idelta]
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
  datlistiids[[idelta]] <- datlistiid
  datlistfulls[[idelta]] <- datlistfull
}

for(idelta in 1:length(ndeltasiid)){
  datlistiid <- datlistiids[[idelta]]
  datlistfull <- datlistfulls[[idelta]]
  ndeltafull <- ndeltasfull[idelta]
  ndeltaiid <- ndeltasiid[idelta]
  nell <- ndeltafull*(ndeltafull + 1)/2
  for(ranef in ranefs){
    datlist <- switch(ranef, iid = datlistiid, full = datlistfull)
    ndelta <- switch(ranef, iid = ndeltaiid, full = ndeltafull)
    for(model in models){
      print(c(model, ranef, idelta))
      lpost <- switch(paste(model, ranef, sep=""),
                      poisiid = poislpostiid,
                      lnormiid = lnormlpostiid,
                      poisfull = poislpostfull,
                      lnormfull = lnormlpostfull)
      npar <- nbeta + ndelta + (1 - nell)*(ranef == "iid") + nell + 1*(model == "lnorm")
      bfgsinit <- rep(0, npar)
      bfgsctr <- 0
      while(abs(lpost(bfgsinit, datlist)) == Inf){
        bfgsctr <- bfgsctr + 1
        bfgsinit <- rnorm(npar)
        if(bfgsctr %% 100){
          print("while loop counter is ", bfgsctr)
        }
      }
      bfgsout <- optim(bfgsinit, lpost, datlist = datlist, method = "BFGS",
                       control=list(fnscale=-1))
      mubfgs <- bfgsout$par
      init <- matrix(runif(npar*nswarm, -1, 1), ncol = nswarm) + mubfgs
      init[,1] <- mubfgs
      psolist[[model]][[ranef]][[idelta]] <-
        pso(niter, nswarm, inertia, cognitive, social, init, nbhd,
            lpost, datlist = datlist)
      bbpsolist[[model]][[ranef]][[idelta]] <-
        bbpso(niter, nswarm, 0, 1, init, nbhd, lpost, Inf, FALSE,
              c(.5,.5), datlist = datlist)
      dipsolist[[model]][[ranef]][[idelta]] <-
        pso(niter, nswarm, inertia, social, cognitive, init, nbhd,
            lpost, datlist = datlist, tune = TRUE, style = "deterministic",
            alpha = alpha, beta = beta)
      atpso1list[[model]][[ranef]][[idelta]] <-
        pso(niter, nswarm, 0.9, cognitive, social, init, nbhd,
            lpost, datlist = datlist, tune = TRUE, style = "adaptive",
            rate = rates[1], ccc = ccc)
      atpso2list[[model]][[ranef]][[idelta]] <-
        pso(niter, nswarm, 0.9, cognitive, social, init, nbhd,
            lpost, datlist = datlist, tune = TRUE, style = "adaptive",
            rate = rates[2], ccc = ccc)
      atbbpso1list[[model]][[ranef]][[idelta]] <-
        bbpso(niter, nswarm, 1, rates[1], init, nbhd, lpost, 1,
              TRUE, c(0.5,0.5), datlist = datlist, ccc = ccc)
      atbbpso2list[[model]][[ranef]][[idelta]] <-
        bbpso(niter, nswarm, 1, rates[2], init, nbhd, lpost, 1,
              TRUE, c(0.5,0.5), datlist = datlist, ccc = ccc)
      save(psolist, file = "psolist.RData")
      save(bbpsolist, file = "bbpsolist.RData")
      save(dipsolist, file = "dipsolist.RData")
      save(atpso1list, file = "atpso1list.RData")
      save(atpso2list, file = "atpso2list.RData")
      save(atbbpso1list, file = "atbbpso1list.RData")
      save(atbbpso2list, file = "atbbpso2list.RData")
    }
  }
}

## i indexes swarm size, j indexes number of random effects

load("psolist.RData")
load("bbpsolist.RData")
load("dipsolist.RData")
load("atbbpso1list.RData")
load("atbbpso2list.RData")
load("atpso1list.RData")
load("atpso2list.RData")

## pois/lnorm full/iid nswarm ndelta
mcmciter <- 10000
niters <- c(0, 100, 500, 1000, 1500, 2000)
ndeltasiid <- c(10, 30, 50)
ndeltasfull <- c(5, 10, 15)
ranefs <- c("iid", "full")
models <- c("pois", "lnorm")
nbeta <- 1
accout <- NULL
df <- 100

{
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
  for(niter in niters){
    for(ranef in ranefs){
      for(model in models){
        print(c(ndeltak, niter, ranef, model))
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
        mu <- psolist[[model]][[ranef]][[ndeltak]]$argmaxes[,niter+1]
        lpbest <- psolist[[model]][[ranef]][[ndeltak]]$maxes[niter+1]
        try({
          out <- indmetrop(mcmciter, lpost, lposthess, mu, mu, 100, lpbest, datlist)
          acc <- mean(out$acc)
          accout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                     nswarm = nswarm, niter = niter,
                                     pso = "PSO", mcmc = "IMH", acc = acc),
                          accout)
        })
        try({
          out <- indmetropwithingibbs(mcmciter, mu, mu, df, lpbest, datlist)
          acc <- mean(out$acc)
          accout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                     nswarm = nswarm, niter = niter,
                                     pso = "PSO", mcmc = "IMHwG", acc = acc),
                          accout)
        })
        mu <- dipsolist[[model]][[ranef]][[ndeltak]]$argmaxes[,niter+1]
        lpbest <- dipsolist[[model]][[ranef]][[ndeltak]]$maxes[niter+1]
        try({
          out <- indmetrop(mcmciter, lpost, lposthess, mu, mu, 100, lpbest, datlist)
          acc <- mean(out$acc)
          accout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                     nswarm = nswarm, niter = niter,
                                     pso = "DI-PSO", mcmc = "IMH", acc = acc),
                          accout)
        })
        try({                  
          out <- indmetropwithingibbs(mcmciter, mu, mu, df, lpbest, datlist)
          acc <- mean(out$acc)
          accout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                     nswarm = nswarm, niter = niter,
                                     pso = "DI-PSO", mcmc = "IMHwG", acc = acc),
                          accout)
        })
        mu <- bbpsolist[[model]][[ranef]][[ndeltak]]$argmaxes[,niter+1]
        lpbest <- bbpsolist[[model]][[ranef]][[ndeltak]]$maxes[niter+1]
        try({                  
          out <- indmetrop(mcmciter, lpost, lposthess, mu, mu, 100, lpbest, datlist)
          acc <- mean(out$acc)
          accout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                     nswarm = nswarm, niter = niter,
                                     pso = "BBPSO", mcmc = "IMH", acc = acc),
                          accout)
        })
        try({                  
          out <- indmetropwithingibbs(mcmciter, mu, mu, df, lpbest, datlist)
          acc <- mean(out$acc)
          accout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                     nswarm = nswarm, niter = niter,
                                     pso = "BBPSO", mcmc = "IMHwG", acc = acc),
                          accout)
        })
        mu <- atbbpso1list[[model]][[ranef]][[ndeltak]]$argmaxes[,niter+1]
        lpbest <- atbbpso1list[[model]][[ranef]][[ndeltak]]$maxes[niter+1]
        try({                          
          out <- indmetrop(mcmciter, lpost, lposthess, mu, mu, 100, lpbest, datlist)
          acc <- mean(out$acc)
          accout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                     nswarm = nswarm, niter = niter,
                                     pso = "AT-BBPSO-1", mcmc = "IMH", acc = acc),
                          accout)
        })
        try({                          
          out <- indmetropwithingibbs(mcmciter, mu, mu, df, lpbest, datlist)
          acc <- mean(out$acc)
          accout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                     nswarm = nswarm, niter = niter,
                                     pso = "AT-BBPSO-1", mcmc = "IMHwG", acc = acc),
                          accout)
        })
        mu <- atbbpso2list[[model]][[ranef]][[ndeltak]]$argmaxes[,niter+1]
        lpbest <- atbbpso2list[[model]][[ranef]][[ndeltak]]$maxes[niter+1]
        try({                                  
          out <- indmetrop(mcmciter, lpost, lposthess, mu, mu, 100, lpbest, datlist)
          acc <- mean(out$acc)
          accout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                     nswarm = nswarm, niter = niter,
                                     pso = "AT-BBPSO-2", mcmc = "IMH", acc = acc),
                          accout)
        })
        try({                                  
          out <- indmetropwithingibbs(mcmciter, mu, mu, df, lpbest, datlist)
          acc <- mean(out$acc)
          accout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                     nswarm = nswarm, niter = niter,
                                     pso = "AT-BBPSO-2", mcmc = "IMHwG", acc = acc),
                          accout)
        })
        mu <- atpso1list[[model]][[ranef]][[ndeltak]]$argmaxes[,niter+1]
        lpbest <- atpso1list[[model]][[ranef]][[ndeltak]]$maxes[niter+1]
        try({
          out <- indmetrop(mcmciter, lpost, lposthess, mu, mu, 100, lpbest, datlist)
          acc <- mean(out$acc)
          accout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                     nswarm = nswarm, niter = niter,
                                     pso = "AT-PSO-1", mcmc = "IMH", acc = acc),
                          accout)
        })
        try({
          out <- indmetropwithingibbs(mcmciter, mu, mu, df, lpbest, datlist)
          acc <- mean(out$acc)
          accout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                     nswarm = nswarm, niter = niter,
                                     pso = "AT-PSO-1", mcmc = "IMHwG", acc = acc),
                          accout)
        })
        mu <- atpso2list[[model]][[ranef]][[ndeltak]]$argmaxes[,niter+1]
        lpbest <- atpso2list[[model]][[ranef]][[ndeltak]]$maxes[niter+1]
        try({
          out <- indmetrop(mcmciter, lpost, lposthess, mu, mu, 100, lpbest, datlist)
          acc <- mean(out$acc)
          accout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                     nswarm = nswarm, niter = niter,
                                     pso = "AT-PSO-2", mcmc = "IMH", acc = acc),
                          accout)
        })
        try({
          out <- indmetropwithingibbs(mcmciter, mu, mu, df, lpbest, datlist)
          acc <- mean(out$acc)
          accout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                     nswarm = nswarm, niter = niter,
                                     pso = "AT-PSO-2", mcmc = "IMHwG", acc = acc),
                          accout)
        })
      }
    }
  }
}
save(accout, file = "accout.RData")
}
