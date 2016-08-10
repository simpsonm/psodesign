source("psofun.R")
source("mcmcfun.R")
library(rstan)
load("popdat/popdat.RData")

nbeta <- 1
ndeltasiid <- c(10, 20, 30)
ndeltasfull <- c(5, 7, 9)

psoiter <- 1000
nswarm <- 50
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496

nbhd <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)})

nsubs <- floor(nrow(popdat)/c(4, 2))

mcmciter <- 10000

### acceptance rate sims
{
  accout <- matrix(0, ncol = 7, nrow = 24)
  colnames(accout) <- c("model", "ndelta", "df", "iidIM", "iidIMwG", "fullIM", "fullIMwG")
  idx <- 0
  for(k in 1:3){
    ndeltafull <- ndeltasfull[k]
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
    psopoisfull <- pso(psoiter, nswarm, inertia, cognitive, social, poispsofullinit, nbhd,
                       poislpostfull, datlistfull)
    poisfullmu <- psopoisfull$argmax
    poisfullinit <- psopoisfull$argmax
    poisfulllpbest <- psopoisfull$max
    psolnormfull <- pso(psoiter, nswarm, inertia, cognitive, social, lnormpsofullinit, nbhd,
                        lnormlpostfull, datlistfull)
    lnormfullmu <- psolnormfull$argmax
    lnormfullinit <- psolnormfull$argmax
    lnormfulllpbest <- psolnormfull$max
    ndeltaiid <- ndeltasiid[k]
    poispsoiidinit <- matrix(runif((ndeltaiid + nbeta + 1)*nswarm, -1, 1), ncol = nswarm)
    lnormpsoiidinit <- matrix(runif((ndeltaiid + nbeta + 1 + 1)*nswarm, -1, 1), ncol = nswarm)
    datlistiid <- list(dat = popdat, nbeta = nbeta, ndelta = ndeltaiid,
                    betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                    aphi = 1, bphi = 1)
    psopoisiid <- pso(psoiter, nswarm, inertia, cognitive, social, poispsoiidinit, nbhd,
                      poislpostiid, datlistiid)
    poisiidmu <- psopoisiid$argmax
    poisiidinit <- psopoisiid$argmax
    poisiidlpbest <- psopoisiid$max
    psolnormiid <- pso(psoiter, nswarm, inertia, cognitive, social, lnormpsoiidinit, nbhd,
                       lnormlpostiid, datlistiid)
    lnormiidmu <- psolnormiid$argmax
    lnormiidinit <- psolnormiid$argmax
    lnormiidlpbest <- psolnormiid$max
    for(df in dfs){
      idx <- idx + 1
      print(c(k, df))
      poisfullind <- indmetrop(mcmciter, poislpostfull, poislposthessfull,
                               poisfullmu, poisfullmu, df,
                               poisfulllpbest, datlistfull)
      poisiidind <- indmetrop(mcmciter, poislpostiid, poislposthessiid,
                              poisiidmu, poisiidmu, df,
                              poisiidlpbest, datlistiid)
      poisiidindgibbsout <- poisiidindgibbs(mcmciter, poisiidmu, poisiidmu, df,
                                            poisiidlpbest, datlistiid)
      poisfullindgibbsout <- poisfullindgibbs(mcmciter, poisfullmu, poisfullmu, df,
                                              poisfulllpbest, datlistfull)
      accout[idx,] <- c("Poisson", k, df,
                        mean(poisiidind$acc[-c(1:1000)]),
                        mean(poisiidindgibbsout$acc[-c(1:1000)]),
                        mean(poisfullind$acc[-c(1:1000)]),
                        mean(poisfullindgibbsout$acc[-c(1:1000)]))
      lnormfullind <- indmetrop(mcmciter, lnormlpostfull, lnormlposthessfull,
                               lnormfullmu, lnormfullmu, df,
                               lnormfulllpbest, datlistfull)
      lnormiidind <- indmetrop(mcmciter, lnormlpostiid, lnormlposthessiid,
                              lnormiidmu, lnormiidmu, df,
                              lnormiidlpbest, datlistiid)
      lnormiidindgibbsout <- lnormiidindgibbs(mcmciter, lnormiidmu, lnormiidmu, df,
                                            lnormiidlpbest, datlistiid)
      lnormfullindgibbsout <- lnormfullindgibbs(mcmciter, lnormfullmu, lnormfullmu, df,
                                                lnormfulllpbest, datlistfull)
      idx <- idx + 1
      accout[idx,] <- c("lognormal", k, df,
                        mean(lnormiidind$acc[-c(1:1000)]),
                        mean(lnormiidindgibbsout$acc[-c(1:1000)]),
                        mean(lnormfullind$acc[-c(1:1000)]),
                        mean(lnormfullindgibbsout$acc[-c(1:1000)]))
    }
  }
  save(accout, file = "accout.RData")
}

{
### poisson sims
  mcmcout <- data.frame(ndelta = NULL, model = NULL, ranef = NULL, alg = NULL, df = NULL,
                        nsub = NULL,
                        neff = NULL, dattime = NULL, inittime = NULL, itertime = NULL,
                        accrate1 = NULL, accrate2 = NULL, lpbest = NULL)
  mcmciter <- 50000
  for(k in 1:3){
    print(paste("k = ", k, sep = ""))
    ndeltafull <- ndeltasfull[k]
    ominvscale <- diag(ndeltafull)
    nell <- ndeltafull*(ndeltafull + 1)/2
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
    fullinit <- rep(0, nbeta + ndeltafull + nell)
    ndeltaiid <- ndeltasiid[k]
    datlistiid <- list(dat = popdat, nbeta = nbeta, ndelta = ndeltaiid,
                       betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                       aphi = 1, bphi = 1)
    iidinit <- rep(0, nbeta + ndeltaiid + 1)
    print("poisfullrwgibbs")    
    inittime <- system.time({
      out <- poisfullrwgibbs(5000, fullinit, datlistfull)
      logrwsds <- out$logrwsds
      rwinit <- out$draws[5000 + 1, ]
    })
    itertime <- system.time({
      out <- poisfullrwgibbs(mcmciter, rwinit, datlistfull, tune = FALSE, logrwsds = logrwsds)
    })
    mcmcout <- rbind(data.frame(ndelta = ndeltasfull[k], model = "Poisson", ranef = "full",
                                alg = "RWwG", df = 0,
                                nsub = nrow(popdat),
                                neff = min(effectiveSize(mcmc(out$draws))),
                                dattime = 0,
                                inittime = inittime[3],
                                itertime = itertime[3],
                                accrate1 = min(apply(out$accs, 2, mean)),
                                accrate2 = max(apply(out$accs, 2, mean)),
                                lpbest = -Inf), mcmcout)
    print("poisfullblockrwgibbs")
    inittime <- system.time({
      out <- poisfullrwgibbs(5000, fullinit, datlistfull)
      logrwsds <- out$logrwsds
      rwinit <- out$draws[5000 + 1, ]
      out <- poisfullrwgibbs(1000, rwinit, datlistfull, tune = FALSE, logrwsds = logrwsds)
      sighatfull <- cov(out$draws[,1:(nbeta + ndeltasfull[k])])
      cholhatfull <- t(chol(sighatfull))
      blockinit <- out$draws[1000 + 1,]
      out <- poisfullblockrwgibbs(1000, blockinit, datlistfull, cholhatfull)
      blockinit <- out$draws[1000 + 1,]
      logrwsd <- out$logrwsd
    })
    itertime <- system.time({
      out <- poisfullblockrwgibbs(mcmciter, blockinit, datlistfull, cholhatfull,
                                  tune = FALSE, logrwsd = logrwsd)
    })
    mcmcout <- rbind(data.frame(ndelta = ndeltasfull[k], model = "Poisson", ranef = "full",
                                alg = "blockRWwG",
                                df = 0,
                                nsub = nrow(popdat),
                                neff = min(effectiveSize(mcmc(out$draws))),
                                dattime = 0,
                                inittime = inittime[3],
                                itertime = itertime[3],
                                accrate1 = mean(out$acc),
                                accrate2 = 0,
                                lpbest = -Inf), mcmcout)
    print("poisiidrwgibbs")
    inittime <- system.time({
      out <- poisiidrwgibbs(5000, iidinit, datlistiid)
      logrwsds <- out$logrwsds
      rwinit <- out$draws[5000 + 1, ]
    })
    itertime <- system.time({
      out <- poisiidrwgibbs(mcmciter, iidinit, datlistiid, tune = FALSE, logrwsds = logrwsds)
    })
    mcmcout <- rbind(data.frame(ndelta = ndeltasiid[k], model = "Poisson", ranef = "iid",
                                alg = "RWwG", df = 0,
                                nsub = nrow(popdat),
                                neff = min(effectiveSize(mcmc(out$draws))),
                                dattime = 0,
                                inittime = inittime[3],
                                itertime = itertime[3],
                                accrate1 = min(apply(out$accs, 2, mean)),
                                accrate2 = max(apply(out$accs, 2, mean)),
                                lpbest = -Inf), mcmcout)
    print("poisiidblockrwgibbs")
    inittime <- system.time({
      out <- poisiidrwgibbs(5000, iidinit, datlistiid)
      logrwsds <- out$logrwsds
      rwinit <- out$draws[5000 + 1, ]
      out <- poisiidrwgibbs(1000, rwinit, datlistiid, tune = FALSE, logrwsds = logrwsds)
      sighatiid <- cov(out$draws[,1:(nbeta + ndeltasiid[k])])
      cholhatiid <- t(chol(sighatiid))
      blockinit <- out$draws[1000 + 1,]
      out <- poisiidblockrwgibbs(1000, blockinit, datlistiid, cholhatiid)
      blockinit <- out$draws[1000 + 1,]
      logrwsd <- out$logrwsd
    })
    itertime <- system.time({
      out <- poisiidblockrwgibbs(mcmciter, blockinit, datlistiid, cholhatiid,
                                 tune = FALSE, logrwsd = logrwsd)
    })
    mcmcout <- rbind(data.frame(ndelta = ndeltasiid[k], model = "Poisson", ranef = "iid",
                                alg = "blockRWwG", df = 0,
                                nsub = nrow(popdat),
                                neff = min(effectiveSize(mcmc(out$draws))),
                                dattime = 0,
                                inittime = inittime[3],
                                itertime = itertime[3],
                                accrate1 = mean(out$acc),
                                accrate2 = 0,
                                lpbest = -Inf), mcmcout)
    save(mcmcout, file = "mcmcout.RData")
### lognormal sims
    print("lognormal sims")
    print(paste("k = ", k, sep = ""))
    fullinit <- rep(0, nbeta + ndeltafull + nell + 1)
    iidinit <- rep(0, nbeta + ndeltaiid + 1 + 1)
    print("lnormfullrwgibbs")    
    inittime <- system.time({
      out <- lnormfullrwgibbs(5000, fullinit, datlistfull)
      logrwsds <- out$logrwsds
      rwinit <- out$draws[5000 + 1, ]
    })
    itertime <- system.time({
      out <- lnormfullrwgibbs(mcmciter, rwinit, datlistfull, tune = FALSE, logrwsds = logrwsds)
    })
    mcmcout <- rbind(data.frame(ndelta = ndeltasfull[k], model = "lognormal", ranef = "full",
                                alg = "RWwG", df = 0,
                                nsub = nrow(popdat),
                                neff = min(effectiveSize(mcmc(out$draws))),
                                dattime = 0,
                                inittime = inittime[3],
                                itertime = itertime[3],
                                accrate1 = min(apply(out$accs, 2, mean)),
                                accrate2 = max(apply(out$accs, 2, mean)),
                                lpbest = -Inf), mcmcout)
    print("lnormfullblockrwgibbs")
    inittime <- system.time({
      out <- lnormfullrwgibbs(5000, fullinit, datlistfull)
      logrwsds <- out$logrwsds
      rwinit <- out$draws[5000 + 1, ]
      out <- lnormfullrwgibbs(1000, rwinit, datlistfull, tune = FALSE, logrwsds = logrwsds)
      sighatfull <- cov(out$draws[,1:(nbeta + ndeltasfull[k])])
      cholhatfull <- t(chol(sighatfull))
      blockinit <- out$draws[1000 + 1,]
      out <- lnormfullblockrwgibbs(1000, blockinit, datlistfull, cholhatfull)
      blockinit <- out$draws[1000 + 1,]
      logrwsd <- out$logrwsd
    })
    itertime <- system.time({
      out <- lnormfullblockrwgibbs(mcmciter, blockinit, datlistfull, cholhatfull,
                                  tune = FALSE, logrwsd = logrwsd)
    })
    mcmcout <- rbind(data.frame(ndelta = ndeltasfull[k], model = "lognormal", ranef = "full",
                                alg = "blockRWwG",
                                df = 0,
                                nsub = nrow(popdat),
                                neff = min(effectiveSize(mcmc(out$draws))),
                                dattime = 0,
                                inittime = inittime[3],
                                itertime = itertime[3],
                                accrate1 = mean(out$acc),
                                accrate2 = 0,
                                lpbest = -Inf), mcmcout)
    print("lnormiidrwgibbs")
    inittime <- system.time({
      out <- lnormiidrwgibbs(5000, iidinit, datlistiid)
      logrwsds <- out$logrwsds
      rwinit <- out$draws[5000 + 1, ]
    })
    itertime <- system.time({
      out <- lnormiidrwgibbs(mcmciter, iidinit, datlistiid, tune = FALSE, logrwsds = logrwsds)
    })
    mcmcout <- rbind(data.frame(ndelta = ndeltasiid[k], model = "lognormal", ranef = "iid",
                                alg = "RWwG", df = 0,
                                nsub = nrow(popdat),
                                neff = min(effectiveSize(mcmc(out$draws))),
                                dattime = 0,
                                inittime = inittime[3],
                                itertime = itertime[3],
                                accrate1 = min(apply(out$accs, 2, mean)),
                                accrate2 = max(apply(out$accs, 2, mean)),
                                lpbest = -Inf), mcmcout)
    print("lnormiidblockrwgibbs")
    inittime <- system.time({
      out <- lnormiidrwgibbs(5000, iidinit, datlistiid)
      logrwsds <- out$logrwsds
      rwinit <- out$draws[5000 + 1, ]
      out <- lnormiidrwgibbs(1000, rwinit, datlistiid, tune = FALSE, logrwsds = logrwsds)
      sighatiid <- cov(out$draws[,1:(nbeta + ndeltasiid[k])])
      cholhatiid <- t(chol(sighatiid))
      blockinit <- out$draws[1000 + 1,]
      out <- lnormiidblockrwgibbs(1000, blockinit, datlistiid, cholhatiid)
      blockinit <- out$draws[1000 + 1,]
      logrwsd <- out$logrwsd
    })
    itertime <- system.time({
      out <- lnormiidblockrwgibbs(mcmciter, blockinit, datlistiid, cholhatiid,
                                 tune = FALSE, logrwsd = logrwsd)
    })
    mcmcout <- rbind(data.frame(ndelta = ndeltasiid[k], model = "lognormal", ranef = "iid",
                                alg = "blockRWwG", df = 0,
                                nsub = nrow(popdat),
                                neff = min(effectiveSize(mcmc(out$draws))),
                                dattime = 0,
                                inittime = inittime[3],
                                itertime = itertime[3],
                                accrate1 = mean(out$acc),
                                accrate2 = 0,
                                lpbest = -Inf), mcmcout)
    print("lnormiidgibbs")
    inittime <- system.time({
      out <- lnormiidgibbs(1000, rep(0, ndeltaiid + nbeta + 2), datlistiid)
      blockinit <- out$draws[1000 + 1,]
    })
    itertime <- system.time({
      out <- lnormiidgibbs(mcmciter, blockinit, datlistiid)
    })
    mcmcout <- rbind(data.frame(ndelta = ndeltasiid[k], model = "lognormal", ranef = "iid",
                                alg = "Gibbs", df = 0,
                                nsub = nrow(popdat),
                                neff = min(effectiveSize(mcmc(out$draws))),
                                dattime = 0,
                                inittime = inittime[3],
                                itertime = itertime[3],
                                accrate1 = 0,
                                accrate2 = 0,
                                lpbest = -Inf), mcmcout)
    print("lnormfullgibbs")
    inittime <- system.time({
      out <- lnormfullgibbs(1000, rep(0, ndeltafull + nbeta + nell + 1), datlistfull)
      blockinit <- out$draws[1000 + 1,]
    })
    itertime <- system.time({
      out <- lnormfullgibbs(mcmciter, blockinit, datlistfull)
    })
    mcmcout <- rbind(data.frame(ndelta = ndeltasfull[k], model = "lognormal", ranef = "full",
                                alg = "Gibbs", df = 0,
                                nsub = nrow(popdat),
                                neff = min(effectiveSize(mcmc(out$draws))),
                                dattime = 0,
                                inittime = inittime[3],
                                itertime = itertime[3],
                                accrate1 = 0,
                                accrate2 = 0,
                                lpbest = -Inf), mcmcout)
    save(mcmcout, file = "mcmcout.RData")
    save(mcmcout, file = "mcmcout.RData")
  }
  save(mcmcout, file = "mcmcout.RData")
}

## IMH & IMHwG sims
psoiter <- 10000
nswarm <- 100
nbhd <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)})
{
### poisson sims
  mcmcout <- data.frame(ndelta = NULL, model = NULL, ranef = NULL, alg = NULL, df = NULL,
                        nsub = NULL,
                        neff = NULL, dattime = NULL, inittime = NULL, itertime = NULL,
                        accrate1 = NULL, accrate2 = NULL, lpbest = NULL)
  mcmciter <- 50000
  df <- 100
  for(k in 1:3){
    print(paste("k = ", k, sep = ""))
    ndeltafull <- ndeltasfull[k]
    ominvscale <- diag(ndeltafull)
    nell <- ndeltafull*(ndeltafull + 1)/2
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
    fullinit <- rep(0, nbeta + ndeltafull + nell)
    ndeltaiid <- ndeltasiid[k]
    datlistiid <- list(dat = popdat, nbeta = nbeta, ndelta = ndeltaiid,
                       betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                       aphi = 1, bphi = 1)
    iidinit <- rep(0, nbeta + ndeltaiid + 1)
    print(paste("k = ", k, ", df = ", df, sep = ""))
    print("poisson sims")
    psoinit <- matrix(runif((ndeltafull + nbeta + ndeltafull*(ndeltafull + 1)/2)*nswarm, -1, 1),
                      ncol = nswarm)
    inittime <- system.time({
      psoout <- pso(psoiter, nswarm, inertia, cognitive, social, psoinit, nbhd,
                    poislpostfull, datlistfull)
      mu <- psoout$argmax
      lpbest <- psoout$max
    })
    print("poisfullind")
    itertime <- system.time({
      out <- indmetrop(mcmciter, poislpostfull, poislposthessfull,
                       mu, mu, df,lpbest, datlistfull)
    })
    ## in the full model using the MIH algorithm the signs of the elements of L are unidentified
    ## so to make effective sample size computations comparable: sign the elements of L
    out$draws2 <- out$draws
    L <- matrix(0, ndeltafull, ndeltafull)
    for(ii in 1:ncol(out$draws)){
      ell <- out$draws[ii, nbeta + ndeltafull + 1:nell]
      L[lower.tri(L, TRUE)] <- ell
      L <- t(chol(tcrossprod(L)))
      out$draws2[ii, nbeta + ndeltafull + 1:nell] <- L[lower.tri(L, TRUE)]
    }
    mcmcout <- rbind(data.frame(ndelta = ndeltasfull[k], model = "Poisson", ranef = "full",
                                alg = "IMH", df = df,
                                nsub = nrow(popdat),
                                neff = min(effectiveSize(mcmc(out$draws2))),
                                dattime = 0,
                                inittime = inittime[3],
                                itertime = itertime[3],
                                accrate1 = mean(out$acc),
                                accrate2 = 0,
                                lpbest = max(out$lpbests)), mcmcout)
    print("poisfullindgibbs")
    itertime <- system.time({
      out <- poisfullindgibbs(mcmciter, mu, mu, df, lpbest, datlistfull)
    })
    mcmcout <- rbind(data.frame(ndelta = ndeltasfull[k], model = "Poisson", ranef = "full",
                                alg = "IMHwGibbs",
                                df = df,
                                nsub = nrow(popdat),
                                neff = min(effectiveSize(mcmc(out$draws))),
                                dattime = 0,
                                inittime = inittime[3],
                                itertime = itertime[3],
                                accrate1 = mean(out$acc),
                                accrate2 = 0,
                                lpbest = max(out$lpbests)), mcmcout)
    ndeltaiid <- ndeltasiid[k]
    psoinit <- matrix(runif((ndeltaiid + nbeta + 1)*nswarm, -1, 1), ncol = nswarm)
    print("poisson iid sims")
    inittime <- system.time({
      psoout <- pso(psoiter, nswarm, inertia, cognitive, social, psoinit, nbhd,
                    poislpostiid, datlistiid)
      mu <- psoout$argmax
      lpbest <- psoout$max
    })
    print("poisiidind")
    itertime <- system.time({
      out <- indmetrop(mcmciter, poislpostiid, poislposthessiid,
                       mu, mu, df, lpbest, datlistiid)
    })
    mcmcout <- rbind(data.frame(ndelta = ndeltasiid[k], model = "Poisson", ranef = "iid",
                                alg = "IMH", df = df,
                                nsub = nrow(popdat),
                                neff = min(effectiveSize(mcmc(out$draws))),
                                dattime = 0,
                                inittime = inittime[3],
                                itertime = itertime[3],
                                accrate1 = mean(out$acc),
                                accrate2 = 0,
                                lpbest = max(out$lpbests)), mcmcout)
    print("poisiidindgibbs")
    itertime <- system.time({
      out <- poisiidindgibbs(mcmciter, mu, mu, df, lpbest, datlistiid)
    })
    mcmcout <- rbind(data.frame(ndelta = ndeltasiid[k], model = "Poisson", ranef = "iid",
                                alg = "IMHwGibbs",
                                df = df,
                                nsub = nrow(popdat),
                                neff = min(effectiveSize(mcmc(out$draws))),
                                dattime = 0,
                                inittime = inittime[3],
                                itertime = itertime[3],
                                accrate1 = mean(out$acc),
                                accrate2 = 0,
                                lpbest = max(out$lpbests)), mcmcout)
    print("lognormal sims")
    ndeltafull <- ndeltasfull[k]
    psoinit <- matrix(runif((ndeltafull + nbeta + ndeltafull*(ndeltafull + 1)/2 + 1)*nswarm,
                            -1, 1), ncol = nswarm)
    inittime <- system.time({
      psoout <- pso(psoiter, nswarm, inertia, cognitive, social, psoinit, nbhd,
                    lnormlpostfull, datlistfull)
      mu <- psoout$argmax
      lpbest <- psoout$max
    })
    print("lnormfullind")
    itertime <- system.time({
      out <- indmetrop(mcmciter, lnormlpostfull, lnormlposthessfull,
                       mu, mu, df,lpbest, datlistfull)
    })
    ## in the full model using the MIH algorithm the signs of the elements of L are unidentified
    ## to make effective sample size computations comparable: sign the elements of L
    out$draws2 <- out$draws
    L <- matrix(0, ndeltafull, ndeltafull)
    for(ii in 1:ncol(out$draws)){
      ell <- out$draws[ii, nbeta + ndeltafull + 1:nell]
      L[lower.tri(L, TRUE)] <- ell
      L <- t(chol(tcrossprod(L)))
      out$draws2[ii, nbeta + ndeltafull + 1:nell] <- L[lower.tri(L, TRUE)]
    }
    mcmcout <- rbind(data.frame(ndelta = ndeltasfull[k], model = "lognormal", ranef = "full",
                                alg = "IMH", df = df,
                                nsub = nrow(popdat),
                                neff = min(effectiveSize(mcmc(out$draws2))),
                                dattime = 0,
                                inittime = inittime[3],
                                itertime = itertime[3],
                                accrate1 = mean(out$acc),
                                accrate2 = 0,
                                lpbest = max(out$lpbests)), mcmcout)
    print("lnormfullindgibbs")
    itertime <- system.time({
      out <- lnormfullindgibbs(mcmciter, mu, mu, df, lpbest, datlistfull)
    })
    mcmcout <- rbind(data.frame(ndelta = ndeltasfull[k], model = "lognormal", ranef = "full",
                                alg = "IMHwGibbs",
                                df = df,
                                nsub = nrow(popdat),
                                neff = min(effectiveSize(mcmc(out$draws))),
                                dattime = 0,
                                inittime = inittime[3],
                                itertime = itertime[3],
                                accrate1 = mean(out$acc),
                                accrate2 = 0,
                                lpbest = max(out$lpbests)), mcmcout)
    ndeltaiid <- ndeltasiid[k]
    psoinit <- matrix(runif((ndeltaiid + nbeta + 1 + 1)*nswarm, -1, 1), ncol = nswarm)
    print("lognormal iid sims")
    inittime <- system.time({
      psoout <- pso(psoiter, nswarm, inertia, cognitive, social, psoinit, nbhd,
                    lnormlpostiid, datlistiid)
      mu <- psoout$argmax
      lpbest <- psoout$max
    })
    print("lnormiidind")
    itertime <- system.time({
      out <- indmetrop(mcmciter, lnormlpostiid, lnormlposthessiid,
                       mu, mu, df, lpbest, datlistiid)
    })
    mcmcout <- rbind(data.frame(ndelta = ndeltasiid[k], model = "lognormal", ranef = "iid",
                                alg = "IMH", df = df,
                                nsub = nrow(popdat),
                                neff = min(effectiveSize(mcmc(out$draws))),
                                dattime = 0,
                                inittime = inittime[3],
                                itertime = itertime[3],
                                accrate1 = mean(out$acc),
                                accrate2 = 0,
                                lpbest = max(out$lpbests)), mcmcout)
    print("lnormiidindgibbs")
    itertime <- system.time({
      out <- lnormiidindgibbs(mcmciter, mu, mu, df, lpbest, datlistiid)
    })
    mcmcout <- rbind(data.frame(ndelta = ndeltasiid[k], model = "lognormal", ranef = "iid",
                                alg = "IMHwGibbs",
                                df = df,
                                nsub = nrow(popdat),
                                neff = min(effectiveSize(mcmc(out$draws))),
                                dattime = 0,
                                inittime = inittime[3],
                                itertime = itertime[3],
                                accrate1 = mean(out$acc),
                                accrate2 = 0,
                                lpbest = max(out$lpbests)), mcmcout)
  }
  mcmcoutimh <- mcmcout
  save(mcmcoutimh, file = "mcmcoutimh.RData")
}



## stan sims
mcmcout <- NULL
nbeta <- 1
ndeltasiid <- c(10, 20, 30)
ndeltasfull <- c(5, 7, 9)
load("popdat/popdat.RData")
z <- popdat[,1]
xmat <- as.matrix(popdat[,2])
smat <- popdat[,-c(1:2)]
aphi <- 1
bphi <- 1
asig <- 1
bsig <- 1
mubeta <- array(rep(0, ncol(xmat)), dim = ncol(xmat))
sigbeta <- 10

mcmciter <- 50000
warmup <- 10000

for(k in 1:3){
  print(paste("k = ", k, sep = ""))
  ndeltafull <- ndeltasfull[k]
  ndeltaiid <- ndeltasiid[k]
  d <- ndeltafull + 1
  Einv <- diag(ndeltafull)
  fullstandat <- list(nobs = length(z), ndelta = ndeltafull, nell = ndeltafull*(ndeltafull + 1)/2,
                      nbeta = ncol(xmat), zobs = z,
                      lzobs = log(z), S = smat[,1:ndeltafull], X = xmat,
                      a_phi = aphi, b_phi = bphi, a_sig = asig, b_sig = bsig,
                      mu_beta = mubeta, sig_beta = sigbeta,
                      d_omega = d, Einv_omega = Einv)
  iidstandat <- list(nobs = length(z), ndelta = ndeltaiid, nell = ndeltaiid*(ndeltaiid + 1)/2,
                     nbeta = ncol(xmat), zobs = z,
                     lzobs = log(z), S = smat[,1:ndeltaiid], X = xmat,
                     a_phi = aphi, b_phi = bphi, a_sig = asig, b_sig = bsig,
                     mu_beta = mubeta, sig_beta = sigbeta,
                     d_omega = d, Einv_omega = Einv)
  print("poisson full init")
  inittime <- system.time({out <- stan(file = "poppoisfull.stan", data = fullstandat, chains = 1,
                                       iter = warmup, warmup = warmup,
                                       control = list(adapt_delta = 0.9, max_treedepth = 11))
  })
  print("poisson full fit")
  itertime <- system.time({out <- stan(fit = out, data = fullstandat, chains = 1,
                                       iter = mcmciter, warmup = warmup,
                                       control = list(adapt_delta = 0.9, max_treedepth = 11))
  })
  print("poisson full init time")
  inittime2 <- system.time({out2 <- stan(fit = out, data = fullstandat, chains = 1,
                                        iter = warmup, warmup = warmup,
                                        control = list(adapt_delta = 0.9, max_treedepth = 11))
  })
  mcmcout <- rbind(data.frame(ndelta = ndeltafull, model = "Poisson", ranef = "full",
                              alg = "Stan", df = 0,
                              nsub = nrow(popdat),
                              neff = min(summary(out, pars = c("beta", "delta", "ell"))$
                                         summary[,9]),
                              dattime = 0,
                              inittime = inittime[3],
                              itertime = itertime[3] - inittime2[3],
                              accrate1 = mean(get_sampler_params(out)[[1]][,1]),
                              accrate2 = length(unique(extract(out)$beta))/mcmciter,
                              lpbest = max(extract(out)$lp__)), mcmcout)
  print("poisson iid init")
  inittime <- system.time({out <- stan(file = "poppoisiid.stan", data = iidstandat, chains = 1,
                                       iter = warmup, warmup = warmup)
  })
  print("poisson iid fit")
  itertime <- system.time({out <- stan(fit = out, data = fullstandat, chains = 1,
                                       iter = mcmciter, warmup = warmup)
  })
  inittime2 <- system.time({out2 <- stan(fit = out, data = iidstandat, chains = 1,
                                        iter = warmup, warmup = warmup)
  })
  print("poisson iid init time")
  mcmcout <- rbind(data.frame(ndelta = ndeltaiid, model = "Poisson", ranef = "iid",
                              alg = "Stan", df = 0,
                              nsub = nrow(popdat),
                              neff = min(summary(out, pars = c("beta", "delta", "sig2"))$
                                         summary[,9]),
                              dattime = 0,
                              inittime = inittime[3],
                              itertime = itertime[3] - inittime2[3],
                              accrate1 = mean(get_sampler_params(out)[[1]][,1]),
                              accrate2 = length(unique(extract(out)$beta))/mcmciter,
                              lpbest = max(extract(out)$lp__)), mcmcout)
  print("lognormal full init")
  inittime <- system.time({out <- stan(file = "poplognormfull.stan", data = fullstandat,
                                       chains = 1, iter = warmup, warmup = warmup,
                                       control = list(adapt_delta = 0.9, max_treedepth = 11))
  })
  print("lognormal full init fit")
  itertime <- system.time({out <- stan(fit = out, data = fullstandat, chains = 1,
                                       iter = mcmciter, warmup = warmup,
                                       control = list(adapt_delta = 0.9, max_treedepth = 11))
  })
  print("lognormal full init init time")
  inittime2 <- system.time({out2 <- stan(fit = out, data = fullstandat,
                                         chains = 1, iter = warmup, warmup = warmup,
                                         control = list(adapt_delta = 0.9, max_treedepth = 11))
  })
  mcmcout <- rbind(data.frame(ndelta = ndeltafull, model = "lognormal", ranef = "full",
                              alg = "Stan", df = 0,
                              nsub = nrow(popdat),
                              neff = min(summary(out, pars = c("beta", "delta", "ell", "phi2"))$
                                         summary[,9]),
                              dattime = 0,
                              inittime = inittime[3],
                              itertime = itertime[3] - inittime2[3],
                              accrate1 = mean(get_sampler_params(out)[[1]][,1]),
                              accrate2 = length(unique(extract(out)$beta))/mcmciter,
                              lpbest = max(extract(out)$lp__)), mcmcout)
  print(" lognormal iid init")
  inittime <- system.time({out <- stan(file = "poplognormiid.stan", data = fullstandat, chains = 1,
                                       iter = warmup, warmup = warmup)
  })
  print(" lognormal iid fit")
  itertime <- system.time({out <- stan(fit = out, data = fullstandat, chains = 1,
                                       iter = mcmciter, warmup = warmup)
  })
  print(" lognormal iid init time")
  inittime2 <- system.time({out2 <- stan(fit = out, data = fullstandat, chains = 1,
                                         iter = warmup, warmup = warmup)
  })
  mcmcout <- rbind(data.frame(ndelta = ndeltaiid, model = "lognormal", ranef = "iid",
                              alg = "Stan", df = 0,
                              nsub = nrow(popdat),
                              neff = min(summary(out, pars = c("beta", "delta", "sig2", "phi2"))$
                                         summary[,9]),
                              dattime = 0,
                              inittime = inittime[3],
                              itertime = itertime[3] - inittime2[3],
                              accrate1 = mean(get_sampler_params(out)[[1]][,1]),
                              accrate2 = length(unique(extract(out)$beta))/mcmciter,
                              lpbest = max(extract(out)$lp__)), mcmcout)
}
mcmcoutstan <- mcmcout
save(mcmcoutstan, file = "mcmcoutstan.RData")
























## stochastic pso MCMC
## basic result:
## the MCMC gets you closer than stochastic PSO got you...
## ...but often not close enough
## ...and numerical problems happen with the fully parameterized version
mcmcout2 <- data.frame(ndelta = NULL, model = NULL, alg = NULL, df = NULL, nsub = NULL,
                      neff = NULL, dattime = NULL, inittime = NULL, itertime = NULL,
                      accrate1 = NULL, accrate2 = NULL, lpbest = NULL)
mcmciter <- 200000

for(k in 1:3){
  ndelta <- ndeltasfull[k]
  ominvscale <- diag(ndelta)
  dattime <- system.time({
    nell <- ndelta*(ndelta + 1)/2
    K <- commutation.matrix(ndelta)
    M <- elimination.matrix(ndelta)
    N2 <- 2*N.matrix(ndelta)
    R <- matrix(0, nell, nell)
    MKprime <- tcrossprod(M, K)
    diag(R)[(ndelta+1)*(1:ndelta) - (1:ndelta)*(1 + 1:ndelta)/2 + 1:ndelta - ndelta] <-
      (ndelta + 2 - 1:ndelta)
    dvelldvell <- - M%*%tcrossprod(kronecker(diag(ndelta), ominvscale), M)
  })
  datlistfull <- list(dat = popdat, nbeta = nbeta, ndelta = ndelta,
                      betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                      aphi = 1, bphi = 1, omdf = ndelta + 1, ominvscale = diag(ndelta),
                      K = K, M = M, N2 = N2, R = R, MKprime = MKprime, dvelldvell = dvelldvell)
  ndelta <- ndeltasiid[k]
  datlistiid <- list(dat = popdat, nbeta = nbeta, ndelta = ndelta,
                      betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                      aphi = 1, bphi = 1)
  for(df in dfs){
    for(nsub in nsubs){
      print(c(k, df, nsub))
      print("full pso")
      ndelta <- ndeltasfull[k]
      psoinit <- matrix(runif((ndelta + nbeta + ndelta*(ndelta + 1)/2)*nswarm, -1, 1),
                        ncol = nswarm)
      inittime <- system.time({
        psoout <- stochpso(psoiter, nswarm, nsub, inertia, cognitive, social,
                           psoinit, nbhd, poislpostfull, datlistfull)
        mu <- psoout$argmax
        lpbest <- psoout$max
      })
      ## itertime <- system.time({
      ##   out <- indmetrop(mcmciter, poislpostfull, poislposthessfull,
      ##                    mu, mu, df, lpbest, datlist)
      ## })
      ## mcmcout2 <- rbind(data.frame(ndelta = ndelta, model = "Poisson", ranef = "full",
      ##                              alg = "IMH", df = df,
      ##                              nsub = nsub,
      ##                              neff = min(effectiveSize(mcmc(out$draws))),
      ##                              dattime = 0,
      ##                              inittime = inittime[3],
      ##                              itertime = itertime[3],
      ##                              accrate1 = mean(out$acc),
      ##                              accrate2 = mean(out$acc[-c(1:100000)]),
      ##                              lpbest = max(out$lpbests)), mcmcout2)
      print("full gibbs")
      itertime <- system.time({
        out <- poisfullindgibbs(mcmciter, mu, mu, df,lpbest, datlistfull)
      })
      mcmcout2 <- rbind(data.frame(ndelta = ndelta, model = "Poisson", ranef = "full",
                                   alg = "IMHwGibbs", df = df,
                                   nsub = nsub,
                                   neff = min(effectiveSize(mcmc(out$draws))),
                                   dattime = 0,
                                   inittime = inittime[3],
                                   itertime = itertime[3],
                                   accrate1 = mean(out$acc),
                                   accrate2 = mean(out$acc[-c(1:100000)]),
                                   lpbest = max(out$lpbests)), mcmcout2)
    }
    print("iid pso")
    ndelta <- ndeltasiid[k]
    psoinit <- matrix(runif((ndelta + nbeta + 1)*nswarm, -1, 1), ncol = nswarm)
    inittime <- system.time({
      psoout <- stochpso(psoiter, nswarm, nsub, inertia, cognitive, social,
                         psoinit, nbhd, poislpostiid, datlistiid)
      mu <- psoout$argmax
      lpbest <- psoout$max
    })
    ## print("iid IMH")
    ## itertime <- system.time({
    ##   out <- indmetrop(mcmciter, poislpostiid, poislposthessiid, mu, mu, df, lpbest, datlist)
    ## })
    ## mcmcout2 <- rbind(data.frame(ndelta = ndelta, model = "Poisson", ranef = "iid",
    ##                              alg = "IMH", df = df,
    ##                              nsub = nsub,
    ##                              neff = min(effectiveSize(mcmc(out$draws))),
    ##                              dattime = 0,
    ##                              inittime = inittime[3],
    ##                              itertime = itertime[3],
    ##                              accrate1 = mean(out$acc),
    ##                              accrate2 = mean(out$acc[-c(1:100000)]),
    ##                              lpbest = max(out$lpbests)), mcmcout2)
    print("iid Gibbs")
    itertime <- system.time({
      out <- poisiidindgibbs(mcmciter, mu, mu, df, lpbest, datlistiid)
    })
    mcmcout2 <- rbind(data.frame(ndelta = ndelta, model = "Poisson", ranef = "iid",
                                 alg = "IMHwGibbs", df = df,
                                 nsub = nsub,
                                 neff = min(effectiveSize(mcmc(out$draws))),
                                 dattime = 0,
                                 inittime = inittime[3],
                                 itertime = itertime[3],
                                 accrate1 = mean(out$acc),
                                 accrate2 = mean(out$acc[-c(1:100000)]),
                                 lpbest = max(out$lpbests)), mcmcout2)
  }
}
save(mcmcout2, file = "mcmcout2.RData")


### IMH w/ joint proposal from prior
### (very poor acceptance rates)
### maybe only do it on the margins... ?
mcmciter <- 50000
mcmcout3 <- data.frame(ndelta = NULL, model = NULL, alg = NULL, df = NULL, nsub = NULL,
                      neff = NULL, dattime = NULL, inittime = NULL, itertime = NULL,
                      accrate1 = NULL, accrate2 = NULL, lpbest = NULL)

{
  for(ndelta in ndeltas){
    ominvscale <- diag(ndelta)
    dattime <- system.time({
      nell <- ndelta*(ndelta + 1)/2
      K <- commutation.matrix(ndelta)
      M <- elimination.matrix(ndelta)
      N2 <- 2*N.matrix(ndelta)
      R <- matrix(0, nell, nell)
      MKprime <- tcrossprod(M, K)
      diag(R)[(ndelta+1)*(1:ndelta) - (1:ndelta)*(1 + 1:ndelta)/2 + 1:ndelta - ndelta] <-
        (ndelta + 2 - 1:ndelta)
      dvelldvell <- - M%*%tcrossprod(kronecker(diag(ndelta), ominvscale), M)
    })
    datlist <- list(dat = popdat, nbeta = nbeta, ndelta = ndelta,
                    betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                    aphi = 1, bphi = 1, omdf = ndelta + 1, ominvscale = diag(ndelta),
                    K = K, M = M, N2 = N2, R = R, MKprime = MKprime, dvelldvell = dvelldvell)
    fullinit <- rep(0, nbeta + ndelta + nell)
    iidinit <- rep(0, nbeta + ndelta + 1)
    itertime <- system.time({
      out <- poisfullindpriorpropgibbs(mcmciter, fullinit, datlistfull)
    })
    mcmcout3 <- rbind(data.frame(ndelta = ndelta, model = "full", alg = "IMHwGibbsProp", df = 0,
                                 nsub = nrow(popdat),
                                 neff = min(effectiveSize(mcmc(out$draws))),
                                 dattime = 0,
                                 inittime = 0,
                                 itertime = itertime[3],
                                 accrate1 = min(apply(out$accs,2,mean)),
                                 accrate2 = max(apply(out$accs,2,mean)),
                                 lpbest = 0), mcmcout3)
    itertime <- system.time({
        out <- poisiidindpriorpropgibbs(mcmciter, iidinit, datlistiid)
    })
    mcmcout3 <- rbind(data.frame(ndelta = ndelta, model = "iid", alg = "IMHwGibbsProp", df = 0,
                                 nsub = nrow(popdat),
                                 neff = min(effectiveSize(mcmc(out$draws))),
                                 dattime = 0,
                                 inittime = 0,
                                 itertime = itertime[3],
                                 accrate1 = min(apply(out$accs,2,mean)),
                                 accrate2 = max(apply(out$accs,2,mean)),
                                 lpbest = 0), mcmcout3)
  }
  save(mcmcout3, file = "mcmcout3.RData")
}



ndelta <- 30

nsub <- floor(nrow(popdat)/2)

df <- 100
ominvscale <- diag(ndelta)
nell <- ndelta*(ndelta + 1)/2
K <- commutation.matrix(ndelta)
M <- elimination.matrix(ndelta)
N2 <- 2*N.matrix(ndelta)
R <- matrix(0, nell, nell)
MKprime <- tcrossprod(M, K)
diag(R)[(ndelta+1)*(1:ndelta) - (1:ndelta)*(1 + 1:ndelta)/2 + 1:ndelta - ndelta] <-
  (ndelta + 2 - 1:ndelta)
dvelldvell <- - M%*%tcrossprod(kronecker(diag(ndelta), ominvscale), M)

datlist <- list(dat = popdat, nbeta = nbeta, ndelta = ndelta,
                betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                aphi = 1, bphi = 1)

psoinitfull <- matrix(runif((ndelta + nbeta + ndelta*(ndelta + 1)/2)*nswarm, -1, 1),
                      ncol = nswarm)
psoinitiid <- matrix(runif((ndelta + nbeta + 1)*nswarm, -1, 1),
                      ncol = nswarm)
psoiter <- 200

psoout <- stochpso(psoiter, nswarm, nsub, inertia, cognitive, social,
                   psoinitiid, nbhd, poislpostiid, datlist)

psoout$realmaxes <- apply(psoout$argmaxes, 2, poislpostiid, datlist)

par(mfrow=c(2,1))
plot(ts(psoout$maxes))
plot(ts(psoout$realmaxes))

mu <- psoout$argmax
lpbest <- poislpostiid(mu, datlist)

mcmciter <- 50000

out1 <- indmetrop(mcmciter, poislpostiid, poislposthessiid,
                  mu, mu, df, lpbest, datlist)

out2 <- poisiidindgibbs(mcmciter, mu, mu, df, lpbest, datlist)


par(mfrow=c(2,1))
plot(ts(out1$lpbests))
plot(ts(out2$lpbests))

c(mean(out1$acc[-c(1:30000)]), mean(out2$acc[-c(1:40000)]))

sum(out1$acc) 


plot(ts(out1$draws[,1]))

nn <- 1000000
idx <- sample(1:nrow(popdat), nn, TRUE)

popdat2 <- popdat[idx,]
popdat2[,1] <- popdat2[,1] + rnorm(nn)



ndelta <- 5
nsub <- floor(nrow(popdat)/10)

df <- 100
datlist <- list(dat = popdat2, nbeta = nbeta, ndelta = ndelta,
                betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                aphi = 1, bphi = 1)

psoinitiid <- matrix(runif((ndelta + nbeta + 1)*nswarm, -1, 1),
                      ncol = nswarm)
psoiter <- 200

psotime <- system.time(psoout <- stochpso(psoiter, nswarm, nsub, inertia, cognitive, social,
                                          psoinitiid, nbhd, poislpostiid, datlist))

mu <- psoout$argmax
lpbest <- poislpostiid(mu, datlist)

mcmciter <- 10000
mcmctime1 <- system.time(out1 <- indmetrop(mcmciter, poislpostiid, poislposthessiid,
                                           mu, mu, df, lpbest, datlist))

mcmctime2 <- system.time(out2 <- poisiidrwgibbs(10000, rep(0, nbeta + ndelta + 1), datlist))
cholhat <- t(chol(cov(out2$draws[-c(1:1000),1:(nbeta + ndelta)])))

mcmctime3 <- system.time(out3 <- poisiidblockrwgibbs(10000, out2$draws[10000+1,], datlist, cholhat))


cbind(mcmctime1, mcmctime2, mcmctime3)

min(effectiveSize(mcmc(out1$draws[-c(1:5000),])))
min(effectiveSize(mcmc(out2$draws[-c(1:5000),])))
min(effectiveSize(mcmc(out3$draws[-c(1:5000),])))

mean(out1$acc)

par(mfrow=c(2,1))
plot(ts(out1$lpbests))
plot(ts(out1$draws[,1]))
