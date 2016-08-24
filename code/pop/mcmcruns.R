source("../psofun.R")
source("popfun.R")
source("mcmcfun.R")
library(rstan)
library(matrixcalc)
load("popdat/popdat.RData")

nbeta <- 1
ndeltasiid <- c(10, 30)
ndeltasfull <- c(5, 10)
ranefs <- c("iid", "full")
models <- c("pois", "lnorm")
ccc <- 0.1
rate <- 0.5
psoiter <- 2000
nswarm <- 50
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496
nbhd <- sapply(1:nswarm, function(x){return( (x + -3:3 - 1)%%nswarm + 1)})
mcmciter <- 50000
warmup <- 10000
df <- 100

datlistiids <- list()
datlistfulls <- list()
standatlistiids <- list()
standatlistfulls <- list()
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
  z <- popdat[,1]
  xmat <- as.matrix(popdat[,2])
  smat <- popdat[,-c(1:2)]
  aphi <- 1
  bphi <- 1
  asig <- 1
  bsig <- 1
  mubeta <- array(rep(0, ncol(xmat)), dim = ncol(xmat))
  sigbeta <- 10
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
  datlistiids[[idelta]] <- datlistiid
  datlistfulls[[idelta]] <- datlistfull
  standatlistiids[[idelta]] <- iidstandat
  standatlistfulls[[idelta]] <- fullstandat
}

mcmcout <- data.frame(ndelta = NULL, model = NULL, ranef = NULL, alg = NULL, df = NULL,
                      neff = NULL, dattime = NULL, inittime = NULL, itertime = NULL,
                      accrate1 = NULL, accrate2 = NULL, lpbest = NULL)

for(idelta in 1:length(ndeltasiid)){
  datlistiid <- datlistiids[[idelta]]
  datlistfull <- datlistfulls[[idelta]]
  standatiid <- standatlistiids[[idelta]]
  standatfull <-   standatlistfulls[[idelta]]
  ndeltafull <- ndeltasfull[idelta]
  ndeltaiid <- ndeltasiid[idelta]
  nell <- ndeltafull*(ndeltafull + 1)/2
  for(ranef in ranefs){
    datlist <- switch(ranef, iid = datlistiid, full = datlistfull)
    standat <- switch(ranef, iid = standatiid, full = standatfull)
    ndelta <- switch(ranef, iid = ndeltaiid, full = ndeltafull)
    for(model in models){
      print(c(model, ranef, idelta))
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
      indmetropwithingibbs <- switch(paste(model, ranef, sep=""),
                                     poisiid = poisiidindgibbs,
                                     poisfull = poisfullindgibbs,
                                     lnormiid = lnormiidindgibbs,
                                     lnormfull = lnormfullindgibbs)
      rwgibbs <- switch(paste(model, ranef, sep=""),
                        poisiid = poisiidrwgibbs,
                        poisfull = poisfullrwgibbs,
                        lnormiid = lnormiidrwgibbs,
                        lnormfull = lnormfullrwgibbs)
      blockrwgibbs <- switch(paste(model, ranef, sep=""),
                        poisiid = poisiidblockrwgibbs,
                        poisfull = poisfullblockrwgibbs,
                        lnormiid = lnormiidblockrwgibbs,
                        lnormfull = lnormfullblockrwgibbs)
      if(model == "pois"){
        stanfile <- paste("pop", model, ranef, ".stan", sep="")
      } else {
        stanfile <- paste("pop", "lognorm", ranef, ".stan", sep="")
      }
      npar <- nbeta + ndelta + (1 - nell)*(ranef == "iid") + nell + 1*(model == "lnorm")
      print("pso")
      inittime <- system.time({
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
        psoinit <- matrix(runif(npar*nswarm, -1, 1), ncol = nswarm) + mubfgs
        psoinit[,1] <- mubfgs
        psoout <- pso(psoiter, nswarm, 0.9, cognitive, social, psoinit, nbhd,
                      lpost, datlist = datlist, tune = TRUE, style = "adaptive",
                      rate = rate, ccc = ccc)
        mu <- psoout$argmax
        lpbest <- psoout$max
      })
      print("IMH")
      itertime <- system.time({
        out <- indmetrop(mcmciter, lpost, lposthess, mu, mu, df,lpbest, datlist)
      })
      ## in the full model using the MIH algorithm the signs of the elements of L are unidentified
      ## so to make effective sample size computations comparable: sign the elements of L
      if(ranef == "full"){
        out$draws2 <- out$draws
        L <- matrix(0, ndelta, ndelta)
        for(ii in 1:ncol(out$draws)){
          ell <- out$draws[ii, nbeta + ndelta + 1:nell]
          L[lower.tri(L, TRUE)] <- ell
          L <- t(chol(tcrossprod(L)))
          out$draws2[ii, nbeta + ndelta + 1:nell] <- L[lower.tri(L, TRUE)]
        }
      }
      mcmcout <- rbind(data.frame(ndelta = ndelta, model = model, ranef = ranef,
                                  alg = "IMH", df = df,
                                  neff = ifelse(ranef == "full",
                                                min(effectiveSize(mcmc(out$draws2))),
                                                min(effectiveSize(mcmc(out$draws)))),
                                  dattime = 0,
                                  inittime = inittime[3],
                                  itertime = itertime[3],
                                  accrate1 = mean(out$acc),
                                  accrate2 = 0,
                                  lpbest = max(out$lpbests)), mcmcout)
      print("IMHwG")
      itertime <- system.time({
        out <- indmetropwithingibbs(mcmciter, mu, mu, df, lpbest, datlist)
      })
      mcmcout <- rbind(data.frame(ndelta = ndelta, model = model, ranef = ranef,
                                  alg = "IMHwGibbs",
                                  df = df,
                                  neff = min(effectiveSize(mcmc(out$draws))),
                                  dattime = 0,
                                  inittime = inittime[3],
                                  itertime = itertime[3],
                                  accrate1 = mean(out$acc),
                                  accrate2 = 0,
                                  lpbest = max(out$lpbests)), mcmcout)
      print("RWwGibbs")    
      inittime <- system.time({
        out <- rwgibbs(warmup, rep(0, npar), datlist)
        logrwsds <- out$logrwsds
        rwinit <- out$draws[warmup + 1, ]
      })
      itertime <- system.time({
        out <- rwgibbs(mcmciter, rwinit, datlist, tune = FALSE, logrwsds = logrwsds)
      })
      mcmcout <- rbind(data.frame(ndelta = ndelta, model = model, ranef = ranef,
                                  alg = "RWwG", df = 0,
                                  neff = min(effectiveSize(mcmc(out$draws))),
                                  dattime = 0,
                                  inittime = inittime[3],
                                  itertime = itertime[3],
                                  accrate1 = min(apply(out$accs, 2, mean)),
                                  accrate2 = max(apply(out$accs, 2, mean)),
                                  lpbest = -Inf), mcmcout)
      print("blockRWwGibbs")
      inittime <- system.time({
        out <- rwgibbs(warmup, rep(0, npar), datlist)
        logrwsds <- out$logrwsds
        rwinit <- out$draws[warmup + 1, ]
        out <- rwgibbs(1000, rwinit, datlist, tune = FALSE, logrwsds = logrwsds)
        sighatfull <- cov(out$draws[,1:(nbeta + ndelta)])
        cholhatfull <- t(chol(sighatfull))
        blockinit <- out$draws[1000 + 1,]
        out <- blockrwgibbs(1000, blockinit, datlist, cholhatfull)
        blockinit <- out$draws[1000 + 1,]
        logrwsd <- out$logrwsd
      })
      itertime <- system.time({
        out <- blockrwgibbs(mcmciter, blockinit, datlist, cholhatfull,
                            tune = FALSE, logrwsd = logrwsd)
      })
      mcmcout <- rbind(data.frame(ndelta = ndelta, model = model, ranef = ranef,
                                  alg = "blockRWwG",
                                  df = 0,
                                  neff = min(effectiveSize(mcmc(out$draws))),
                                  dattime = 0,
                                  inittime = inittime[3],
                                  itertime = itertime[3],
                                  accrate1 = mean(out$acc),
                                  accrate2 = 0,
                                  lpbest = -Inf), mcmcout)
      print("Gibbs")
      if(model == "lnorm"){
        gibbs <- switch(paste(model, ranef, sep=""),
                        lnormiid = lnormiidgibbs,
                        lnormfull = lnormfullgibbs)
        inittime <- system.time({
          out <- gibbs(warmup, rep(0, npar), datlist)
          blockinit <- out$draws[warmup + 1,]
        })
        itertime <- system.time({
          out <- gibbs(mcmciter, blockinit, datlist)
        })
        mcmcout <- rbind(data.frame(ndelta = ndelta, model = model, ranef = ranef,
                                    alg = "Gibbs", df = 0,
                                    neff = min(effectiveSize(mcmc(out$draws))),
                                    dattime = 0,
                                    inittime = inittime[3],
                                    itertime = itertime[3],
                                    accrate1 = 0,
                                    accrate2 = 0,
                                    lpbest = -Inf), mcmcout)
      }
      print("Stan")
      inittime0 <- system.time({
        out <- stan(file = stanfile, data = standat, chains = 0,
                    iter = 0, warmup = 0,
                    control = list(adapt_delta = 0.9, max_treedepth = 11))
      })
      inittime <- system.time({
        out <- stan(fit = out, data = standat, chains = 1,
                    iter = warmup, warmup = warmup,
                    control = list(adapt_delta = 0.9, max_treedepth = 11))
      })
      itertime <- system.time({
        out <- stan(fit = out, data = standat, chains = 1,
                    iter = mcmciter + warmup, warmup = warmup,
                    control = list(adapt_delta = 0.9, max_treedepth = 11))
      })
      parnames <- switch(paste(model, ranef, sep=""),
                         poisiid = c("beta", "delta", "sig2"),
                         poisfull = c("beta", "delta", "ell"),
                         lnormiid = c("beta", "delta", "sig2", "phi2"),
                         lnormfull = c("beta", "delta", "ell", "phi2"))
      mcmcout <- rbind(data.frame(ndelta = ndelta, model = model, ranef = ranef,
                                  alg = "Stan", df = 0,
                                  neff = min(summary(out, pars = parnames)$
                                             summary[,9]),
                                  dattime = 0,
                                  inittime = inittime[3] + inittime0[3],
                                  itertime = itertime[3] - inittime[3],
                                  accrate1 = mean(get_sampler_params(out)[[1]][,1]),
                                  accrate2 = length(unique(extract(out)$beta))/mcmciter,
                                  lpbest = max(extract(out)$lp__)), mcmcout)
      print("saving")
      write(mcmcout, file = "mcmcout.csv", row.names=FALSE)
    }
  }
}
