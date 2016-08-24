source("../psofun.R")
source("popfun.R")
source("mcmcfun.R")
load("popdat/popdat.RData")
library(MCMCpack)


nbeta <- 1
ndeltasiid <- c(10, 30)
ndeltasfull <- c(5, 15)
niter <- 2000
nswarm <- 50
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496
nbhdnames <- c("ring-1", "ring-3", "global")
rates <- c(0.3, 0.5)
dfs <- c(3, 5)
ccc <- c(0.1)
alpha <- .2*niter
beta <- 1
ranefs <- c("iid", "full")
models <- c("pois", "lnorm")
accpsoout <- NULL
nbhd <- list()
nbhd[[1]] <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)}) ## ring-1
nbhd[[2]] <- sapply(1:nswarm, function(x){return( (x + -3:3 - 1)%%nswarm + 1)}) ## ring-3
nbhd[[3]] <- matrix(1:nswarm, ncol=nswarm, nrow=nswarm) ## global
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

## psoid <- 1
## {
##   for(idelta in 1:length(ndeltasiid)){
##     datlistiid <- datlistiids[[idelta]]
##     datlistfull <- datlistfulls[[idelta]]
##     ndeltafull <- ndeltasfull[idelta]
##     ndeltaiid <- ndeltasiid[idelta]
##     nell <- ndeltafull*(ndeltafull + 1)/2
##     for(ranef in ranefs){
##       datlist <- switch(ranef, iid = datlistiids[[idelta]], full = datlistfulls[[idelta]])
##       ndelta <- switch(ranef, iid = ndeltaiid, full = ndeltafull)
##       for(model in models){
##         lpost <- switch(paste(model, ranef, sep=""),
##                         poisiid = poislpostiid,
##                         lnormiid = lnormlpostiid,
##                         poisfull = poislpostfull,
##                         lnormfull = lnormlpostfull)
##         npar <- nbeta + ndelta + (1 - nell)*(ranef == "iid") + nell + 1*(model == "lnorm")
##         init <- matrix(runif(npar*nswarm, -100, 100), ncol = nswarm)
##         for(m in 1:3){
##           cat("idelta = ")
##           cat(idelta)
##           cat(", ")
##           cat(model)
##           cat(" ")
##           cat(ranef)
##           cat(" ")
##           cat(ndelta)
##           cat(", ")
##           cat("nbhd = ")
##           cat(nbhdnames[m])
##           cat("\n")
##           print("PSO")
##           psotemp <- pso(niter, nswarm, inertia, cognitive, social, init, nbhd[[m]],
##                          lpost, datlist = datlist)
##           psotempout <- data.frame(model = model, ranef = ranef, ndelta = ndelta,
##                                    algorithm = "PSO",
##                                    nbhd = nbhdnames[m],
##                                    iteration = rep(0:niter, each = npar),
##                                    maxes = rep(psotemp$maxes, each = npar),
##                                    pars = c(psotemp$argmaxes), psoid = psoid)
##           psoid <- psoid + 1
##           accpsoout <- rbind(accpsoout, psotempout)
##           print("BBPSO-MC")
##           psotemp <- bbpso(niter, nswarm, 0, 1, init, nbhd[[m]], lpost, Inf, FALSE,
##                            0, datlist = datlist)
##           psotempout <- data.frame(model = model, ranef = ranef, ndelta = ndelta,
##                                    algorithm = "BBPSO-MC",
##                                    nbhd = nbhdnames[m],
##                                    iteration = rep(0:niter, each = npar),
##                                    maxes = rep(psotemp$maxes, each = npar),
##                                    pars = c(psotemp$argmaxes), psoid = psoid)
##           psoid <- psoid + 1          
##           accpsoout <- rbind(accpsoout, psotempout)            
##           print("BBPSOxp-MC")
##           psotemp <- bbpso(niter, nswarm, 0, 1, init, nbhd[[m]], lpost, Inf, FALSE,
##                            .5, datlist = datlist)
##           psotempout <- data.frame(model = model, ranef = ranef, ndelta = ndelta,
##                                    algorithm = "BBPSOxp-MC",
##                                    nbhd = nbhdnames[m],
##                                    iteration = rep(0:niter, each = npar),
##                                    maxes = rep(psotemp$maxes, each = npar),
##                                    pars = c(psotemp$argmaxes), psoid = psoid)
##           psoid <- psoid + 1          
##           accpsoout <- rbind(accpsoout, psotempout)
##           print("DI-PSO")
##           psotemp <- pso(niter, nswarm, inertia, social, cognitive, init, nbhd[[m]],
##                          lpost, datlist = datlist, tune = TRUE, style = "deterministic",
##                          alpha = alpha, beta = beta)
##           psotempout <- data.frame(model = model, ranef = ranef, ndelta = ndelta,
##                                    algorithm = "DI-PSO",
##                                    nbhd = nbhdnames[m],
##                                    iteration = rep(0:niter, each = npar),
##                                    maxes = rep(psotemp$maxes, each = npar),
##                                    pars = c(psotemp$argmaxes), psoid = psoid)
##           psoid <- psoid + 1          
##           accpsoout <- rbind(accpsoout, psotempout)
##           print("AT-PSO & AT-BBPSO-MC & AT-BBPSOxp-MC")
##           for(rate in rates){
##             psotemp <- pso(niter, nswarm, 0.9, cognitive, social, init, nbhd[[m]],
##                            lpost, datlist = datlist, tune = TRUE, style = "adaptive",
##                            rate = rate, ccc = ccc)
##             psotempout <- data.frame(model = model, ranef = ranef, ndelta = ndelta,
##                                      algorithm = paste("AT-PSO", rate, ccc, sep="-"),
##                                      nbhd = nbhdnames[m],
##                                      iteration = rep(0:niter, each = npar),
##                                      maxes = rep(psotemp$maxes, each = npar),
##                                      pars = c(psotemp$argmaxes), psoid = psoid)
##             psoid <- psoid + 1            
##             for(df in dfs){
##               print(paste(c(rate, df)))
##               accpsoout <- rbind(accpsoout, psotempout)
##               psotemp <- bbpso(niter, nswarm, 1, rate, init, nbhd[[m]], lpost, df,
##                                TRUE, 0, datlist = datlist, ccc = ccc)
##               psotempout <- data.frame(model = model, ranef = ranef, ndelta = ndelta,
##                                        algorithm =
##                                          paste("AT-BBPSO-MC", df, rate, ccc, sep="-"),
##                                        nbhd = nbhdnames[m],
##                                        iteration = rep(0:niter, each = npar),
##                                        maxes = rep(psotemp$maxes, each = npar),
##                                        pars = c(psotemp$argmaxes), psoid = psoid)
##               psoid <- psoid + 1              
##               accpsoout <- rbind(accpsoout, psotempout)
##               psotemp <- bbpso(niter, nswarm, 1, rate, init, nbhd[[m]], lpost, df,
##                                TRUE, 0.5, datlist = datlist, ccc = ccc)
##               psotempout <- data.frame(model = model, ranef = ranef, ndelta = ndelta,
##                                        algorithm =
##                                          paste("AT-BBPSOxp-MC", df, rate, ccc, sep="-"),
##                                        nbhd = nbhdnames[m],
##                                        iteration = rep(0:niter, each = npar),
##                                        maxes = rep(psotemp$maxes, each = npar),
##                                        pars = c(psotemp$argmaxes), psoid = psoid)
##               psoid <- psoid + 1              
##               accpsoout <- rbind(accpsoout, psotempout)
##             }
##           }
##           write.csv(accpsoout, file = "accpsoout.csv", row.names=FALSE)
##         }
##       }
##     }
##   }
## }

accpsoout <- read.csv("accpsoout.csv")

accout <- NULL
niters <- c(1000, 1500, 2000)
mcmciter <- 10000
df <- 100
for(id in 1:max(accpsoout$psoid)){
  cat("\n PSOID = ")
  cat(id)
  cat(", niter = ")
  psooutid <- subset(accpsoout, psoid == id)
  model <- psooutid$model[1]
  ranef <- psooutid$ranef[1]
  nbhd <- psooutid$nbhd[1]
  alg <- psooutid$algorithm[1]
  ndelta <- psooutid$ndelta[1]
  idelta <- ifelse(ndelta %in% c(5, 10), 1, 2)
  nell <- ndelta*(ndelta + 1)/2
  datlist <- switch(ranef, iid = datlistiids[[idelta]], full = datlistfulls[[idelta]])
  lpost <- switch(paste(model, ranef, sep=""),
                  poisiid = poislpostiid,
                  lnormiid = lnormlpostiid,
                  poisfull = poislpostfull,
                  lnormfull = lnormlpostfull)
  npar <- nbeta + ndelta + (1 - nell)*(ranef == "iid") + nell + 1*(model == "lnorm")
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
  for(niter in niters){
    cat(niter)
    cat(" ")
    mu <- subset(psooutid, iteration == niter)$pars
    lpbest <- subset(psooutid, iteration == niter)$maxes[1]
    try({
      out <- indmetrop(mcmciter, lpost, lposthess, mu, mu, 100, lpbest, datlist)
      acc <- mean(out$acc)
      accout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                 nswarm = nswarm, niter = niter, pso = alg, nbhd = nbhd,
                                 mcmc = "IMH", acc = acc, psoid = id),
                      accout)
    })
    try({
      out <- indmetropwithingibbs(mcmciter, mu, mu, df, lpbest, datlist)
      acc <- mean(out$acc)
      accout <- rbind(data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                 nswarm = nswarm, niter = niter, pso = alg, nbhd = nbhd,
                                 mcmc = "IMHwG", acc = acc, psoid = id),
                      accout)
    })
    write.csv(accout, file = "accout.csv", row.names=FALSE)
  }
}
