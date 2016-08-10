source("../psofun.R")
source("popfun.R")
source("mcmcfun.R")
load("popdat/popdat.RData")

load("psoout.RData")
psooutold2 <- psoout
save(psooutold2, file = "psooutold2.RData")

nbeta <- 1
nrep <- 20
ndeltasiid <- c(10, 30)
ndeltasfull <- c(5, 15)
niter <- 1000
nswarm <- 50
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496
nbhdnames <- c("ring-1", "ring-3")
rates <- c(0.3, 0.5)
cccs <- c(0.1, 0.00001)
alpha <- .2*niter
beta <- 1
ranefs <- c("iid", "full")
models <- c("pois", "lnorm")
psoout <- NULL
nbhd <- list()
nbhd[[1]] <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)}) ## ring-1
nbhd[[2]] <- sapply(1:nswarm, function(x){return( (x + -3:3 - 1)%%nswarm + 1)}) ## ring-3
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


{
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
        for(m in 1:2){
          for(psoinit in c("zero", "bfgs")){
            for(rep in 1:nrep){
              if(psoinit == "zero"){
                init <- matrix(runif(npar*nswarm, -100, 100), ncol = nswarm)
              } else {
                init <- matrix(runif(npar*nswarm, -1, 1), ncol = nswarm) + mubfgs
                init[,1] <- mubfgs
              }
              cat("idelta = ")
              cat(idelta)
              cat(", ")
              cat(model)
              cat(" ")
              cat(ranef)
              cat(" ")
              cat(ndelta)
              cat(", ")
              cat("nbhd = ")
              cat(nbhdnames[m])
              cat(", psoinit = ")
              cat(psoinit)
              cat(", rep = ")
              cat(rep)
              cat("\n")
              print("PSO")
              ## psotemp <- pso(niter, nswarm, inertia, cognitive, social, init, nbhd[[m]],
              ##                lpost, datlist = datlist)
              ## psotempout <- data.frame(model = model, ranef = ranef, ndelta = ndelta,
              ##                          algorithm = "PSO", init = psoinit,
              ##                          nbhd = nbhdnames[m], rep = rep,
              ##                          iteration = 0:niter,
              ##                          maxes = psotemp$maxes)
              ## psoout <- rbind(psoout, psotempout)
              ## print("BBPSOxp-MC")
              ## psotemp <- bbpso(niter, nswarm, 0, 1, init, nbhd[[m]], lpost, Inf, FALSE,
              ##               c(.5,.5), datlist = datlist)
              ## psotempout <- data.frame(model = model, ranef = ranef, ndelta = ndelta,
              ##                          algorithm = "BBPSOxp-MC", init = psoinit,
              ##                          nbhd = nbhdnames[m], rep = rep,
              ##                          iteration = 0:niter,
              ##                          maxes = psotemp$maxes)
              ## psoout <- rbind(psoout, psotempout)
              ## print("DI-PSO")
              ## psotemp <- pso(niter, nswarm, inertia, social, cognitive, init, nbhd[[m]],
              ##                lpost, datlist = datlist, tune = TRUE, style = "deterministic",
              ##                alpha = alpha, beta = beta)
              ## psotempout <- data.frame(model = model, ranef = ranef, ndelta = ndelta,
              ##                          algorithm = "DI-PSO", init = psoinit,
              ##                          nbhd = nbhdnames[m], rep = rep,
              ##                          iteration = 0:niter,
              ##                          maxes = psotemp$maxes)
              ## psoout <- rbind(psoout, psotempout)
              print("AT-PSO & AT-BBPSOxp-MC")
              for(rate in rates){
                for(ccc in cccs){
                  print(paste(c(rate, ccc)))
                  psotemp <- pso(niter, nswarm, 0.9, cognitive, social, init, nbhd[[m]],
                                 lpost, datlist = datlist, tune = TRUE, style = "adaptive",
                                 rate = rate, ccc = ccc)
                  psotempout <- data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                           algorithm = paste("AT-PSO", rate, ccc, sep="-"),
                                           init = psoinit,
                                           nbhd = nbhdnames[m], rep = rep,
                                           iteration = 0:niter,
                                           maxes = psotemp$maxes)
                  psoout <- rbind(psoout, psotempout)
                  psotemp <- bbpso(niter, nswarm, 1, rate, init, nbhd[[m]], lpost, 1,
                                TRUE, c(0.5,0.5), datlist = datlist, ccc = ccc)
                  psotempout <- data.frame(model = model, ranef = ranef, ndelta = ndelta,
                                           algorithm =
                                             paste("AT-BBPSOxp-MC", 1, rate, ccc, sep="-"),
                                           init = psoinit,
                                           nbhd = nbhdnames[m], rep = rep,
                                           iteration = 0:niter,
                                           maxes = psotemp$maxes)
                  psoout <- rbind(psoout, psotempout)
                }
              }
              save(psoout, file = "psoout.RData")
            }
          }
        }
      }
    }
  }
}


psooutnew2 <- psoout

psoout <- rbind(psooutnew2, subset(psooutold2, !(algorithm %in% levels(psooutnew2$algorithm))))

save(psooutnew2, file = "psooutnew2.RData")
save(psoout, file = "psoout.RData")
