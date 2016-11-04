library(parallel)
library(doParallel)
library(foreach)
library(sp)
library(maptools)
source("../psofun.R")
source("krigingfun.R")
load("datlist.Rdata")
datlist$sppoly <- SpatialPolygons(list(b=Polygons(list(a=datlist$poly), "a")))

ncores <- 20  ## how many cores do we want to request
registerDoParallel(ncores)

nswarm <- 40
niter <- 1000
nrep <- 1
inertias <- c(0.7298, 1/(log(2)*2))
cognitives <- c(1.496, log(2) + 1/2)
socials <- c(1.496, log(2) + 1/2)
nnbors <- c(3, 40)
alpha <- 0.2*niter
beta <- 2
rates <- c(0.3, 0.5)
ccc <- 0.1
df <- 1
pcuts <- c(0, 0.5)
sig0 <- 1
inertia0 <- 1.2

ndesign <- 20
lower <- rep(apply(datlist$poly@coords, 2, min), each = ndesign)
upper <- rep(apply(datlist$poly@coords, 2, max), each = ndesign)

time <- 0:niter

parsets <- 1:2
CFs <- c("CF", "nCF")
objnames <- c("sig2fuk.mean", "sig2fuk.max")

specs <- c(outer(c(outer(c(outer(paste(c("CI", "DI", "AT1", "AT2"), "PSO", sep="-"), parsets, paste, sep = "-")),
                        c(outer(CFs, nnbors, paste, sep="-")), paste, sep="-"),
                  c(outer(c(outer(paste(c("AT1", "AT2"), "BBPSO", sep="-"), parsets, paste, sep = "-")),
                          c(outer(CFs, nnbors, paste, sep="-")), paste, sep="-"))), objnames, paste, sep="-"))
set.seed(234132)
seeds <- rnorm(length(specs))

psowrap <- function(i, datlist, specs, seeds){
  spec <- specs[i]
  set.seed(seeds[i])
  splt <- strsplit(spec, "-")[[1]]
  style <- splt[1]
  alg <- splt[2]
  parset <- as.numeric(splt[3])
  CF <- splt[4]=="CF"
  nnbor <- as.numeric(splt[5])
  objname <- splt[6]
  obj <- switch(objname, sig2fuk.mean = sig2fuk.mean, sig2fuk.max = sig2fuk.max)
  repl <- 1
  if(alg == "BBPSO"){
    rate <- ifelse(style=="AT1", rates[1], rates[2])
    pcut <- pcuts[parset]
    init <- replicate(nswarm, c(spsample(datlist$poly, ndesign, "random")@coords))
    temp <- sbbpso(niter, nswarm, nnbor, sig0, obj, lower, upper,
                   pcut = pcut, CF = CF, AT = TRUE, rate = rate, df = df, ccc = 0.1,
                   init = init, boundaryfun = movetoboundary, datlist = datlist)
    algid <- paste("BBPSO", parset, ifelse(CF, "CF", "notCF"), style, sep = "-")
    tempdat <- data.frame(obj = objname, logpost = temp[["values"]],
                          time = time, algid = algid,
                          type = "BBPSO", parset = parset, CF = CF,
                          style = style, nbhd = nnbor, rep = repl,
                          inertias = temp$sigs)
    temppar <- data.frame(obj = objname, logpost = temp[["value"]],
                          algid = algid, type = "BBPSO", parset = parset, CF = CF,
                          style = style, nbhd = nnbor, rep = repl,
                          parid = 1:(ndesign*2), par = temp[["par"]])
  } else {
    rate <- ifelse(style=="AT1", rates[1], rates[2])
    c.in <- ifelse(style=="CI", inertias[parset], inertia0)
    c.co <- cognitives[parset]
    c.so <- socials[parset]
    init <- replicate(nswarm, c(spsample(datlist$poly, ndesign, "random")@coords))
    temp <- spso(niter, nswarm, nnbor, c.in, c.co, c.so, obj, lower, upper,
                 style = substr(style, 1, 2), CF = CF, alpha = alpha, beta = beta,
                 rate = rate, ccc = ccc,
                 init = init, boundaryfun = movetoboundary, datlist = datlist)
    algid <- paste("PSO", parset, ifelse(CF, "CF", "notCF"), style, sep = "-")
    tempdat <- data.frame(obj = objname, logpost = temp[["values"]],
                          time = time, algid = algid,
                          type = "PSO", parset = parset, CF = CF,
                          style = style, nbhd = nnbor, rep = repl,
                          inertias = temp$inertias)
    temppar <- data.frame(obj = objname, logpost = temp[["value"]],
                          algid = algid, type = "PSO", parset = parset, CF = CF,
                          style = style, nbhd = nnbor, rep = repl,
                          parid = 1:(ndesign*2), par = temp[["par"]])
  }
  out <- list(values = tempdat, pars = temppar)
  save(out, file = paste(spec, ".RData", sep = ""))
  print(paste("Spec ", i, " finished. Spec: ", spec, sep=""))
  return(out)
}

psoouts <- foreach(i=1:length(specs), .packages = c("sp", "maptools", "rgeos", "mnormt")) %dopar%
  psowrap(i, datlist, specs, seeds)

save(psoouts, file = "psoouts.RData")

stopImplicitCluster()

## psoout <- NULL
## parout <- NULL
## cat("\n")
## for(repl in 1:nrep){
##   for(nnbor in nnbors){
##     for(CF in c(TRUE, FALSE)){
##       for(parset in 1:2){
##         for(objnum in 1:2){
##           if(objnum == 1){
##             obj <- sig2fuk.mean
##             objname <- "sig2fuk.mean"
##           } else {
##             obj <- sig2fuk.max
##             objname <- "sig2fuk.max"
##           }
##           cat("rep = "); cat(repl); cat("; obj = "); cat(objname); cat("; nnbor = "); cat(nnbor)
##           cat("; CF = "); cat(CF); cat("; parset = "); cat(parset); cat("\n")
##           for(style in c("AT1", "AT2")){
##             rate <- ifelse(style=="AT1", rates[1], rates[2])
##             pcut <- pcuts[parset]
##             init <- replicate(nswarm, c(spsample(datlist$poly, ndesign, "random")@coords))
##             temp <- sbbpso(niter, nswarm, nnbor, sig0, obj, lower, upper,
##                            pcut = pcut, CF = CF, AT = TRUE, rate = rate, df = df, ccc = 0.1,
##                            init = init, boundaryfun = movetoboundary, datlist = datlist)
##             algid <- paste("BBPSO", parset, ifelse(CF, "CF", "notCF"), style, sep = "-")
##             tempdat <- data.frame(obj = objname, logpost = temp[["values"]],
##                                   time = time, algid = algid,
##                                   type = "BBPSO", parset = parset, CF = CF,
##                                   style = style, nbhd = nnbor, rep = repl,
##                                   inertias = temp$sigs)
##             psoout <- rbind(psoout, tempdat)
##             temppar <- data.frame(obj = objname, logpost = temp[["value"]],
##                                   algid = algid, type = "BBPSO", parset = parset, CF = CF,
##                                   style = style, nbhd = nnbor, rep = repl,
##                                   parid = 1:(ndesign*2), par = temp[["par"]])
##             parout <- rbind(parout, temppar)
##           }
##           for(style in c("CI", "DI", "AT1", "AT2")){
##             rate <- ifelse(style=="AT1", rates[1], rates[2])
##             c.in <- ifelse(style=="CI", inertias[parset], inertia0)
##             c.co <- cognitives[parset]
##             c.so <- socials[parset]
##             init <- replicate(nswarm, c(spsample(datlist$poly, ndesign, "random")@coords))
##             temp <- spso(niter, nswarm, nnbor, c.in, c.co, c.so, obj, lower, upper,
##                          style = substr(style, 1, 2), CF = CF, alpha = alpha, beta = beta,
##                          rate = rate, ccc = ccc,
##                          init = init, boundaryfun = movetoboundary, datlist = datlist)
##             algid <- paste("PSO", parset, ifelse(CF, "CF", "notCF"), style, sep = "-")
##             tempdat <- data.frame(obj = objname, logpost = temp[["values"]],
##                                   time = time, algid = algid,
##                                   type = "PSO", parset = parset, CF = CF,
##                                   style = style, nbhd = nnbor, rep = repl,
##                                   inertias = temp$inertias)
##             psoout <- rbind(psoout, tempdat)
##             temppar <- data.frame(obj = objname, logpost = temp[["value"]],
##                                   algid = algid, type = "PSO", parset = parset, CF = CF,
##                                   style = style, nbhd = nnbor, rep = repl,
##                                   parid = 1:(ndesign*2), par = temp[["par"]])
##             parout <- rbind(parout, temppar)
##           }
##           write.csv(psoout, file = "psosimsout.csv", row.names=FALSE)
##           write.csv(parout, file = "parsimsout.csv", row.names=FALSE)
##         }
##       }
##     }
##   }
## }


## set.seed(324280)

## nbatches <- c(1,2)
## nchrome <- 2
## nrun <- ndesign
## mutvars <- c(1,2)
## mutrates <- c(1/100, 1/10)
## nexnbors <- c(5, 10)
## ncand <- 2000
## parout2 <- NULL
## gaout <- NULL
## exout <- NULL
## for(repl in 1:nrep){
##   print("rep")
##   print(repl)
##   for(objnum in 1:2){
##     print("obj")
##     print(objnum)
##     if(objnum == 1){
##       obj <- sig2fuk.mean
##       objname <- "sig2fuk.mean"
##     } else {
##       obj <- sig2fuk.max
##       objname <- "sig2fuk.max"
##     }
##     for(nbatch in nbatches){
##       for(mutvar in mutvars){
##         for(mutrate in mutrates){
##           print("GA")
##           print(c(nbatch, mutvar, mutrate))
##           init <- replicate(nswarm, c(spsample(datlist$poly, ndesign, "random")@coords))
##           temp <- ga(niter/nbatch, nbatch, floor(nswarm/2), nchrome, nrun, mutvar, mutrate,
##                      lower, upper, obj,
##                      init = init, boundaryfun = movetoboundary, datlist = datlist)
##           algid <- paste("GA", nbatch, mutrate, mutvar, sep="-")
##           tempdat <- data.frame(obj = objname, logpost = temp[["values"]],
##                                 time = time, algid = algid, type = "GA",
##                                 nbatch = nbatch, mutrate = mutrate, mutvar = mutvar,
##                                 rep = repl)
##           gaout <- rbind(gaout, tempdat)
##           temppar <- data.frame(obj = objname, logpost = temp[["value"]],
##                                 algid = algid, type = "GA", parset = NA, CF = NA,
##                                 style = NA, nbhd = NA, rep = repl,
##                                 parid = 1:(ndesign*2), par = temp[["par"]])
##           parout2 <- rbind(parout2, temppar)
##         }
##       }
##     }
##     for(nexnbor in nexnbors){
##       print("EX")
##       ### this needs to be rewritten because of exch not having a fixed niter
##       temp <- exch(ncand, obj, datlist$poly@coords, nexnbor, ndesign, datlist = datlist)
##       algid <- paste("EX", nexnbor, sep="-")
##       tempdat <- data.frame(obj = objname, logpost = temp[["values"]],
##                             objcount = temp$objcount,
##                             time = 1:length(temp$values), algid = algid,
##                             type = "EX", nnbor = nexnbor, rep = repl)
##       exout <- rbind(exout, tempdat)
##       temppar <- data.frame(obj = objname, logpost = temp[["value"]],
##                             algid = algid, type = "EX", parset = NA, CF = NA,
##                             style = NA, nbhd = NA, rep = repl,
##                             parid = 1:(ndesign*2), par = temp[["par"]])
##       parout2 <- rbind(parout2, temppar)
##     }
##     write.csv(exout, file = "exsimsout.csv", row.names=FALSE)
##     write.csv(gaout, file = "gasimsout.csv", row.names=FALSE)
##     write.csv(parout2, file = "parsimsout2.csv", row.names=FALSE)
##   }
## }




