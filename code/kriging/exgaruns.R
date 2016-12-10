#! /usr/bin/env Rscript
library(parallel)
library(doParallel)
library(foreach)
library(sp)
library(maptools)
library(rgeos)
library(mnormt)
source("../psofun.R")
source("krigingfun.R")
load("datlist.Rdata")
datlist$sppoly <- SpatialPolygons(list(b=Polygons(list(a=datlist$poly), "a")))

ncores <- 4
registerDoParallel(ncores)

niter <- 2000
time <- 0:niter

nswarm <- 40
ndesign <- 100
lower <- rep(apply(datlist$poly@coords, 2, min), each = ndesign)
upper <- rep(apply(datlist$poly@coords, 2, max), each = ndesign)

nbatches <- c(1,2)
nchrome <- 2
nrun <- ndesign
mutvars <- c(1,2)
mutrates <- c(1/100, 1/10)
nexnbors <- c(5, 10)
ncand <- 2000
parout2 <- NULL
gaout <- NULL
exout <- NULL

objnames <- c("sig2fuk.mean", "sig2fuk.max")

## define specs
gaspecs <- c(outer(c(outer(objnames, nbatches, paste, sep="-")),
                 c(outer(mutvars, mutrates, paste, sep="-")), paste, sep = "-"))
exspecs <- c(outer(objnames, nexnbors, paste, sep = "-"))
set.seed(324280)
gaseeds <- rnorm(length(gaspecs))
exseeds <- rnorm(length(exspecs))

gawrap <- function(i, datlist, specs, seeds){
  spec <- specs[i]
  set.seed(seeds[i])
  splt <- strsplit(spec, "-")[[1]]
  objname <- splt[1]
  obj <- switch(objname, sig2fuk.mean = sig2fuk.mean, sig2fuk.max = sig2fuk.max)
  nbatch <- as.numeric(splt[2])
  mutvar <- as.numeric(splt[3])
  mutrate <- as.numeric(splt[4])
  repl <- 1
  init <- replicate(floor(nswarm/2), c(spsample(datlist$poly, ndesign, "random")@coords))
  temp <- ga(niter/nbatch, nbatch, floor(nswarm/2), nchrome, nrun, mutvar, mutrate,
             lower, upper, obj,
             init = init, boundaryfun = movetoboundary, datlist = datlist)
  algid <- paste("GA", nbatch, mutrate, mutvar, sep="-")
  tempdat <- data.frame(obj = objname, logpost = temp[["values"]],
                        time = time, algid = algid, type = "GA",
                        nbatch = nbatch, mutrate = mutrate, mutvar = mutvar,
                        rep = repl)
  temppar <- data.frame(obj = objname, logpost = temp[["value"]],
                        algid = algid, type = "GA", parset = NA, CF = NA,
                        style = NA, nbhd = NA, rep = repl,
                        parid = 1:(ndesign*2), par = temp[["par"]])
  out <- list(values = tempdat, pars = temppar)
  save(out, file = paste("GA-", spec, ".RData", sep = ""))
  print(paste("Spec ", i, " finished. Spec: GA-", spec, sep=""))
  return(out)
}

exwrap <- function(i, datlist, specs, seeds){
  spec <- specs[i]
  set.seed(seeds[i])
  splt <- strsplit(spec, "-")[[1]]
  objname <- splt[1]
  obj <- switch(objname, sig2fuk.mean = sig2fuk.mean, sig2fuk.max = sig2fuk.max)
  nexnbor <- splt[2]
  temp <- exch(ncand, obj, datlist$poly@coords, nexnbor, ndesign, datlist = datlist)
  algid <- paste("EX", nexnbor, sep="-")
  repl <- 1
  tempdat <- data.frame(obj = objname, logpost = temp[["values"]],
                        objcount = temp$objcount,
                        time = 1:length(temp$values), algid = algid,
                        type = "EX", nnbor = nexnbor, rep = repl)
  temppar <- data.frame(obj = objname, logpost = temp[["value"]],
                        algid = algid, type = "EX", parset = NA, CF = NA,
                        style = NA, nbhd = NA, rep = repl,
                        parid = 1:(ndesign*2), par = temp[["par"]])
  out <- list(values = tempdat, pars = temppar)
  save(out, file = paste("EX-", spec, ".RData", sep = ""))
  print(paste("Spec ", i, " finished. Spec: EX-", spec, sep=""))
  return(out)
}

gaouts <- foreach(i=1:length(gaspecs), .packages = c("sp", "maptools", "rgeos", "mnormt")) %dopar%
  gawrap(i, datlist, gaspecs, gaseeds)

save(gaouts, file = "gaouts.RData")

exouts <- foreach(i=1:length(exspecs), .packages = c("sp", "maptools", "rgeos", "mnormt")) %dopar%
  exwrap(i, datlist, exspecs, exseeds)

save(exouts, file = "exouts.RData")


stopImplicitCluster()
