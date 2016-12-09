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

ncores <- detectCores() - 4
registerDoParallel(ncores)

time <- 0:niter

nswarm <- 40
niter <- 2000
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
set.seed(324280)
seeds <- rnorm(length(specs))

for(objnum in 1:2){
  print("obj")
  print(objnum)
  if(objnum == 1){
    obj <- sig2fuk.mean
    objname <- "sig2fuk.mean"
  } else {
    obj <- sig2fuk.max
    objname <- "sig2fuk.max"
  }
  for(nbatch in nbatches){
    for(mutvar in mutvars){
      for(mutrate in mutrates){
        print("GA")
        print(c(nbatch, mutvar, mutrate))
        init <- replicate(nswarm, c(spsample(datlist$poly, ndesign, "random")@coords))
        temp <- ga(niter/nbatch, nbatch, floor(nswarm/2), nchrome, nrun, mutvar, mutrate,
                   lower, upper, obj,
                   init = init, boundaryfun = movetoboundary, datlist = datlist)
        algid <- paste("GA", nbatch, mutrate, mutvar, sep="-")
        tempdat <- data.frame(obj = objname, logpost = temp[["values"]],
                              time = time, algid = algid, type = "GA",
                              nbatch = nbatch, mutrate = mutrate, mutvar = mutvar,
                              rep = repl)
        gaout <- rbind(gaout, tempdat)
        temppar <- data.frame(obj = objname, logpost = temp[["value"]],
                              algid = algid, type = "GA", parset = NA, CF = NA,
                              style = NA, nbhd = NA, rep = repl,
                              parid = 1:(ndesign*2), par = temp[["par"]])
        parout2 <- rbind(parout2, temppar)
      }
    }
  }
  for(nexnbor in nexnbors){
    print("EX")
### this needs to be rewritten because of exch not having a fixed niter
    temp <- exch(ncand, obj, datlist$poly@coords, nexnbor, ndesign, datlist = datlist)
    algid <- paste("EX", nexnbor, sep="-")
    tempdat <- data.frame(obj = objname, logpost = temp[["values"]],
                          objcount = temp$objcount,
                          time = 1:length(temp$values), algid = algid,
                          type = "EX", nnbor = nexnbor, rep = repl)
    exout <- rbind(exout, tempdat)
    temppar <- data.frame(obj = objname, logpost = temp[["value"]],
                          algid = algid, type = "EX", parset = NA, CF = NA,
                          style = NA, nbhd = NA, rep = repl,
                          parid = 1:(ndesign*2), par = temp[["par"]])
    parout2 <- rbind(parout2, temppar)
  }
  write.csv(exout, file = "exsimsout.csv", row.names=FALSE)
  write.csv(gaout, file = "gasimsout.csv", row.names=FALSE)
  write.csv(parout2, file = "parsimsout2.csv", row.names=FALSE)
}




