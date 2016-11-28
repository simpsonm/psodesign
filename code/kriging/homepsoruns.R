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

nswarm <- 40
niter <- 2000
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

ndesign <- 100
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

homepsoouts <- foreach(i=1:40, .packages = c("sp", "maptools", "rgeos", "mnormt")) %dopar%
  psowrap(i, datlist, specs, seeds)

save(homepsoouts, file = "homepsoouts.RData")

stopImplicitCluster()
