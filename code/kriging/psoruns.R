library(sp)
library(maptools)
source("../psofun.R")
source("krigingfun.R")
load("datlist.Rdata")


nswarm <- 40
niter <- 500
nrep <- 5
inertias <- c(0.7298, 1/(log(2)*2))
cognitives <- c(1.496, log(2) + 1/2)
socials <- c(1.496, log(2) + 1/2)
nnbors <- c(1, 3, 40)
alpha <- 0.2*niter
beta <- 2
rates <- c(0.3, 0.5)
ccc <- 0.1
df <- 1
pcuts <- c(0, 0.5)
sig0 <- 1
inertia0 <- 1.2

ndesign <- 5
lower <- rep(apply(datlist$poly@coords, 2, min), each = ndesign)
upper <- rep(apply(datlist$poly@coords, 2, max), each = ndesign)

time <- 0:niter

set.seed(3453)


psoout <- NULL
cat("\n")
for(repl in 1:nrep){
  for(nnbor in nnbors){
    for(CF in c(TRUE, FALSE)){
      for(parset in 1:2){
        for(objnum in 1:2){
          if(objnum == 1){
            obj <- sig2fuk.mean
            objname <- "sig2fuk.mean"
          } else {
            obj <- sig2fuk.max
            objname <- "sig2fuk.max"
          }
          cat("rep = "); cat(repl); cat("; obj = "); cat(objname); cat("; nnbor = "); cat(nnbor)
          cat("; CF = "); cat(CF); cat("; parset = "); cat(parset); cat("\n")
          for(style in c("AT1", "AT2")){
            system.time({
            rate <- ifelse(style=="AT1", rates[1], rates[2])
            pcut <- pcuts[parset]
            temp <- sbbpso(niter, nswarm, nnbor, sig0,
                           pcut = pcut, CF = CF, AT = TRUE, rate = rate, df = df, ccc = 0.1,
                           obj = obj, datlist = datlist)
            algid <- paste("BBPSO", parset, ifelse(CF, "CF", "notCF"), style, sep = "-")
            tempdat <- data.frame(obj = objname, logpost = temp[["values"]],
                                  argnorm = apply(temp[["pars"]], 2, normvec),
                                  time = time, algid = algid,
                                  type = "BBPSO", parset = parset, CF = CF,
                                  style = style, nbhd = nnbor, rep = repl,
                                  inertias = temp$sigs)
            psoout <- rbind(psoout, tempdat)
            })
          }
          for(style in c("CI", "DI", "AT1", "AT2")){
            rate <- ifelse(style=="AT1", rates[1], rates[2])
            c.in <- ifelse(style=="CI", inertias[parset], inertia0)
            c.co <- cognitives[parset]
            c.so <- socials[parset]
            temp <- spso(niter, nswarm, nnbor, c.in, c.co, c.so, obj, lower, upper,
                         style = substr(style, 1, 2), CF = CF, alpha = alpha, beta = beta,
                         rate = rate, ccc = ccc, datlist = datlist)
            algid <- paste("PSO", parset, ifelse(CF, "CF", "notCF"), style, sep = "-")
            tempdat <- data.frame(obj = objname, logpost = temp[["values"]],
                                  argnorm = apply(temp[["pars"]], 2, normvec),
                                  time = time, algid = algid,
                                  type = "PSO", parset = parset, CF = CF,
                                  style = style, nbhd = nnbor, rep = repl,
                                  inertias = temp$inertias)
            psoout <- rbind(psoout, tempdat)
          }
          write.csv(psoout, file = "psosimsout.csv", row.names=FALSE)
        }
      }
    }
  }
}


nbatches <- c(1,2)
nchrome <- 2
nrun <- ndesign
mutvars <- c(1,2)
mutrates <- c(1/100, 1/10)
nexnbors <- c(1, 3, 5)
ncand <- 2000
gaout <- NULL
exout <- NULL
for(repl in 1:nrep){
  for(objnum in 1:2){
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
          temp <- ga(niter/nbatch, nbatch, floor(nswarm/2), nchrome, nrun, mutvar, mutrate,
                     lower, upper, obj, datlist=datlist)
          algid <- paste("GA", nbatch, mutrate, mutvar, sep="-")
          tempdat <- data.frame(obj = objname, logpost = temp[["values"]],
                                argnorm = apply(temp[["pars"]], 2, normvec),
                                time = time, algid = algid, type = "GA",
                                nbatch = nbatch, mutrate = mutrate, mutvar = mutvar,
                                rep = repl)
          gaout <- rbind(gaout, tempdat)
        }
      }
    }
    for(nexnbor in nexnbors){
      ### this needs to be rewritten because of exch not having a fixed niter
      temp <- exch(ncand, obj, datlist$poly@coords, nexnbor, ndesign, datlist = datlist)
      algid <- paste("EX", nexnbor, sep="-")
      tempdat <- data.frame(obj = objname, logpost = temp[["values"]],
                            argnorm = apply(temp[["pars"]], 2, normvec),
                            objcount = temp$objcount,
                            time = 1:temp$count, algid = algid,
                            type = "GA", nnbor = nexnbor, rep = repl)

      exout <- rbind(exout, tempdat)
    }
    write.csv(exout, file = "exsimsout.csv", row.names=FALSE)
    write.csv(gaout, file = "gasimsout.csv", row.names=FALSE)
  }
}




