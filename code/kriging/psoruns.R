library(sp)
library(maptools)
source("../psofun.R")
source("krigingfun.R")
load("datlist.Rdata")

nrep <- 10
niter <- 500
nswarm <- 20
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496
nbhdnames <- c("ring-1", "ring-3", "global")
rates <- c(0.5)
dfs <- c(5)
ccc <- c(0.1)
alpha <- .2*niter
beta <- 1
psoout <- NULL
nbhd <- list()
nbhd[[1]] <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)}) ## ring-1
nbhd[[2]] <- sapply(1:nswarm, function(x){return( (x + -3:3 - 1)%%nswarm + 1)}) ## ring-3
nbhd[[3]] <- matrix(1:nswarm, ncol=nswarm, nrow=nswarm) ## global

ndesign <- 5
npar <- 2*ndesign
inits <- list()
inits[[1]] <- replicate(nswarm, c(spsample(datlist$poly, ndesign, "random")@coords))
idxs <- replicate(nswarm, sample(1:nrow(datlist$tt), ndesign))
inits2 <- matrix(0, npar, nswarm)
for(i in 1:nswarm){
  inits2[,i] <- c(datlist$tt[idxs[,i],])
}
inits[[2]] <- inits2

system.time({
for(m in 1:3){
  for(i in 1:2){
    for(j in 1:2){
      if(j == 1){
        obj <- negsig2sk.mean
      } else {
        obj <- negsig2sk.min
      }
      for(rep in 1:nrep){
        cat("nbhd = ")
        cat(nbhdnames[m])
        cat(", init = ")
        cat(i)
        cat(", obj = ")
        cat(j)
        cat(", rep = ")
        cat(rep)
        cat("\n")
        print("PSO")
        psotemp <- pso(niter, nswarm, inertia, cognitive, social, inits[[i]], nbhd[[m]],
                       obj, datlist = datlist)
        psotempout <- data.frame(algorithm = "PSO", obj = j,
                                 init = i, nbhd = nbhdnames[m], rep = rep,
                                 iteration = 0:niter, maxes = psotemp$maxes)
        psoout <- rbind(psoout, psotempout)
        print("BBPSO-MC")
        psotemp <- bbpso(niter, nswarm, 0, 1, inits[[i]], nbhd[[m]], obj, Inf, FALSE,
                         0, datlist = datlist)
        psotempout <- data.frame(algorithm = "BBPSO-MC", obj = j,
                                 init = i, nbhd = nbhdnames[m], rep = rep,
                                 iteration = 0:niter, maxes = psotemp$maxes)
        psoout <- rbind(psoout, psotempout)            
        print("BBPSOxp-MC")
        psotemp <- bbpso(niter, nswarm, 0, 1, inits[[i]], nbhd[[m]], obj, Inf, FALSE,
                         .5, datlist = datlist)
        psotempout <- data.frame(algorithm = "BBPSOxp-MC", obj = j,
                                 init = i, nbhd = nbhdnames[m], rep = rep,
                                 iteration = 0:niter, maxes = psotemp$maxes)
        psoout <- rbind(psoout, psotempout)
        print("DI-PSO")
        psotemp <- pso(niter, nswarm, inertia, social, cognitive, inits[[i]], nbhd[[m]],
                       obj, datlist = datlist, tune = TRUE, style = "deterministic",
                       alpha = alpha, beta = beta)
        psotempout <- data.frame(algorithm = "DI-PSO", obj = j,
                                 init = i, nbhd = nbhdnames[m], rep = rep,
                                 iteration = 0:niter, maxes = psotemp$maxes)
        psoout <- rbind(psoout, psotempout)
        print("AT-PSO & AT-BBPSO-MC & AT-BBPSOxp-MC")
        for(rate in rates){
          psotemp <- pso(niter, nswarm, 0.9, cognitive, social, inits[[i]], nbhd[[m]],
                         obj, datlist = datlist, tune = TRUE, style = "adaptive",
                         rate = rate, ccc = ccc)
          psotempout <- data.frame(algorithm = paste("AT-PSO", rate, ccc, sep="-"),
                                   obj = j,
                                   init = i, nbhd = nbhdnames[m], rep = rep,
                                   iteration = 0:niter, maxes = psotemp$maxes)
          psoout <- rbind(psoout, psotempout)
          for(df in dfs){
            print(paste(c(rate, df)))
            psoout <- rbind(psoout, psotempout)
            psotemp <- bbpso(niter, nswarm, 1, rate, inits[[i]], nbhd[[m]], obj, df,
                             TRUE, 0, datlist = datlist, ccc = ccc)
            psotempout <- data.frame(algorithm = paste("AT-BBPSO-MC", df, rate, ccc, sep="-"),
                                     obj = j,
                                     init = i, nbhd = nbhdnames[m], rep = rep,
                                     iteration = 0:niter, maxes = psotemp$maxes)
            psoout <- rbind(psoout, psotempout)
            psotemp <- bbpso(niter, nswarm, 1, rate, inits[[i]], nbhd[[m]], obj, df,
                             TRUE, 0.5, datlist = datlist, ccc = ccc)
            psotempout <- data.frame(algorithm = paste("AT-BBPSOxp-MC", df, rate, ccc, sep="-"),
                                     obj = j,
                                     init = i, nbhd = nbhdnames[m], rep = rep,
                                     iteration = 0:niter, maxes = psotemp$maxes)
            psoout <- rbind(psoout, psotempout)
          }
        }
        write.csv(psoout, file = "psoout.csv", row.names=FALSE)
      }
    }
  }
}
})
