source("../psofun.R")
source("electionfun.R")
source("mcmcfun.R")
load("datlistsmall.RData")
load("datlistplus.RData")

nrep <- 5
niter <- 10000
nswarm <- 100
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496
nbhdnames <- c("ring-1", "ring-3", "global")
rates <- c(0.5)
ccc <- c(0.1)
dfs <- 5
alpha <- .2*niter
beta <- 1
pollpsoout <- NULL
nbhd <- list()
nbhd[[1]] <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)}) ## ring-1
nbhd[[2]] <- sapply(1:nswarm, function(x){return( (x + -3:3 - 1)%%nswarm + 1)}) ## ring-3
nbhd[[3]] <- matrix(1:nswarm, ncol=nswarm, nrow=nswarm) ## global
models <- c("small", "poll")

for(model in models){
  lpost <- switch(model,
                  small = gelmanlpost,
                  poll = gelmanpluslpost)
  datlist <- switch(model,
                    small = datlistsmall,
                    poll = datlistplus)
  npar <- 80 + 9*(model == "poll")
  bfgsout <- optim(rep(0,npar), lpost, datlist = datlist,
                   control=list(fnscale=-1, reltol = .Machine$double.eps, maxit = 10000),
                   method = "BFGS")
  mufbgs <- bfgsout$par
  for(m in 1:2){
    for(rep in 1:nrep){
      init <- matrix(runif(npar*nswarm, -1, 1), ncol = nswarm) + bfgsout$par
      init[,1] <- bfgs$par
      cat(model)
      cat(" ")
      cat("nbhd = ")
      cat(nbhdnames[m])
      cat(", rep = ")
      cat(rep)
      cat("\n")
      print("PSO")
      psotemp <- pso(niter, nswarm, inertia, cognitive, social, init, nbhd[[m]],
                     lpost, datlist = datlist)
      psotempout <- data.frame(model = model,
                               algorithm = "PSO",
                               nbhd = nbhdnames[m], rep = rep,
                               iteration = 0:niter,
                               maxes = psotemp$maxes)
      pollpsoout <- rbind(pollpsoout, psotempout)
      print("BBPSO-MC")
      temp <- bbpso(niter, nswarm, 0, 1, init, nbhd[[m]], lpost, Inf, FALSE,
                    0, datlist = datlist)
      psotempout <- data.frame(model = model,
                               algorithm = "BBPSO-MC",
                               nbhd = nbhdnames[m], rep = rep,
                               iteration = 0:niter,
                               maxes = psotemp$maxes)
      pollpsoout <- rbind(pollpsoout, psotempout)
      print("BBPSOxp-MC")
      temp <- bbpso(niter, nswarm, 0, 1, init, nbhd[[m]], lpost, Inf, FALSE,
                    .5, datlist = datlist)
      psotempout <- data.frame(model = model,
                               algorithm = "BBPSOxp-MC",
                               nbhd = nbhdnames[m], rep = rep,
                               iteration = 0:niter,
                               maxes = psotemp$maxes)
      pollpsoout <- rbind(pollpsoout, psotempout)
      print("DI-PSO")
      psotemp <- pso(niter, nswarm, inertia, social, cognitive, init, nbhd[[m]],
                     lpost, datlist = datlist, tune = TRUE, style = "deterministic",
                     alpha = alpha, beta = beta)
      psotempout <- data.frame(model = model,
                               algorithm = "DI-PSO",
                               nbhd = nbhdnames[m], rep = rep,
                               iteration = 0:niter,
                               maxes = psotemp$maxes)
      pollpsoout <- rbind(pollpsoout, psotempout)
      print("AT-PSO & AT-BBPSO-MC & AT-BBPSOxp-MC")
      for(rate in rates){
        psotemp <- pso(niter, nswarm, 0.9, cognitive, social, init, nbhd[[m]],
                       lpost, datlist = datlist, tune = TRUE, style = "adaptive",
                       rate = rate, ccc = ccc)
        psotempout <- data.frame(model = model,
                                 algorithm = paste("AT-PSO", rate, ccc, sep="-"),
                                 nbhd = nbhdnames[m], rep = rep,
                                 iteration = 0:niter,
                                 maxes = psotemp$maxes)
        pollpsoout <- rbind(pollpsoout, psotempout)
        for(df in dfs){
          print(paste(c(rate, df)))
          psotemp <- bbpso(niter, nswarm, 1, rate, init, nbhd[[m]], lpost, df,
                        TRUE, 0, datlist = datlist, ccc = ccc)
          psotempout <- data.frame(model = model,
                                   algorithm =
                                     paste("AT-BBPSO-MC", df, rate, ccc, sep="-"),
                                   nbhd = nbhdnames[m], rep = rep,
                                   iteration = 0:niter,
                                   maxes = psotemp$maxes)
          pollpsoout <- rbind(pollpsoout, psotempout)
          psotemp <- bbpso(niter, nswarm, 1, rate, init, nbhd[[m]], lpost, df,
                        TRUE, 0.5, datlist = datlist, ccc = ccc)
          psotempout <- data.frame(model = model,
                                   algorithm =
                                     paste("AT-BBPSOxp-MC", df, rate, ccc, sep="-"),
                                   nbhd = nbhdnames[m], rep = rep,
                                   iteration = 0:niter,
                                   maxes = psotemp$maxes)
          pollpsoout <- rbind(pollpsoout, psotempout)
          
        }
      }
      write.csv(pollpsoout, file = "pollpsoout.csv", row.names=FALSE)
    }
  }
}



