source("psofun.R")
source("electionfun.R")
source("electionmcmcfun.R")
load("datlistsmall.RData")
load("datlistplus.RData")

nrep <- 5
niter <- 10000
nswarm <- 100
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496
nbhdnames <- c("ring-1", "ring-3")
rates <- c(0.5)
cccs <- c(0.1)
alpha <- .2*niter
beta <- 1
pollpsoout <- NULL
nbhd <- list()
nbhd[[1]] <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)}) ## ring-1
nbhd[[2]] <- sapply(1:nswarm, function(x){return( (x + -3:3 - 1)%%nswarm + 1)}) ## ring-3
models <- c("small", "poll")
for(model in models){
  lpost <- switch(model,
                  small = gelmanlpost,
                  poll = gelmanpluslpost)
  datlist <- switch(model,
                    small = datlistsmall,
                    poll = datlistplus)
  npar <- 80 + 9*(model == "poll")
  for(m in 1:2){
    for(rep in 1:nrep){
      init <- matrix(runif(npar*nswarm, -100, 100), ncol = nswarm)
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
      print("BBPSOxp-MC")
      temp <- bbpso(niter, nswarm, 0, 1, init, nbhd[[m]], lpost, Inf, FALSE,
                    c(.5,.5), datlist = datlist)
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
      print("AT-PSO & AT-BBPSOxp-MC")
      for(rate in rates){
        for(ccc in cccs){
          print(paste(c(rate, ccc)))
          psotemp <- pso(niter, nswarm, 0.9, cognitive, social, init, nbhd[[m]],
                         lpost, datlist = datlist, tune = TRUE, style = "adaptive",
                         rate = rate, ccc = ccc)
          psotempout <- data.frame(model = model,
                                   algorithm = paste("AT-PSO", rate, ccc, sep="-"),
                                   nbhd = nbhdnames[m], rep = rep,
                                   iteration = 0:niter,
                                   maxes = psotemp$maxes)
          pollpsoout <- rbind(pollpsoout, psotempout)
          temp <- bbpso(niter, nswarm, 1, rate, init, nbhd[[m]], lpost, 1,
                        TRUE, c(0.5,0.5), datlist = datlist, ccc = ccc)
          psotempout <- data.frame(model = model,
                                   algorithm =
                                     paste("AT-BBPSOxp-MC", 1, rate, ccc, sep="-"),
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



