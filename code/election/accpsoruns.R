source("../psofun.R")
source("electionfun.R")
source("mcmcfun.R")
load("datlistsmall.RData")
load("datlistplus.RData")

niter <- 20000
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
accpollpsoout <- NULL
nbhd <- list()
nbhd[[1]] <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)}) ## ring-1
nbhd[[2]] <- sapply(1:nswarm, function(x){return( (x + -3:3 - 1)%%nswarm + 1)}) ## ring-3
nbhd[[3]] <- matrix(1:nswarm, ncol=nswarm, nrow=nswarm) ## global
models <- c("small", "poll")

psoid <- 1
for(model in models){
  lpost <- switch(model,
                  small = gelmanlpost,
                  poll = gelmanpluslpost)
  datlist <- switch(model,
                    small = datlistsmall,
                    poll = datlistplus)
  npar <- 80 + 9*(model == "poll")
  for(m in 1:3){
    init <- matrix(runif(npar*nswarm, -100, 100), ncol = nswarm)
    cat(model)
    cat(" ")
    cat("nbhd = ")
    cat(nbhdnames[m])
    cat("\n")
    print("PSO")
    psotemp <- pso(niter, nswarm, inertia, cognitive, social, init, nbhd[[m]],
                   lpost, datlist = datlist)
    psotempout <- data.frame(model = model,
                             algorithm = "PSO",
                             nbhd = nbhdnames[m],
                             iteration = rep(0:niter, each = npar),
                             maxes = rep(psotemp$maxes, each = npar),
                             pars = c(psotemp$argmaxes), psoid = psoid)
    psoid <- psoid + 1    
    accpollpsoout <- rbind(accpollpsoout, psotempout)
    print("BBPSO-MC")
    temp <- bbpso(niter, nswarm, 0, 1, init, nbhd[[m]], lpost, Inf, FALSE,
                  0, datlist = datlist)
    psotempout <- data.frame(model = model,
                             algorithm = "BBPSO-MC",
                             nbhd = nbhdnames[m],
                             iteration = rep(0:niter, each = npar),
                             maxes = rep(psotemp$maxes, each = npar),
                             pars = c(psotemp$argmaxes), psoid = psoid)
    psoid <- psoid + 1    
    accpollpsoout <- rbind(accpollpsoout, psotempout)
    print("BBPSOxp-MC")
    temp <- bbpso(niter, nswarm, 0, 1, init, nbhd[[m]], lpost, Inf, FALSE,
                  .5, datlist = datlist)
    psotempout <- data.frame(model = model,
                             algorithm = "BBPSOxp-MC",
                             nbhd = nbhdnames[m],
                             iteration = rep(0:niter, each = npar),
                             maxes = rep(psotemp$maxes, each = npar),
                             pars = c(psotemp$argmaxes), psoid = psoid)
    psoid <- psoid + 1    
    accpollpsoout <- rbind(accpollpsoout, psotempout)
    print("DI-PSO")
    psotemp <- pso(niter, nswarm, inertia, social, cognitive, init, nbhd[[m]],
                   lpost, datlist = datlist, tune = TRUE, style = "deterministic",
                   alpha = alpha, beta = beta)
    psotempout <- data.frame(model = model,
                             algorithm = "DI-PSO",
                             nbhd = nbhdnames[m],
                             iteration = rep(0:niter, each = npar),
                             maxes = rep(psotemp$maxes, each = npar),
                             pars = c(psotemp$argmaxes), psoid = psoid)
    psoid <- psoid + 1
    accpollpsoout <- rbind(accpollpsoout, psotempout)
    print("AT-PSO & AT-BBPSO-MC & AT-BBPSOxp-MC")
    for(rate in rates){
      psotemp <- pso(niter, nswarm, 0.9, cognitive, social, init, nbhd[[m]],
                     lpost, datlist = datlist, tune = TRUE, style = "adaptive",
                     rate = rate, ccc = ccc)
      psotempout <- data.frame(model = model,
                               algorithm = paste("AT-PSO", rate, ccc, sep="-"),
                               nbhd = nbhdnames[m],
                               iteration = rep(0:niter, each = npar),
                               maxes = rep(psotemp$maxes, each = npar),
                               pars = c(psotemp$argmaxes), psoid = psoid)
      psoid <- psoid + 1      
      accpollpsoout <- rbind(accpollpsoout, psotempout)
      for(df in dfs){
        print(paste(c(rate, df)))
        psotemp <- bbpso(niter, nswarm, 1, rate, init, nbhd[[m]], lpost, df,
                         TRUE, 0, datlist = datlist, ccc = ccc)
        psotempout <- data.frame(model = model,
                                 algorithm =
                                   paste("AT-BBPSO-MC", df, rate, ccc, sep="-"),
                                 nbhd = nbhdnames[m],
                                 iteration = rep(0:niter, each = npar),
                                 maxes = rep(psotemp$maxes, each = npar),
                                 pars = c(psotemp$argmaxes), psoid = psoid)
        psoid <- psoid + 1        
        accpollpsoout <- rbind(accpollpsoout, psotempout)
        psotemp <- bbpso(niter, nswarm, 1, rate, init, nbhd[[m]], lpost, df,
                         TRUE, 0.5, datlist = datlist, ccc = ccc)
        psotempout <- data.frame(model = model,
                                 algorithm =
                                   paste("AT-BBPSOxp-MC", df, rate, ccc, sep="-"),
                                 nbhd = nbhdnames[m],
                                 iteration = rep(0:niter, each = npar),
                                 maxes = rep(psotemp$maxes, each = npar),
                                 pars = c(psotemp$argmaxes), psoid = psoid)
        psoid <- psoid + 1        
        accpollpsoout <- rbind(accpollpsoout, psotempout)
        
      }
    }
    write.csv(accpollpsoout, file = "accpollpsoout.csv", row.names=FALSE)
  }
}



accout <- NULL
niters <- c(2000, 4000, 8000, 10000)
mcmciter <- 10000
df <- 100
for(id in 1:max(accpsoout$psoid)){
  cat("\n PSOID = ")
  cat(id)
  cat(", niter = ")
  psooutid <- subset(accpollpsoout, psoid == id)
  model <- psooutid$model[1]
  nbhd <- psooutid$nbhd[1]
  alg <- psooutid$algorithm[1]
  lpost <- switch(as.character(model),
                  small = gelmanlpost,
                  poll = gelmanpluslpost)
  lposthess <- switch(as.character(model),
                  small = gelmanlposthess,
                  poll = gelmanpluslposthess)
  datlist <- switch(as.character(model),
                    small = datlistsmall,
                    poll = datlistplus)
  npar <- 80 + 9*(as.character(model) == "poll")
  indmetropwithingibbs <- switch(as.character(model),
                                 small = gelmanindwithingibbs,
                                 poll = gelmanplusindwithingibbs)
  for(niter in niters){
    cat(niter)
    cat(" ")
    mu <- subset(psooutid, iteration == niter)$pars
    lpbest <- subset(psooutid, iteration == niter)$maxes[1]
    try({
      out <- indmetrop(mcmciter, lpost, lposthess, mu, mu, 100, lpbest, datlist)
      acc <- mean(out$acc)
      accout <- rbind(data.frame(model = model, 
                                 nswarm = nswarm, niter = niter, pso = alg, nbhd = nbhd,
                                 mcmc = "IMH", acc = acc, psoid = id),
                      accout)
    })
    try({
      out <- indmetropwithingibbs(mcmciter, mu, mu, df, datlist)
      acc <- mean(out$acc)
      accout <- rbind(data.frame(model = model,
                                 nswarm = nswarm, niter = niter, pso = alg, nbhd = nbhd,
                                 mcmc = "IMHwG", acc = acc, psoid = id),
                      accout)
    })
    write.csv(accout, file = "accout.csv", row.names=FALSE)
  }
}
  
