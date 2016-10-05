source("../psofun.R")
source("testfuns.R")

library(pso)

np <- 10
nswarm <- 40
niter <- 1000
inertia <- 1/(log(2)*2)
cognitive <- log(2) + 1/2
social <- log(2) + 1/2
lower <- rep(-100, np)
upper <- rep(100, np)
nnbor <- round((1-(1-1/nswarm)^3)*nswarm,0)

test1 <- psoptim(par = rep(NA, np), fn = negfwrap, opt = 1, lower = -100, upper = 100,
                control = list(maxit = niter, s = nswarm, vectorize = FALSE))

test2 <- spso2007(niter, nswarm, nnbor, inertia, cognitive, social, negfwrap, lower,
                  upper, opt = 1)


test2 <- psoptim(par = rep(NA, np), fn = negfwrap, opt = 1, lower = -100, upper = 100,
                 control = list(maxit = niter, s = nswarm, p = 1, vectorize = TRUE,
                                type = "SPSO2011"))

test3 <- psoptim(par = rep(NA, np), fn = negfwrap, opt = 1, lower = -100, upper = 100,
                 control = list(maxit = niter, s = nswarm, p = 1, vectorize = FALSE,
                                type = "SPSO2011"))

test4 <- psoptim(par = rep(NA, np), fn = negfwrap, opt = 1, lower = -100, upper = 100,
                 control = list(maxit = niter, s = nswarm, p = 1, vectorize = TRUE,
                                w = inertia, c.p = cognitive, c.g = social))

test5 <- psoptim(par = rep(NA, np), fn = negfwrap, opt = 1, lower = -100, upper = 100,
                 control = list(maxit = niter, s = nswarm, p = 1, vectorize = TRUE,
                                hybrid = TRUE))

test6 <- psoptim(par = rep(NA, np), fn = negfwrap, opt = 1, lower = -100, upper = 100,
                 control = list(maxit = niter, s = nswarm, p = 1, vectorize = TRUE,
                                w = c(1, 0.3)))

temp <- pso(niter, nswarm, inertia, cognitive, social, inits, nbhd[[2]], fwrap, opt=1)

temp2 <- pso(niter, nswarm, 1, cognitive, social, inits, nbhd[[2]], fwrap, opt=n,
            tune = TRUE, style = "adaptive", rate = rates[2])


inits <- matrix(runif(nswarm*np, -100, 100), ncol = nswarm)
nswarm <- 20
np <- 20
niter <- 500
nrep <- 50
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496
sig0 <- 1
rates <- c(0.1, .3, 0.5, 0.7)
betas <- c(1, 2, 4)
alphas <- c(0.1, 0.2, 0.4)*niter
dfs <- c(1, 3, 5, Inf)

nbhd <- list()
nbhd[[1]] <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)}) ## ring-1
nbhd[[2]] <- matrix(1:nswarm, ncol=nswarm, nrow=nswarm) ## global
nbhd[[3]] <- sapply(1:nswarm, function(x){return( (x + -3:3 - 1)%%nswarm + 1)}) ## ring-3
time <- 0:niter
nbhdnames <- c("ring-1", "global", "ring-3")

normvec <- function(x){
  sqrt(mean(x^2))
}

set.seed(342)

initfun <- function(n, nswarm, np){
  inits <- switch(n,
         matrix(runif(nswarm*np, 50, 100), ncol = nswarm),
         matrix(runif(nswarm*np, 50, 100), ncol = nswarm),
         matrix(runif(nswarm*np, 50, 100), ncol = nswarm),
         matrix(runif(nswarm*np, 15, 30), ncol = nswarm),
         matrix(runif(nswarm*np, 2.56, 5.12), ncol = nswarm),
         matrix(runif(nswarm*np, -500, -250), ncol = nswarm),
         matrix(runif(nswarm*np, 300, 600), ncol = nswarm),
         matrix(runif(nswarm*np, 16, 32), ncol = nswarm))
  return(inits)
}



## PSO
psoout <- NULL
for(n in c(1:8)[-c(3,6)]){
  for(m in 1:3){
    for(l in 1:nrep){
      print(c(n,m,l))
      inits <- initfun(n, nswarm, np)
      print("PSO")
      temp <- pso(niter, nswarm, inertia, cognitive, social, inits, nbhd[[m]], fwrap, opt=n)
      psoout <- rbind(psoout, data.frame(obj = n, logpost = temp[["maxes"]],
                                         argnorm = apply(temp[["argmaxes"]], 2, normvec),
                                         time = time, type = "PSO", nbhd = nbhdnames[m], rep = l,
                                         inertias = rep(0, niter+1)))
      print("BBPSO-MC")
      temp <- bbpso(niter, nswarm, 0, 1, inits, nbhd[[m]], fwrap, Inf, FALSE, 0, opt = n)
      psoout <- rbind(psoout, data.frame(obj = n, logpost= temp[["maxes"]],
                                         argnorm = apply(temp[["argmaxes"]], 2, normvec),
                                         time = time, type = "BBPSO-MC",
                                         nbhd = nbhdnames[m], rep = l,
                                         inertias = rep(0, niter+1)))
      print("BBPSOxp-MC")
      temp <- bbpso(niter, nswarm, 0, 1, inits, nbhd[[m]], fwrap, Inf, FALSE, .5, opt = n)
      psoout <- rbind(psoout, data.frame(obj = n, logpost= temp[["maxes"]],
                                         argnorm = apply(temp[["argmaxes"]], 2, normvec),
                                         time = time, type = "BBPSOxp-MC",
                                         nbhd = nbhdnames[m], rep = l,
                                         inertias = rep(0, niter+1)))
      for(i in 1:length(dfs)){
        for(j in 1:length(rates)){
          print(paste("AT-BBPSO-MC", dfs[i], rates[j], sep="-"))
          temp <- bbpso(niter, nswarm, sig0, rates[j], inits, nbhd[[m]], fwrap, dfs[i],
                        TRUE, 0,  opt = n)
          psoout <- rbind(psoout, data.frame(obj = n, logpost= temp[["maxes"]],
                                             argnorm = apply(temp[["argmaxes"]], 2, normvec),
                                             time = time,
                                             type = paste("AT-BBPSO-MC", dfs[i], rates[j], sep="-"),
                                             nbhd = nbhdnames[m], rep = l,
                                             inertias = rep(0, niter+1)))
          print(paste("AT-BBPSOxp-MC", dfs[i], rates[j], sep="-"))
          temp <- bbpso(niter, nswarm, sig0, rates[j], inits, nbhd[[m]], fwrap, dfs[i],
                        TRUE, .5,  opt = n)
          psoout <- rbind(psoout,
                          data.frame(obj = n, logpost= temp[["maxes"]],
                                     argnorm = apply(temp[["argmaxes"]], 2, normvec),
                                     time = time,
                                     type = paste("AT-BBPSOxp-MC", dfs[i], rates[j], sep="-"),
                                     nbhd = nbhdnames[m], rep = l,
                                     inertias = rep(0, niter+1)))
        }
      }
      for(i in 1:length(alphas)){
        for(j in 1:length(betas)){
          print(paste("DI-PSO", alphas[i], betas[j], sep="-"))
          temp <- pso(niter, nswarm, inertia, cognitive, social, inits, nbhd[[m]], fwrap, opt=n,
                      tune = TRUE, style = "deterministic", alpha = alphas[i], beta = betas[j])
          psoout <- rbind(psoout, data.frame(obj = n, logpost = temp[["maxes"]],
                                             argnorm = apply(temp[["argmaxes"]], 2, normvec),
                                             time = time,
                                             type = paste("DI-PSO", alphas[i], betas[j], sep="-"),
                                             nbhd = nbhdnames[m], rep = l,
                                             inertias = c(0,temp$inertias)))
        }
      }
      for(i in 1:length(rates)){
        print(paste("AT-PSO", rates[i], sep="-"))
        temp <- pso(niter, nswarm, 1, cognitive, social, inits, nbhd[[m]], fwrap, opt=n,
                    tune = TRUE, style = "adaptive", rate = rates[i])
        psoout <- rbind(psoout, data.frame(obj = n, logpost = temp[["maxes"]],
                                           argnorm = apply(temp[["argmaxes"]], 2, normvec),
                                           time = time,
                                           type = paste("AT-PSO", rates[i], sep="-"),
                                           nbhd = nbhdnames[m], rep = l,
                                           inertias = c(0,temp$inertias)))
      }
    }
  }
  write.csv(psoout, file = "psoout.csv", row.names=FALSE)
}
