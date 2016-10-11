source("../psofun.R")
source("testfuns.R")
library(geoR)
library(rgdal)

ncand <- 10000
coshape <- readOGR(dsn="../kriging/shape", layer="05000")

poly <- coshape@polygons[[which(coshape$NAME == "Harris")]]@Polygons[[1]]@coords













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
gniter <- 100
nbatch <- 10
nchrome <- 1
nrun <- 10
mutvar1 <- 3
mutrate1 <- 0.1
mutvar2 <- 1
mutrate2 <- 1
glower <- rep(-100, nchrome*nrun)
gupper <- rep( 100, nchrome*nrun)

fnum <- 2
test1 <- spso(niter, nswarm, nnbor, inertia, cognitive, social, negfwrap, lower,
              upper, style = "CI", opt = fnum)
test2 <- spso(niter, nswarm, nnbor, c(.2*niter, 2), cognitive, social, negfwrap, lower,
              upper, style = "DI", opt = fnum)
test3 <- spso(niter, nswarm, nnbor, c(1.2, 0.3, 0.1), cognitive, social, negfwrap, lower,
              upper, style = "AT", opt = fnum)
temp1 <- spso(niter, nswarm, nnbor, inertia, cognitive, social, negfwrap, lower,
              upper, style = "CI", CF = TRUE, opt = fnum)
temp2 <- spso(niter, nswarm, nnbor, c(.2*niter, 2), cognitive, social, negfwrap, lower,
              upper, style = "DI", CF = TRUE, opt = fnum)
temp3 <- spso(niter, nswarm, nnbor, c(1.2, 0.3, .1), cognitive, social, negfwrap, lower,
              upper, style = "AT", CF = TRUE, opt = fnum)
tt1 <- sbbpso(niter, nswarm, nnbor, 1, obj = negfwrap, opt = fnum)
tt2 <- sbbpso(niter, nswarm, nnbor, 1, CF = TRUE, obj = negfwrap, opt = fnum)
tt3 <- sbbpso(niter, nswarm, nnbor, 1, AT=TRUE, obj = negfwrap, opt = fnum)
tt4 <- sbbpso(niter, nswarm, nnbor, 1, AT=TRUE, CF = TRUE, obj = negfwrap, opt = fnum)
ga11 <- ga(gniter, nbatch, nswarm, nchrome, nrun, mutvar1, mutrate1, glower, gupper, negfwrap,
           opt = fnum)
ga12 <- ga(gniter, nbatch, nswarm, nchrome, nrun, mutvar1, mutrate2, glower, gupper, negfwrap,
           opt = fnum)
ga21 <- ga(gniter, nbatch, nswarm, nchrome, nrun, mutvar2, mutrate1, glower, gupper, negfwrap,
           opt = fnum)
ga22 <- ga(gniter, nbatch, nswarm, nchrome, nrun, mutvar2, mutrate2, glower, gupper, negfwrap,
           opt = fnum)
outs <- matrix(c(test1$value, test2$value, test3$value, temp1$value, temp2$value, temp3$value,
                 tt1$value, tt2$value, tt3$value, tt4$value,
                 ga11$value, ga12$value, ga21$value, ga22$value), ncol = 1)
rownames(outs) <- c(paste("2007", c("CI", "DI", "AT"), sep="-"),
                    paste("2011", c("CI", "DI", "AT"), sep="-"),
                    paste("BBPSO", c("-", "-CF", "-AT", "-CF-AT"), sep=""),
                    paste("ga", c(11, 12, 21, 22), sep="-"))
round(outs, 5)

n <- 10
p <- rnorm(n)
g <- rnorm(n) + p
w <- 0.9
sig <- w*tcrossprod(p - g) + (1-w)*diag((p-g)^2)
diag(1/abs(p-g))%*%sig%*%diag(1/abs(p-g))
det(sig)
try(cholsig <- chol(sig))




str(gatest)
nbhd <- list()
nbhd[[1]] <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)}) ## ring-1
nbhd[[2]] <- matrix(1:nswarm, ncol=nswarm, nrow=nswarm) ## global
nbhd[[3]] <- sapply(1:nswarm, function(x){return( (x + -3:3 - 1)%%nswarm + 1)}) ## ring-3

inits <- matrix(runif(nswarm*np, -100, 100), ncol = nswarm)

old <- pso(niter, nswarm, 1, cognitive, social, inits, nbhd[[2]], fwrap, opt=1,
           tune = TRUE, style = "adaptive", rate = 0.3)



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
