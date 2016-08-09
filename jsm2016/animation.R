library(ggplot2)

pso <- function(niter, nswarm, inertia, cognitive, social, init, nbhd, obj, ...,
                tune = FALSE, style = "adaptive",
                rate = 0.3, ccc = 0.1, alpha = 1, beta = 0.1){
  np <- nrow(init)
  nnbor <- nrow(nbhd)
  x <- init
  v <- matrix(runif(np*nswarm, -5, 5), ncol=nswarm)
  pbest <- x
  gbest <- x
  nbest <- x
  pbestval <- apply(x, 2, obj, ...)
  nbestval <- pbestval
  gbestvals <- rep(0, niter)
  gbestvals[1] <- max(pbestval)
  gbests <- matrix(0, ncol = niter + 1, nrow = np)
  gbests[,1] <- pbest[,which.max(pbestval)]
  lam <- log(inertia)
  inertias <- rep(0, niter)
  if(tune == TRUE & style == "deterministic"){
    inertia <- 1/(1 + (1/alpha)^beta)
  }
  history <- array(0, c(np, nswarm, niter + 1))
  history[,,1] <- x
  for(iter in 1:niter){
    inertias[iter] <- inertia
    rand1 <- matrix(runif(nswarm), nrow = np, ncol = nswarm, byrow = FALSE)
    rand2 <- matrix(runif(nswarm), nrow = np, ncol = nswarm, byrow = FALSE)
    v <- inertia*v + cognitive*rand1*(pbest - x) + social*rand2*(nbest - x)
    x <- x + v
    history[,,iter + 1] <- x
    vals <- apply(x, 2, obj, ...)
    better <- vals > pbestval
    idx <- which(better)   ## if new position is better
    pbest[,idx] <- x[,idx]         ## update best x,
    pbestval[idx] <- vals[idx]      ## and best value
    gbest <- which.max(pbestval)    ## global best
    gbestvalue <- max(pbestval)
    gbestvals[iter + 1] <- gbestvalue
    gbests[,iter + 1] <- pbest[,gbest]
    maxidxs <- apply(matrix(pbestval[nbhd], ncol = nswarm), 2, which.max)  ## update nbhd bests
    idxs <- nbhd[(1:nswarm - 1)*nnbor + maxidxs]
    nbest <- pbest[,idxs]
    nbestval <- pbestval[idxs]
    if(tune){
      if(style == "adaptive"){
        acc <- mean(better)
        if(acc > rate){
          lam <- lam + ccc
        } else {
          lam <- lam - ccc
        }
        inertia <- exp(lam)
      } else if(style == "deterministic"){
        inertia <- 1/(1 + ((iter + 1)/alpha)^beta)
      } else {
        cat("\n")
        cat("error: style must be \"adaptive\" or \"deterministic\"\n")
        return(NULL)
      }
    }
  }
  outlist <- list(argmax = pbest[,gbest],  max = gbestvalue, pos = x, maxes = gbestvals, argmaxes = gbests, inertias = inertias, history = history)
  return(outlist)
}


niter <- 25
nswarm <- 50
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496
nbhd <- sapply(1:nswarm, function(x){return( (x + -3:3 - 1)%%nswarm + 1)}) ## ring-3

set.seed(3429)

alpha <- 1.4
beta <- 2.5

z <- rbeta(100, alpha, beta)

betaloglike <- function(pars, z){
  if(is.null(nrow(pars))){
    out <- sum(dbeta(z, exp(pars[1]), exp(pars[2]), log = TRUE))
  } else {
    out <- rep(0, nrow(pars))
    for(i in 1:nrow(pars)){
      out[i] <- sum(dbeta(z, exp(pars[i,1]), exp(pars[i,2]), log = TRUE))
    }
  }
  return(out)
}

### normal MLE (bfgs)
mleout <- optim(c(1,1), betaloglike, z = z, method = "BFGS",
                control=list(fnscale=-1))


### pso
init <- matrix(runif(nswarm*2, -10, 10), ncol = nswarm)
psoout <- pso(niter, nswarm, inertia, cognitive, social, init, nbhd, betaloglike, z = z)

psodf <- data.frame(iter=0, part=1:nswarm, alpha=psoout$history[1,,1], beta=psoout$history[2,,1])
for(i in 1:niter){
  psodf <- rbind(psodf, data.frame(iter=i, part=1:nswarm, alpha=psoout$history[1,,i], beta=psoout$history[2,,i]))
}

ht <- 4
wd <- 4
for(i in 0:niter){
  animplot <- ggplot(data=subset(psodf, iter==i), aes(x=alpha, y=beta)) + geom_point() +
    xlim(c(-10,10)) + ylim(c(-10,10)) +
    labs(x = expression(log(alpha)), y = expression(log(beta))) +
    geom_point(aes(x=psoout$argmax[1], y=psoout$argmax[2]), color = "red")
  ggsave(paste("animplot", i, ".png", sep=""), animplot, width = wd, height = ht)
}


