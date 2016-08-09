pso <- function(niter, nswarm, inertia, cognitive, social, init, nbhd, obj, ...){
  np <- nrow(init)
  nnbor <- nrow(nbhd)
  x <- init
  v <- matrix(rep(0, np*nswarm), ncol=nswarm)
  pbest <- x
  gbest <- x
  nbest <- x
  pbestval <- obj(x, ...)
  nbestval <- pbestval
  gbestvals <- rep(0, niter)
  gbestvals[1] <- max(pbestval)
  gbests <- matrix(0, ncol = niter + 1, nrow = np)
  gbests[,1] <- pbest[,which.max(pbestval)]
  for(iter in 1:niter){
    vals <- obj(x, ...)
    idx <- which(vals > pbestval)   ## if new position is better
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
    rand1 <- matrix(runif(nswarm), nrow = np, ncol = nswarm, byrow = FALSE)
    rand2 <- matrix(runif(nswarm), nrow = np, ncol = nswarm, byrow = FALSE)
    v <- inertia*v + cognitive*rand1*(pbest - x) + social*rand2*(nbest - x)
    x <- x + v
  }
  outlist <- list(argmax = pbest[,gbest],  max = gbestvalue, pos = x, maxes = gbestvals, argmaxes = gbests)
  return(outlist)
}

bbpso <- function(niter, nswarm, sig, eta, rate, maxfail, init, nbhd, obj, df, tune, ...){
  np <- nrow(init)
  nnbor <- nrow(nbhd)
  x <- init
  pbest <- x
  nbest <- x
  pbestval <- obj(x, ...)
  nbestval <- pbestval
  gbest <- which.max(pbestval)
  gbestvals <- rep(0, niter + 1)
  lams <- rep(0, niter)
  gbestvals[1] <- max(pbestval)
  gbests <- matrix(0, ncol = niter + 1, nrow = np)
  gbests[,1] <- pbest[,which.max(pbestval)]
  lam <- log(sig)
  nfail <- rep(0, nswarm)
  mutate <- nfail
  for(iter in 1:niter){
    maxidxs <- apply(matrix(pbestval[nbhd], ncol = nswarm), 2, which.max)
    idxs <- nbhd[(1:nswarm - 1)*nnbor + maxidxs]
    nbest <- pbest[,idxs]
    nbestval <- pbestval[idxs]
    for(i in 1:nswarm){
      if(pbestval[i] < nbestval[i] & nfail[i] < maxfail){
        cholS <- sig*diag(abs(pbest[,i] - pbest[,gbest]))
        x[,i] <- (pbest[,i] + nbest[,i])/2 + cholS%*%rt(np, df=df)
        mutate[i] <- 0
      } else {
        x[,i] <- pbest[,i]*(1 + rcauchy(np, eta))
        mutate[i] <- 1
      }
    }
    vals <- obj(x, ...)
    better <- vals > pbestval
    idx <- which(better)
    nfail[idx] <- 0
    nfail[-idx] <- (nfail[-idx] + 1)*(1 - mutate[-idx]) + 0
    if(tune){
      nsig <- sum(1-mutate)
      if(nsig > 9){
        acc <- sum(better*(1-mutate)) / nsig
        if( acc > rate){
          lam <- lam + 0.01
        } else {
          lam <-  lam - 0.01
        }
        sig <- exp(lam)
      }
    }
    pbest[,idx] <- x[,idx]
    pbestval[idx] <- vals[idx]   
    gbest <- which.max(pbestval)
    gbests[,iter + 1] <- pbest[,gbest]
    gbestvalue <- pbestval[gbest]
    gbestvals[iter + 1] <- gbestvalue
    lams[iter] <- lam
    ##print(round(c(iter, acc, exp(lam), gbestvalue),3))
  }
  outlist <- list(argmax = pbest[,gbest],  max = gbestvalue, pos = x, maxes = gbestvals, argmaxes = gbests, lambdas = lams)
  return(outlist) 
}
