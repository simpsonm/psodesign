library(mnormt)

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

bbpso <- function(niter, nswarm, sig, rate, init, nbhd, obj, df, tune, pcut, ...){
  np <- nrow(init)
  nnbor <- nrow(nbhd)
  x <- init
  pbest <- x
  nbest <- x
  pbestval <- obj(x, ...)
  nbestval <- pbestval
  gbestval <- max(pbestval)
  gbestvals <- rep(0, niter + 1)
  gbestvals[1] <- max(pbestval)
  gbests <- matrix(0, ncol = niter + 1, nrow = np)
  gbestidx <- which.max(pbestval)
  gbests[,1] <- pbest[,gbestidx]
  gbest <- pbest[,gbestidx]
  lams <- rep(0, niter)
  lam <- log(sig)
  for(iter in 1:niter){
    maxidxs <- apply(matrix(pbestval[nbhd], ncol = nswarm), 2, which.max)
    idxs <- nbhd[(1:nswarm - 1)*nnbor + maxidxs]
    nbest <- pbest[,idxs]
    nbestval <- pbestval[idxs]
    for(i in 1:nswarm){
      if(pbestval[i] < nbestval[i]){
        sds <- abs(pbest[,i] - gbest)
        sds[sds == 0] <- 0.001
        temp <- rmt(1, (pbest[,i] + nbest[,i])/2, diag(sds), df)
        u <- runif(np)
        x[,i] <- ifelse(u > pcut[1], temp, pbest[,i])
      } else {
        idxs <- sample(2:nswarm, 3)
        idxs[idxs <= i] <- idxs[idxs <= i] - 1
        temp <- pbest[,idxs[1]] + (pbest[,idxs[2]] - pbest[,idxs[3]])/2
        u <- runif(np)
        x[,i] <- ifelse(u > pcut[2], temp, gbest)
      }
    }
    vals <- obj(x, ...)
    better <- vals > pbestval
    idx <- which(better)
    if(tune){
      acc <- mean(better[-gbestidx]) ## only use particles that used the sd
      if( acc > rate){
        lam <- lam + 0.01
      } else {
        lam <-  lam - 0.01
      }
      sig <- exp(lam)
    }
    pbest[,idx] <- x[,idx]
    pbestval[idx] <- vals[idx]
    if(max(pbestval) > gbestval){
      gbestidx <- which.max(pbestval)
      gbest <- pbest[,gbestidx]
      gbestval <- max(pbestval)
    }
    gbests[,iter + 1] <- gbest
    gbestvals[iter + 1] <- gbestval
    lams[iter] <- lam
  }
  outlist <- list(argmax = gbest,  max = gbestval, pos = x, maxes = gbestvals, argmaxes = gbests, lambdas = lams)
  return(outlist) 
}
