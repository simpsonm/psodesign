## stochastic versions
library(mnormt)

## rmt in mnormt library has a bug! use this instead
rmtfixed <- function(n = 1, mean = rep(0, d), S, df = Inf, sqrt = NULL){
  sqrt.S <- if (is.null(sqrt)) 
              chol(S)
            else sqrt
  d <- if (is.matrix(sqrt.S)) 
         ncol(sqrt.S)
       else 1
  x <- if (df == Inf) 
         1
       else rchisq(n, df)/df
  z <- rmnorm(n, rep(0, d), sqrt = sqrt.S)
  mean <- outer(rep(1, n), as.vector(matrix(mean, d)))
  drop(mean + z/sqrt(x))
}

pso <- function(niter, nswarm, inertia, cognitive, social, init, nbhd, obj, ...,
                tune = FALSE, style = "adaptive",
                rate = 0.3, ccc = 0.1, alpha = 1, beta = 0.1){
  np <- nrow(init)
  nnbor <- nrow(nbhd)
  x <- init
  v <- matrix(rep(0, np*nswarm), ncol=nswarm)
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
  for(iter in 1:niter){
    inertias[iter] <- inertia
    rand1 <- matrix(runif(nswarm), nrow = np, ncol = nswarm, byrow = FALSE)
    rand2 <- matrix(runif(nswarm), nrow = np, ncol = nswarm, byrow = FALSE)
    v <- inertia*v + cognitive*rand1*(pbest - x) + social*rand2*(nbest - x)
    x <- x + v
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
  outlist <- list(argmax = pbest[,gbest],  max = gbestvalue, pos = x, maxes = gbestvals, argmaxes = gbests, inertias = inertias)
  return(outlist)
}

stochpso <- function(niter, nswarm, nsub, inertia, cognitive, social, init, nbhd, obj, datlist){
  np <- nrow(init)
  nnbor <- nrow(nbhd)
  datlisttemp <- datlist
  ndat <- nrow(datlist$dat)
  x <- init
  v <- matrix(rep(0, np*nswarm), ncol=nswarm)
  pbest <- x
  gbest <- x
  nbest <- x
  pbestval <- rep(0, nswarm)
  vals <- pbestval
  if(nsub < ndat){
    for(i in 1:nswarm){
      datlisttemp$dat <- datlist$dat[sample(ndat, nsub),]
      pbestval[i] <- obj(x[,i], datlisttemp)
    }
  } else {
    pbestval <- apply(x, 2, obj, datlisttemp)
  }
  nbestval <- pbestval
  gbestvals <- rep(0, niter)
  gbestvals[1] <- max(pbestval)
  gbests <- matrix(0, ncol = niter + 1, nrow = np)
  gbests[,1] <- pbest[,which.max(pbestval)]
  for(iter in 1:niter){
    rand1 <- matrix(runif(nswarm), nrow = np, ncol = nswarm, byrow = FALSE)
    rand2 <- matrix(runif(nswarm), nrow = np, ncol = nswarm, byrow = FALSE)
    v <- inertia*v + cognitive*rand1*(pbest - x) + social*rand2*(nbest - x)
    x <- x + v
    if(nsub < ndat){
      for(i in 1:nswarm){
        datlisttemp$dat <- datlist$dat[sample(ndat, nsub),]
        vals[i] <- obj(x[,i], datlisttemp)
      }
    } else {
      vals <- apply(x, 2, obj, datlisttemp)
    }
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
  }
  outlist <- list(argmax = pbest[,gbest],  max = gbestvalue, pos = x, maxes = gbestvals, argmaxes = gbests)
  return(outlist)
}


bbpso <- function(niter, nswarm, sig, rate, init, nbhd, obj, df, tune, pcut, ..., ccc = 0.1){
  np <- nrow(init)
  nnbor <- nrow(nbhd)
  x <- init
  pbest <- x
  nbest <- x
  pbestval <- apply(x, 2, obj, ...)
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
        temp <- rmtfixed(1, (pbest[,i] + nbest[,i])/2, diag(sds), df)
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
    vals <- apply(x, 2, obj, ...)
    better <- vals > pbestval
    idx <- which(better)
    if(tune){
      acc <- mean(better[-gbestidx]) ## only use particles that used the sd
      if( acc > rate){
        lam <- lam + ccc
      } else {
        lam <-  lam - ccc
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


stochbbpso <- function(niter, nswarm, nsub, sig, rate, init, nbhd, obj, df, tune, pcut, datlist){
  np <- nrow(init)
  nnbor <- nrow(nbhd)
  datlisttemp <- datlist
  ndat <- nrow(datlist$dat)
  x <- init
  pbest <- x
  nbest <- x
  pbestval <- rep(0, nswarm)
  vals <- pbestval
  if(nsub < ndat){
    for(i in 1:nswarm){
      datlisttemp$dat <- datlist$dat[sample(ndat, nsub),]
      pbestval[i] <- obj(x[,i], datlisttemp)
    }
  } else {
    pbestval <- apply(x, 2, obj, datlisttemp)
  }
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
        temp <- rmtfixed(1, (pbest[,i] + nbest[,i])/2, diag(sds), df)
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
    if(nsub < ndat){
      for(i in 1:nswarm){
        datlisttemp$dat <- datlist$dat[sample(ndat, nsub),]
        vals[i] <- obj(x[,i], datlisttemp)
      }
    } else {
      vals <- apply(x, 2, obj, datlisttemp)
    }
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
