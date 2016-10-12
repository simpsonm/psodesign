library(mnormt)
library(sp)  ## needed for points in polygon function

exch <- function(ncand, obj, poly, nnbor, npoints, start = "rand", ...){
  ndim <- ceiling(sqrt(ncand))
  ncand <- ndim^2
  mins <- apply(poly, 2, min)
  maxes <- apply(poly, 2, max)
  by.x <- (maxes[1] - mins[1])/ndim
  by.y <- (maxes[2] - mins[2])/ndim
  cand2 <- expand.grid(1:ndim, 1:ndim)
  cand <- expand.grid(seq(from = mins[1], by = by.x, length.out = ndim),
                      seq(from = mins[2], by = by.y, length.out = ndim))
  ## remove all candidate points not actually in a polygon
  checks <- apply(cand, 1, function(x, poly) {
    point.in.polygon(x[1], x[2], poly[,1], poly[,2])}, poly = poly)
  grid <- as.matrix(cand[checks==1,])
  grid2 <- as.matrix(cand2[checks==1,])
  ngrid <- nrow(grid)
  ## make distance matrix
  DST <- matrix(0, ngrid, ngrid)
  for(i in 2:ngrid){
    for(j in 1:(i-1)){
      DST[i,j] <- sum(abs(grid2[i,] - grid2[j,]))
      DST[j,i] <- DST[i,j]
    }
  }
  ## start, if necessary
  if(start == "rand"){
    IDstart <- sample(1:ngrid, npoints)
    Dstart <- grid[IDstart,]
  }
  Dbest <- Dstart
  IDbest <- IDstart
  Qbest <- obj(c(Dstart), ...)
  ## main while loop
  improve <- TRUE
  D <- Dstart
  ID <- IDstart
  vals <- Qbest
  count <- 1
  while(improve){
    for(i in 1:npoints){
      IDnb <- order(DST[ID[i],])[-1]
      Ddist <- DST[ID[i],IDnb]
      Ddist <- Ddist[which(Ddist <= nnbor)]
      for(j in 1:length(Ddist)){
        Dtemp <- D
        Dtemp[i,] <- grid[IDnb[j],]
        Qtemp <- obj(c(Dtemp), ...)
        if(Qtemp < Qbest){
          Qbest <- Qtemp
          D <- Dtemp
          ID[i] <- IDnb[j]
          Dbest <- D
        }
      }
    }
    vals <- c(vals, Qbest)
    count <- count + 1
    improve <- vals[count] < vals[count - 1]
  }
  return(list(par = Dbest, value = Qbest, vals = vals, count = count))
}

ga <- function(niter, nbatch, nswarm, nchrome, nrun, mutvar, mutrate, lower, upper, obj, ...){
  npar <- nchrome*nrun
  x <- replicate(nswarm, runif(npar, lower, upper))
  vals <- apply(x, 2, obj, ...)
  valsout <- matrix(0, ncol = nswarm, nrow = niter*nbatch+1)
  valsout[1,] <- vals
  for(bat in 1:nbatch){
    for(iter in 1:niter){
      ## crossover
      xranks <- rank(vals)
      parents <- replicate(nswarm, sample(1:nswarm, 2, prob = 1/xranks))
      genes <- matrix(rbinom(nrun*nswarm, 1, 0.5), nrow=nrun)[rep(1:nrun, nchrome),]
      xcross <- genes*x[,parents[1,]] + (1-genes)*x[,parents[2,]]
      ## mutation
      mutate <- matrix(rbinom(nrun*nswarm, 1, exp(-mutrate*iter)), nrow=nrun)[rep(1:nrun, nchrome),]
      z <- (x - lower)/(upper - lower)
      d <- log(z/(1-z)) + mutvar*exp(-mutrate*iter)*runif(npar*nswarm, -0.5, 0.5)
      u <- lower + (upper - lower)/(1 + exp(-d))
      xmut <- x*(1-mutate) + u*mutate
      crossvals <- apply(xcross, 2, obj, ...)
      mutvals <- apply(xmut, 2, obj, ...)
      fullvals <- c(vals, crossvals, mutvals)
      fullorder <- order(fullvals)
      fullx <- cbind(x, xcross, xmut)
      vals <- fullvals[fullorder[1:nswarm]]
      x <- fullx[,fullorder[1:nswarm]]
      valsout[(bat-1)*niter + iter + 1,] <- vals
    }
  }
  out <- list(par = x[,1], value = vals[1], parpop = x, valpops = valsout)
  return(out)
}

sbbpso <- function(niter, nswarm, nnbor, sig, pcut=0.5, CF=FALSE, AT=FALSE, rate=0.3, df=1,
                   ccc = 0.1, obj, ...){
  npar <- length(lower) ## dimension of search space
  if(nswarm < 1){
    nswarm <- 40  ## automatic setting of nswarm
  }
  sigs <- rep(sig, nswarm + 1)
  if(AT){
    logsig <- log(sig)
  }
  ## initialize positions, velocities, and nbhds
  x <- replicate(nswarm, runif(npar, lower, upper)) 
  inform <- matrix(replicate(nswarm, sample(1:nswarm, nnbor, replace = TRUE)), nrow = nnbor)
  nbhd <- lapply(1:nswarm, function(x) unique(c(which(inform == x, TRUE)[,2],x)))
  ## initialize pbest, gbest, and nbest stuff
  pbest <- x
  nbest <- x
  pbestval <- apply(x, 2, obj, ...)
  gbest <- which.min(pbestval)
  nbestval <- pbestval
  gbestvals <- rep(0, niter)
  gbestvals[1] <- pbestval[gbest]
  gbests <- matrix(0, ncol = niter + 1, nrow = npar)
  gbests[,1] <- pbest[,gbest]
  better <- rep(0, nswarm)
  for(iter in 1:niter){
    partid <- sample(1:nswarm, nswarm)  ## random update order
    for(id in 1:nswarm){
      idx <- partid[id]
      ## update nbhd best
      nminidx <- nbhd[[idx]][which.min(pbestval[nbhd[[idx]]])]
      nbestval[idx] <- pbestval[nminidx]
      nbest[,idx] <- pbest[,nminidx]
      ## bbpso update
      if(pbestval[idx] > nbestval[idx]){
        ## if personal best is worse than nbhd best, compute SDs
        if(CF){
          ## coordinate-free way -> constant SD across coords
          sds <- rep(sqrt(crossprod(pbest[,idx] - gbest))/npar, npar)*sig
        } else {
          ## standard BBPSO SD -> different SD across coords
          sds <- abs(pbest[,idx] - gbest)*sig
        }
        sd0 <- which(sds == 0)
        sd1 <- which(sds > 0)
        if(length(sd1) > 0){
          ## for SDs which are positive, do the usual bbpso draw
          temp <- rmtfixed(1, (pbest[sd1,idx] + nbest[sd1,idx])/2, diag(sds[sd1]), df)
          u <- runif(length(sd1))
          ## include a 1 - pcut chance of just going to the personal best coordinate
          x[sd1,idx] <- ifelse(u > pcut, temp, pbest[sd1,idx])
        }
        if(length(sd0) > 0){
          ## for SDs which are negative, do a mutation from the entire swarm
          id0s <- sample(2:nswarm, 3)
          id0s[id0s <= idx] <- id0s[id0s <= idx] - 1
          x[sd0,idx] <- pbest[sd0,id0s[1]] + (pbest[sd0,id0s[2]] - pbest[sd0,id0s[3]])/2
        }
      } else {
        ## if personal best is the same as nbhd best, do a mutation from the entire swarm
        id0s <- sample(2:nswarm, 3)
        id0s[id0s <= idx] <- id0s[id0s <= idx] - 1
        x[,idx] <- pbest[,id0s[1]] + (pbest[,id0s[2]] - pbest[,id0s[3]])/2
      }
      ## compute new value and update pbests if its an improvement
      newval <- obj(x[,idx], ...)
      if(newval < pbestval[idx]){
        pbestval[idx] <- newval
        pbest[,idx] <- x[,idx]
        better[idx] <- 1
      } else {
        better[idx] <- 0
      }      
    }
    ## update global best
    gbest <- which.min(pbestval)
    gbestvalue <- min(pbestval)
    gbestvals[iter + 1] <- gbestvalue
    gbests[,iter + 1] <- pbest[,gbest]
    ## if no gbest improvements, create new nbhds
    if(gbestvalue >= gbestvals[iter]){
      nbhd <- lapply(1:nswarm, function(x) unique(c(which(inform == x, TRUE)[,2],x)))
    }
    ## update sigma if AT
    if(AT){
      logsig <- logsig + ccc*(mean(better) - rate)/2
      sig <- exp(logsig)
      sigs[iter + 1] <- sig
    }
  }
  outlist <- list(par = pbest[,gbest], value = gbestvalue, pos = x, 
                  values = gbestvals, pars = gbests, sigs = sigs)
  return(outlist)
}

spso <- function(niter, nswarm, nnbor, inertia, cognitive, social, obj, lower, upper,
                 style = "CI", CF = FALSE, alpha = .2*niter, beta = 2, rate = 0.3, ccc = 0.1, ...){
  npar <- length(lower) ## dimension of search space
  if(nswarm < 1){
    nswarm <- 10 + 2*sqrt(npar)  ## automatic setting of nswarm
  }
  ## initialization of non-constant inertia stuff
  DI <- style == "DI"
  AT <- style == "AT"
  if(DI){
    inertia <- 1/(1 + (1/alpha)^beta)
  } else if(AT){
    loginertia <- log(inertia)
  }
  inertias <- rep(inertia, niter + 1)  
  ## initialize positions, velocities, and nbhds
  x <- replicate(nswarm, runif(npar, lower, upper)) 
  v <- replicate(nswarm, runif(npar, lower, upper)) - x  
  inform <- matrix(replicate(nswarm, sample(1:nswarm, nnbor, replace = TRUE)), nrow = nnbor)
  nbhd <- lapply(1:nswarm, function(x) unique(c(which(inform == x, TRUE)[,2],x)))
  ## initialize pbest, gbest, and nbest stuff
  pbest <- x
  nbest <- x
  pbestval <- apply(x, 2, obj, ...)
  gbest <- which.min(pbestval)
  nbestval <- pbestval
  gbestvals <- rep(0, niter)
  gbestvals[1] <- pbestval[gbest]
  gbests <- matrix(0, ncol = niter + 1, nrow = npar)
  gbests[,1] <- pbest[,gbest]
  better <- rep(0, nswarm)
  for(iter in 1:niter){
    partid <- sample(1:nswarm, nswarm)  ## random update order
    for(id in 1:nswarm){
      idx <- partid[id]
      ## update nbhd best
      nminidx <- nbhd[[idx]][which.min(pbestval[nbhd[[idx]]])]
      nbestval[idx] <- pbestval[nminidx]
      nbest[,idx] <- pbest[,nminidx]
      ## update velocity
      if(CF){
        ## coordinate free update of velocity
        ## center of gravity
        G <- x[,idx] + cognitive*(pbest[,idx] - x[,idx])/3
        if(nbestval[idx] < pbestval[idx]){
          ## only include this part of the update if not at nbhd best
          G <- G + social*(nbest[,idx] - x[,idx])/3
        }
        Gxdist <- sqrt(sum((x[,idx] - G)^2))
        ## draw from sphere
        if(npar == 1){
          xprime <- G + runif(1, -Gxdist, Gxdist)
        } else if(npar == 2){
          angle <- runif(1, 0, 2*pi)
          radius <- runif(1, 0, Gxdist)
          xprime <- radius*c(cos(angle), sin(angle)) + G
        } else {
          angles <- c(runif(npar - 2, 0, pi), runif(1, 0, 2*pi))
          sines <- sin(angles)
          cosines <- cos(angles)
          radius <- runif(1, 0, Gxdist)
          trigs <- c(cosines[1], cumprod(sines)[-(npar - 1)]*cosines[-1], cumprod(sines)[npar - 1])
          xprime <- G + radius*trigs
        }
        v[,idx] <- inertia*v[,idx] + xprime - x[,idx]
      } else {
        ## standard update of velocity
        ## update velocity
        v[,idx] <- inertia*v[,idx] + runif(npar, 0, cognitive)*(pbest[,idx] - x[,idx])
        if(nbestval[idx] < pbestval[idx]){
          ## only include this part of the update if not at nbhd best
          v[,idx] <- v[,idx] + runif(npar, 0, social)*(nbest[,idx] - x[,idx])
        }
      }
      ## udpate position
      x[,idx] <- x[,idx] + v[,idx]
      ## if position invalid, put it on the boundary and set velocity to 0
      toosmall <- which(x[,idx] < lower)
      if(length(toosmall)>0){
        x[toosmall, idx] <- lower[toosmall]
        v[toosmall, idx] <- -0.5*v[toosmall, idx] 
      }
      toolarge <- which(x[,idx] > upper)
      if(length(toolarge)>0){
        x[toolarge, idx] <- upper[toolarge]
        v[toolarge, idx] <- -0.5*v[toolarge, idx] 
      }
      ## compute new value and update pbests if its an improvement
      newval <- obj(x[,idx], ...)
      if(newval < pbestval[idx]){
        pbestval[idx] <- newval
        pbest[,idx] <- x[,idx]
        better[idx] <- 1
      } else {
        better[idx] <- 0
      }      
    }
    ## update global best
    gbest <- which.min(pbestval)
    gbestvalue <- min(pbestval)
    gbestvals[iter + 1] <- gbestvalue
    gbests[,iter + 1] <- pbest[,gbest]
    ## if no gbest improvements, create new nbhds
    if(gbestvalue >= gbestvals[iter]){
      nbhd <- lapply(1:nswarm, function(x) unique(c(which(inform == x, TRUE)[,2],x)))
    }
    ## update inertia if DI or AT
    if(DI){
      inertia <- 1/(1 + ((iter + 1)/alpha)^beta)
      inertias[iter + 1] <- inertia
    } else if(AT){
      loginertia <- loginertia + ccc*(mean(better) - rate)
      inertia <- exp(loginertia)
      inertias[iter + 1] <- inertia
    }
  }
  outlist <- list(par = pbest[,gbest], value = gbestvalue, pos = x, vel = v,
                  values = gbestvals, pars = gbests, inertias = inertias)
  return(outlist)
}

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
  maxdiff <- max(apply(x, 1, max) -  apply(x, 1, min))
  v <- matrix(runif(np*nswarm, -maxdiff/2, maxdiff/2), ncol=nswarm)
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
        lam <- lam + ccc*(acc - rate)
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
        sds <- abs(pbest[,i] - gbest)*sig
        sd0 <- which(sds == 0)
        sd1 <- which(sds > 0)
        if(length(sd1) > 0){
          temp <- rmtfixed(1, (pbest[sd1,i] + nbest[sd1,i])/2, diag(sds[sd1]), df)
          u <- runif(length(sd1))
          x[sd1,i] <- ifelse(u > pcut, temp, pbest[sd1,i])
        }
        if(length(sd0) > 0){
          idxs <- sample(2:nswarm, 3)
          idxs[idxs <= i] <- idxs[idxs <= i] - 1
          x[sd0,i] <- pbest[sd0,idxs[1]] + (pbest[sd0,idxs[2]] - pbest[sd0,idxs[3]])/2
        }
      } else {
        idxs <- sample(2:nswarm, 3)
        idxs[idxs <= i] <- idxs[idxs <= i] - 1
        x[,i] <- pbest[,idxs[1]] + (pbest[,idxs[2]] - pbest[,idxs[3]])/2
      }
    }
    vals <- apply(x, 2, obj, ...)
    better <- vals > pbestval
    idx <- which(better)
    if(tune){
      acc <- mean(better[-gbestidx]) ## only use particles that used the sd
      lam <- lam + ccc*(acc - rate)/2
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
