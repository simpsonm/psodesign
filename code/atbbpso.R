library(mnormt)


## bare-bones PSO
## minimizes some objective function obj(x, ...) in x.
##
## warning: very little input validation
##
## niter:        number of iterations
## obj:          objective function to be minimized. has the form:
##               obj(x, ...) where x is a vector in R^D, and ... is other arguments
## lower:        vector of lower bounds for each dimension of x
## upper:        vector of upper bounds for each dimension of x
##               if obj's domain is complicated, lower & upper should bound the domain
## nswarm:       size of the swarm
## nnbor:        number of particles informed in stochastic star topology
##               if nnbor >= nswarm, use the global topology instead
## sig:          scale parameter. only the initial value if AT==TRUE
## pcut:         set = 0.5 for BBPSOxp, set = 0 for BBPSO
## CF:           if TRUE use the coordinate-free version of the algorithm
## AT:           if TRUE adaptively tune sigma
## rate:         desired improvement rate (only used if AT==TRUE)
##               should be strictly between 0 and 1, (0.3 - 0.5 works well)
## df:           degrees of freedom of t kernel. should be >= 1
##               df large (>30) is approximately a normal kernel and works well when AT=FALSE
##               df=1 works well when AT=TRUE
## ccc:          tuning parameter c (only used if AT==TRUE)
##               controls how fast ATBBPSO adapts ( ccc = 0.1 works well)
## init:         optional initial position of the swarm
##               should be a npar x nswarm matrix where npar = length(x)
##               (each column corresponds to the position of a particle)
## boundaryfun:  optional function for moving particle to the nearest point 
##               on the boundary if it is outside of it.
##               should take a vector x and other objective function arguments
##               and return a vector the same length as x
##               e.g.: x <- boundaryfun(x, ...)
## ... :         optional arguments passed to obj() and boundaryfun()
##
## returns a list with the named elements:
## par :         vector of best location found
## value :       scalar value of objective function at par
## pos :         position of the entire swarm at the end of the algorithm
##               each column is a distinct particle
## values :      time series of best objective function values found every iteration
## pars :        time series of best locations found every iteration
##               each column corresponds to a separate iteration
## sigs :        time series of the value of sigma every iteration (for AT variants)
atbbpso <- function(niter, obj, lower, upper, nswarm = 40, nnbor=3, sig=1, pcut=0.5,
                    CF=FALSE, AT=FALSE, rate=0.3, df=1,
                    ccc = 0.1, init = NULL, boundaryfun = NULL, ...){
  npar <- length(lower) ## dimension of search space
  if(!is.vector(lower) | !is.vector(upper) | length(lower) != length(upper))
    stop("lower and upper must be vectors of the same length")
  sigs <- rep(sig, nswarm + 1)
  if(AT){
    logsig <- log(sig)
  }
  ## initialize positions, velocities, and nbhds
  if(is.null(init)){
    x <- replicate(nswarm, runif(npar, lower, upper))
  } else {
    x <- init
  }
  nobfun <- is.null(boundaryfun)   ## check to see if boundary function is supplied.
  if(nnbor < nswarm){
    inform <- matrix(replicate(nswarm, sample(1:nswarm, nnbor, replace = TRUE)), nrow = nnbor)
    nbhd <- lapply(1:nswarm, function(x) unique(c(which(inform == x, TRUE)[,2],x)))
  } else {
    nbhd <- lapply(1:nswarm, function(x) 1:nswarm)
  }
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
        if(length(sd0) > 0){
          ## for SDs which are negative, do a mutation from the entire swarm
          id0s <- sample(2:nswarm, 3)
          id0s[id0s <= idx] <- id0s[id0s <= idx] - 1
          x[sd0,idx] <- pbest[sd0,id0s[1]] + (pbest[sd0,id0s[2]] - pbest[sd0,id0s[3]])/2
        }
        if(length(sd1) > 0){
          ## for SDs which are positive, do the usual bbpso draw
          temp <- rmt(1, (pbest[sd1,idx] + nbest[sd1,idx])/2, diag(sds[sd1]), df)
          u <- runif(length(sd1))
          ## include a 1 - pcut chance of just going to the personal best coordinate
          x[sd1,idx] <- ifelse(u > pcut, temp, pbest[sd1,idx])
        }
      } else {
        ## if personal best is the same as nbhd best, do a mutation from the entire swarm
        id0s <- sample(2:nswarm, 3)
        id0s[id0s <= idx] <- id0s[id0s <= idx] - 1
        x[,idx] <- pbest[,id0s[1]] + (pbest[,id0s[2]] - pbest[,id0s[3]])/2
      }
      ## move particle to the boundary of the space, if necessary
      if(nobfun){
        ## if no boundary function, use the given bounding box
        toosmall <- which(x[,idx] < lower)
        if(length(toosmall)>0){
          x[toosmall, idx] <- lower[toosmall]
        }
        toolarge <- which(x[,idx] > upper)
        if(length(toolarge)>0){
          x[toolarge, idx] <- upper[toolarge]
        }
      } else {
        ## if boundary function exists, use it to move outliers to the boundary
        x[,idx] <- boundaryfun(x[,idx], ...)
      }
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
    ## if no gbest improvements, create new nbhds (assuming not gbest nbhd)
    if(nnbor < nswarm & gbestvalue >= gbestvals[iter]){
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
