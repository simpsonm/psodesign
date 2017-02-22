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
## inertia:      inertia parameter
## cognitive:    cognitive correction factor
## social:       social correction factor
## style:        "CI" : constant inertia (i.e., standard PSO)
##               "DI" : deterministic inertia
##               "AT" : adaptively tune
## CF:           if TRUE use the coordinate-free version of the algorithm
## rate:         desired improvement rate (only used if style=="AT")
##               should be strictly between 0 and 1, (0.3 - 0.5 works well)
## ccc:          tuning parameter c (only used if style=="AT")
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
## returns a list with the named elements:
## par :         vector of best location found
## value :       scalar value of objective function at par
## pos :         position of the entire swarm at the end of the algorithm
##               each column is a distinct particle
## vel :         velocity of the entire swarm at the end of the algorithm
##               each column is the velocity of a distinct particle
## values :      time series of best objective function values found every iteration
## pars :        time series of best locations found every iteration
##               each column corresponds to a separate iteration
## inertias :    time series of the value of omega (inertia) every iteration (for AT/DI variants)
atpso <- function(niter, obj, lower, upper, nswarm = 40, nnbor = 3, inertia = 0.7298,
                  cognitive = 1.496, social = 1.496, style = "CI", CF = FALSE,
                  alpha = .2*niter, beta = 2, rate = 0.3, ccc = 0.1,
                  init = NULL, boundaryfun = NULL, ...){
  npar <- length(lower) ## dimension of search space
  if(!is.vector(lower) | !is.vector(upper) | length(lower) != length(upper))
    stop("lower and upper must be vectors of the same length")  
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
  if(is.null(init)){
    x <- replicate(nswarm, runif(npar, lower, upper))
  } else {
    x <- init
  }
  nobfun <- is.null(boundaryfun)   ## check to see if boundary function is supplied.
  v <- replicate(nswarm, runif(npar, lower, upper)) - x
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
      ## update velocity
      if(CF){
        ## coordinate free update of velocity
        ## center of gravity
        if(nbestval[idx] < pbestval[idx]){
          G <- x[,idx] + cognitive*(pbest[,idx] - x[,idx])/3 + social*(nbest[,idx] - x[,idx])/3
        } else {
          ## only use cognitive component if nbhd best = p best
          G <- x[,idx] + cognitive*(pbest[,idx] - x[,idx])/2
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
      if(nobfun){
        ## if no boundary function, use the given bounding box
        ## if position invalid, put it on the boundary and set velocity to -0.5*velocity
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
      } else {
        ## if boundary function exists, use it to move outliers to the boundary
        newx <- boundaryfun(x[,idx], ...)
        moved <- which(x[,idx] != newx)
        if(length(moved)>0){
          x[,idx] <- newx
          v[moved, idx] <- -0.5*v[moved, idx]
        }
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
    if(nnbor < nswarm & gbestvalue >= gbestvals[iter]){
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
