library(MCMCpack)
library(mnormt)
library(Matrix)
#source("popfun.R")

dmtcholprec <- function (x, mn, Rprec, df, LOG = TRUE){
  d <- length(x)
  X <- x - mn
  RX <- Rprec%*%X
  out <- -(df + d)/2*log(1 + crossprod(RX)/df)
  if(LOG)
    return(out)
  else
    return(exp(out))
}

rwishchol <- function (v, S){
  if (!is.matrix(S)) 
    S <- matrix(S)
  if (nrow(S) != ncol(S)) {
    stop(message = "S not square in rwish().\n")
  }
  if (v < nrow(S)) {
    stop(message = "v is less than the dimension of S in rwish().\n")
  }
  p <- nrow(S)
  CC <- chol(S)
  Z <- matrix(0, p, p)
  diag(Z) <- sqrt(rchisq(p, v:(v - p + 1)))
  if (p > 1) {
    pseq <- 1:(p - 1)
    Z[rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] <-
      rnorm(p * (p - 1)/2)
  }
  R <- Z%*%CC
  return(t(R))
}

## fixed multivariate t random number generation
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

## single move rw within gibbs for poisson data model, iid random effects
poisiidrwgibbs <- function(niter, init, datlist, tune = TRUE, logrwsds = NULL,
                           rwc = 0.1, rwtarget = 0.44){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  if(is.null(logrwsds)){
    logrwsds <- rep(2.38, nbeta + ndelta)
  }
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  sig2a <- datlist$sig2a
  sig2b <- datlist$sig2b
  z <- dat[, 1]
  xmat <- dat[, 1 + 1:nbeta]
  smat <- dat[, 1 + nbeta + 1:ndelta]
  beta <-  as.matrix(init[1:nbeta])
  xbeta <- xmat%*%beta
  delta <- as.matrix(init[nbeta + 1:ndelta])
  sdelta <- smat%*%delta
  y <- xbeta + sdelta
  npars <- length(init)
  draws <- matrix(0, ncol = npars, nrow = niter + 1)
  colnames(draws) <- c(paste("beta", 1:nbeta, sep = "."),
                       paste("delta", 1:ndelta, sep = "."), "sig2")
  draws[1,] <- init
  accs <- matrix(0, nrow = niter, ncol = ndelta + nbeta)
  sig2atilde <- ndelta/2 + sig2a
  datpart <- sum(z*y     - exp(y))
  mhratios <- matrix(0, nrow = niter, ncol = nbeta + ndelta)
  for(iter in 1:niter){
    sig2 <- 1/rgamma(1, sig2atilde, crossprod(delta)/2 + sig2b)
    betaprop <- beta
    sdelta <- smat%*%delta
    for(i in 1:nbeta){
      betaprop[i] <- beta[i] + exp(logrwsds[i])*rnorm(1)
      xbetaprop <- xmat%*%betaprop
      yprop <- xbetaprop + sdelta
      datpartprop <- sum(z*yprop - exp(yprop))
      lanum <-   datpartprop - (betaprop[i] - betamn)^2 / (2*betavar)
      ladenom <- datpart     - (beta[i]     - betamn)^2 / (2*betavar)
      u <- runif(1)
      mhratios[iter, i] <- exp(min(0, lanum - ladenom))
      if(u < mhratios[iter, i]){
        beta[i] <- betaprop[i]
        datpart <- datpartprop
        accs[iter,i] <- 1
        y <- yprop
      } else {
        betaprop[i] <- beta[i]
      }
    }
    deltaprop <- delta
    xbeta <- xmat%*%beta
    for(i in 1:ndelta){
      deltaprop[i] <- delta[i] + exp(logrwsds[nbeta + i])*rnorm(1)
      sdeltaprop <- smat%*%deltaprop
      yprop <- xbeta + sdeltaprop
      datpartprop <- sum(z*yprop - exp(yprop))
      u <- runif(1)
      lanum <-   datpartprop - deltaprop[i]^2 / (2*sig2)
      ladenom <- datpart     - delta[i]^2     / (2*sig2)
      mhratios[iter, nbeta + i] <- exp(min(0, lanum - ladenom))
      if(u < mhratios[iter, nbeta + i]){
        delta[i] <- deltaprop[i]
        datpart <- datpartprop
        accs[iter,nbeta + i] <- 1
        y <- yprop
      } else {
        deltaprop[i] <- delta[i]
      }
    }
    draws[iter + 1,] <- c(beta, delta, sig2)
    if(tune){
      logrwsds <- logrwsds + rwc*(mhratios[iter,] - rwtarget)/2
    }
  }
  out <- list(draws = draws, accs = accs, mhratios = mhratios, logrwsds = logrwsds)
  return(out)
}

## single move rw within gibbs for poisson data model, iid random effects
lnormiidrwgibbs <- function(niter, init, datlist, tune = TRUE, logrwsds = NULL,
                            rwc = 0.1, rwtarget = 0.44){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  if(is.null(logrwsds)){
    logrwsds <- rep(2.38, nbeta + ndelta)
  }
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  aphi <- datlist$aphi
  bphi <- datlist$bphi
  sig2a <- datlist$sig2a
  sig2b <- datlist$sig2b
  z <- dat[, 1]
  lz <- log(z)
  ndat <- length(z)
  xmat <- dat[, 1 + 1:nbeta]
  smat <- dat[, 1 + nbeta + 1:ndelta]
  beta <-  as.matrix(init[1:nbeta])
  xbeta <- xmat%*%beta
  delta <- as.matrix(init[nbeta + 1:ndelta])
  sdelta <- smat%*%delta
  y <- xbeta + sdelta
  npars <- length(init)
  draws <- matrix(0, ncol = npars, nrow = niter + 1)
  colnames(draws) <- c(paste("beta", 1:nbeta, sep = "."),
                       paste("delta", 1:ndelta, sep = "."), "sig2", "phi2")
  draws[1,] <- init
  accs <- matrix(0, nrow = niter, ncol = ndelta + nbeta)
  sig2atilde <- ndelta/2 + sig2a
  phi2atilde <- ndat/2 + aphi
  datpart <- crossprod(lz - y)/2
  mhratios <- matrix(0, nrow = niter, ncol = nbeta + ndelta)
  for(iter in 1:niter){
    sig2 <- 1/rgamma(1, sig2atilde, crossprod(delta)/2 + sig2b)
    phi2 <- 1/rgamma(1, phi2atilde, datpart + bphi)
    betaprop <- beta
    sdelta <- smat%*%delta
    for(i in 1:nbeta){
      betaprop[i] <- beta[i] + exp(logrwsds[i])*rnorm(1)
      xbetaprop <- xmat%*%betaprop
      yprop <- xbetaprop + sdelta
      datpartprop <- crossprod(lz - yprop)/2
      lanum <-   - datpartprop/phi2 - (betaprop[i] - betamn)^2 / (2*betavar)
      ladenom <- - datpart/phi2     - (beta[i]     - betamn)^2 / (2*betavar)
      mhratios[iter,i] <- exp(min(0, lanum - ladenom))
      u <- runif(1)
      if(u < mhratios[iter,i]){
        beta[i] <- betaprop[i]
        datpart <- datpartprop
        accs[iter,i] <- 1
        y <- yprop
      } else {
        betaprop[i] <- beta[i]
      }
    }
    deltaprop <- delta
    xbeta <- xmat%*%beta
    for(i in 1:ndelta){
      deltaprop[i] <- delta[i] + exp(logrwsds[nbeta + i])*rnorm(1)
      sdeltaprop <- smat%*%deltaprop
      yprop <- xbeta + sdeltaprop
      datpartprop <- crossprod(lz - yprop)/2
      u <- runif(1)
      lanum <-   - datpartprop/phi2 - deltaprop[i]^2 / (2*sig2)
      ladenom <- - datpart/phi2     - delta[i]^2     / (2*sig2)
      mhratios[iter,nbeta + i] <- exp(min(0, lanum - ladenom))
      if(u < mhratios[iter,nbeta + i]){
        delta[i] <- deltaprop[i]
        datpart <- datpartprop
        accs[iter,nbeta + i] <- 1
        y <- yprop
      } else {
        deltaprop[i] <- delta[i]
      }
    }
    draws[iter + 1,] <- c(beta, delta, sig2, phi2)
    if(tune){
      logrwsds <- logrwsds + rwc*(mhratios[iter,] - rwtarget)/2
    }
  }
  out <- list(draws = draws, accs = accs, mhratios = mhratios, logrwsds = logrwsds)
  return(out)
}

## block rw within gibbs for poisson data model, iid random effects
poisiidblockrwgibbs <- function(niter, init, datlist, cholhat, tune = TRUE, logrwsd = NULL,
                                rwc = 0.1, rwtarget = 0.234){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  if(is.null(logrwsd)){
    logrwsd = log(2.38/sqrt(nbeta + ndelta))
  }
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  sig2a <- datlist$sig2a
  sig2b <- datlist$sig2b
  z <- dat[, 1]
  xmat <- as.matrix(dat[, 1 + 1:nbeta])
  smat <- as.matrix(dat[, 1 + nbeta + 1:ndelta])
  beta <-  as.matrix(init[1:nbeta])
  xbeta <- xmat%*%beta
  delta <- as.matrix(init[nbeta + 1:ndelta])
  sdelta <- smat%*%delta
  y <- xbeta + sdelta
  betadelta <- c(beta, delta)
  npars <- length(init)
  draws <- matrix(0, ncol = npars, nrow = niter + 1)
  colnames(draws) <- c(paste("beta", 1:nbeta, sep = "."),
                       paste("delta", 1:ndelta, sep = "."), "sig2")
  draws[1,] <- init
  accs <- rep(0, niter)
  sighat <- tcrossprod(cholhat)
  sig2atilde <- ndelta/2 + sig2a
  delpart <- crossprod(delta)/2
  datbetpart <- sum(z*y - exp(y)) - crossprod(beta - betamn)/(2*betavar)
  mhratios <- rep(0, niter)
  for(iter in 1:niter){
    sig2 <- 1/rgamma(1, sig2atilde, crossprod(delta)/2 + sig2b)
    betadeltaprop <- betadelta + cholhat%*%rnorm(nbeta + ndelta, 0, exp(logrwsd))
    betaprop <- betadeltaprop[1:nbeta]
    deltaprop <- betadeltaprop[nbeta + 1:ndelta]
    yprop <- xmat%*%betaprop + smat%*%deltaprop
    delpartprop <- crossprod(deltaprop)/2
    datbetpartprop <- sum(z*yprop - exp(yprop)) - crossprod(betaprop - betamn)/(2*betavar)
    lanum <-   - delpartprop/sig2 + datbetpartprop
    ladenom <- - delpart/sig2     + datbetpart
    u <- runif(1)
    mhratios[iter] <- exp(min(0, lanum - ladenom))
    if(u < mhratios[iter]){
      accs[iter] <- 1
      betadelta <- betadeltaprop
      delta <- betadeltaprop[nbeta + 1:ndelta]
      delpart <- delpartprop
      datbetpart <- datbetpartprop
      y <- yprop
    }
    draws[iter + 1,] <- c(betadelta, sig2)
    if(tune){
      logrwsd <- logrwsd + rwc*(mhratios[iter] - rwtarget)/2
    }    
  }
  out <- list(draws = draws, accs = accs, mhratios = mhratios, logrwsd = logrwsd)
  return(out)
}

## block rw within gibbs for lognormal data model, iid random effects
lnormiidblockrwgibbs <- function(niter, init, datlist, cholhat, tune = TRUE, logrwsd = NULL,
                                rwc = 0.1, rwtarget = 0.234){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  if(is.null(logrwsd)){
    logrwsd = log(2.38/sqrt(nbeta + ndelta))
  }
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  sig2a <- datlist$sig2a
  sig2b <- datlist$sig2b
  aphi <- datlist$aphi
  bphi <- datlist$bphi
  z <- dat[, 1]
  lz <- log(z)
  ndat <- length(z)
  xmat <- as.matrix(dat[, 1 + 1:nbeta])
  smat <- as.matrix(dat[, 1 + nbeta + 1:ndelta])
  beta <-  as.matrix(init[1:nbeta])
  xbeta <- xmat%*%beta
  delta <- as.matrix(init[nbeta + 1:ndelta])
  sdelta <- smat%*%delta
  y <- xbeta + sdelta
  betadelta <- c(beta, delta)
  npars <- length(init)
  draws <- matrix(0, ncol = npars, nrow = niter + 1)
  colnames(draws) <- c(paste("beta", 1:nbeta, sep = "."),
                       paste("delta", 1:ndelta, sep = "."), "sig2", "phi2")
  draws[1,] <- init
  accs <- rep(0, niter)
  sighat <- tcrossprod(cholhat)
  sig2atilde <- ndelta/2 + sig2a
  phi2atilde <- ndat/2 + aphi
  delpart <- crossprod(delta)/2
  betpart <- crossprod(beta - betamn)/(2*betavar)
  datpart <- crossprod(lz - y)/2
  mhratios <- rep(0, niter)
  for(iter in 1:niter){
    sig2 <- 1/rgamma(1, sig2atilde, delpart + sig2b)
    phi2 <- 1/rgamma(1, phi2atilde, datpart + bphi)
    betadeltaprop <- betadelta + cholhat%*%rnorm(nbeta + ndelta, 0, exp(logrwsd))
    betaprop <- betadeltaprop[1:nbeta]
    deltaprop <- betadeltaprop[nbeta + 1:ndelta]
    yprop <- xmat%*%betaprop + smat%*%deltaprop
    delpartprop <- crossprod(deltaprop)/2
    betpartprop <- crossprod(betaprop - betamn)/(2*betavar)
    datpartprop <- crossprod(lz - yprop)/2
    lanum <-   - delpartprop/sig2 - betpartprop - datpartprop/phi2
    ladenom <- - delpart/sig2     - betpart     - datpart/phi2
    u <- runif(1)
    mhratios[iter] <- exp(min(0, lanum - ladenom))
    if(u < mhratios[iter]){
      accs[iter] <- 1
      betadelta <- betadeltaprop
      delta <- betadeltaprop[nbeta + 1:ndelta]
      delpart <- delpartprop
      betpart <- betpartprop
      datpart <- datpartprop
      y <- yprop
    }
    draws[iter + 1,] <- c(betadelta, sig2, phi2)
    if(tune){
      logrwsd <- logrwsd + rwc*(mhratios[iter] - rwtarget)/2
    }    
  }
  out <- list(draws = draws, accs = accs, mhratios = mhratios, logrwsd = logrwsd)
  return(out)
}

## single move rw within gibbs for poisson data model, fully correlated random effects
poisfullrwgibbs <- function(niter, init, datlist, tune = TRUE, logrwsds = NULL, 
                            rwc = 0.1, rwtarget = 0.44){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  omdf <- datlist$omdf
  ominvscale <- datlist$ominvscale
  if(is.null(logrwsds)){
    logrwsds <- rep(2.38, nbeta + ndelta)
  }
  z <- dat[, 1]
  xmat <- dat[, 1 + 1:nbeta]
  smat <- dat[, 1 + nbeta + 1:ndelta]
  beta <-  as.matrix(init[1:nbeta])
  xbeta <- xmat%*%beta
  delta <- as.matrix(init[nbeta + 1:ndelta])
  sdelta <- smat%*%delta
  y <- xbeta + sdelta
  npars <- length(init)
  draws <- matrix(0, ncol = npars, nrow = niter + 1)
  ellnums <- apply(which(lower.tri(diag(ndelta), TRUE), TRUE), 1, paste, collapse = ".")
  colnames(draws) <- c(paste("beta", 1:nbeta, sep = "."),
                       paste("delta", 1:ndelta, sep = "."),
                       paste("ell", ellnums, sep = "."))
  draws[1,] <- init
  accs <- matrix(0, nrow = niter, ncol = ndelta + nbeta)
  datpart <- sum(z*y - exp(y))
  mhratios <- matrix(0, nrow = niter, ncol = nbeta + ndelta)
  for(iter in 1:niter){
    L <- rwishchol(omdf + 1, chol2inv(chol(tcrossprod(delta) + ominvscale)))
    ell <- L[lower.tri(L, TRUE)]
    betaprop <- beta
    sdelta <- smat%*%delta
    for(i in 1:nbeta){
      betaprop[i] <- beta[i] + exp(logrwsds[i])*rnorm(1)
      xbetaprop <- xmat%*%betaprop
      yprop <- xbetaprop + sdelta
      datpartprop <- sum(z*yprop - exp(yprop))
      lanum <-   datpartprop - (betaprop[i] - betamn)^2 / (2*betavar)
      ladenom <- datpart     - (beta[i]     - betamn)^2 / (2*betavar)
      u <- runif(1)
      mhratios[iter, i] <- exp(min(0, lanum - ladenom))
      if(u < mhratios[iter, i]){
        beta[i] <- betaprop[i]
        datpart <- datpartprop
        accs[iter,i] <- 1
        y <- yprop
      } else {
        betaprop[i] <- beta[i]
      }
    }
    deltaprop <- delta
    xbeta <- xmat%*%beta
    for(i in 1:ndelta){
      deltaprop[i] <- delta[i] + exp(logrwsds[nbeta + i])*rnorm(1)
      sdeltaprop <- smat%*%deltaprop
      yprop <- xbeta + sdeltaprop
      datpartprop <- sum(z*yprop - exp(yprop))
      u <- runif(1)
      lanum <-   datpartprop - tcrossprod(crossprod(deltaprop, L))/2
      ladenom <- datpart     - tcrossprod(crossprod(delta, L))/2
      mhratios[iter, nbeta + i] <- exp(min(0, lanum - ladenom))
      if(u < mhratios[iter, nbeta + i]){
        delta[i] <- deltaprop[i]
        datpart <- datpartprop
        accs[iter,nbeta + i] <- 1
        y <- yprop
      } else {
        deltaprop[i] <- delta[i]
      }
    }
    draws[iter + 1,] <- c(beta, delta, ell)
    if(tune){
      logrwsds <- logrwsds + rwc*(mhratios[iter,] - rwtarget)/2
    }
  }
  out <- list(draws = draws, accs = accs, mhratios = mhratios, logrwsds = logrwsds)
  return(out)
}

## single move rw within gibbs for poisson data model, full random effects
lnormfullrwgibbs <- function(niter, init, datlist, tune = TRUE, logrwsds = NULL,
                           rwc = 0.1, rwtarget = 0.44){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  aphi <- datlist$aphi
  bphi <- datlist$bphi
  omdf <- datlist$omdf
  ominvscale <- datlist$ominvscale
  if(is.null(logrwsds)){
    logrwsds <- rep(2.38, nbeta + ndelta)
  }
  z <- dat[, 1]
  lz <- log(z)
  ndat <- length(z)
  xmat <- dat[, 1 + 1:nbeta]
  smat <- dat[, 1 + nbeta + 1:ndelta]
  beta <-  as.matrix(init[1:nbeta])
  xbeta <- xmat%*%beta
  delta <- as.matrix(init[nbeta + 1:ndelta])
  sdelta <- smat%*%delta
  y <- xbeta + sdelta
  npars <- length(init)
  draws <- matrix(0, ncol = npars, nrow = niter + 1)
  ellnums <- apply(which(lower.tri(diag(ndelta), TRUE), TRUE), 1, paste, collapse = ".")
  colnames(draws) <- c(paste("beta", 1:nbeta, sep = "."),
                       paste("delta", 1:ndelta, sep = "."),
                       paste("ell", ellnums, sep = "."), "phi2")
  draws[1,] <- init
  accs <- matrix(0, nrow = niter, ncol = ndelta + nbeta)
  phi2atilde <- ndat/2 + aphi
  datpart <- crossprod(lz - y)/2
  mhratios <- matrix(0, nrow = niter, ncol = nbeta + ndelta)
  for(iter in 1:niter){
    L <- rwishchol(omdf + 1, chol2inv(chol(tcrossprod(delta) + ominvscale)))
    ell <- L[lower.tri(L, TRUE)]
    phi2 <- 1/rgamma(1, phi2atilde, datpart + bphi)
    betaprop <- beta
    sdelta <- smat%*%delta
    for(i in 1:nbeta){
      betaprop[i] <- beta[i] + exp(logrwsds[i])*rnorm(1)
      xbetaprop <- xmat%*%betaprop
      yprop <- xbetaprop + sdelta
      datpartprop <- crossprod(lz - yprop)/2
      lanum <-   - datpartprop/phi2 - (betaprop[i] - betamn)^2 / (2*betavar)
      ladenom <- - datpart/phi2     - (beta[i]     - betamn)^2 / (2*betavar)
      mhratios[iter,i] <- exp(min(0, lanum - ladenom))
      u <- runif(1)
      if(u < mhratios[iter,i]){
        beta[i] <- betaprop[i]
        y <- yprop
        datpart <- datpartprop
        accs[iter,i] <- 1
      } else {
        betaprop[i] <- beta[i]
      }
    }
    deltaprop <- delta
    xbeta <- xmat%*%beta
    for(i in 1:ndelta){
      deltaprop[i] <- delta[i] + exp(logrwsds[nbeta + i])*rnorm(1)
      sdeltaprop <- smat%*%deltaprop
      yprop <- xbeta + sdeltaprop
      u <- runif(1)
      datpartprop <- crossprod(lz - yprop)/2
      lanum <-   - datpartprop/phi2 - tcrossprod(crossprod(deltaprop, L))/2
      ladenom <- - datpart/phi2     - tcrossprod(crossprod(delta, L))/2
      mhratios[iter,nbeta + i] <- exp(min(0, lanum - ladenom))
      if(u < mhratios[iter,nbeta + i]){
        delta[i] <- deltaprop[i]
        y <- yprop
        datpart <- datpartprop
        accs[iter,nbeta + i] <- 1
      } else {
        deltaprop[i] <- delta[i]
      }
    }
    draws[iter + 1,] <- c(beta, delta, ell, phi2)
    if(tune){
      logrwsds <- logrwsds + rwc*(mhratios[iter,] - rwtarget)/2
    }
  }
  out <- list(draws = draws, accs = accs, mhratios = mhratios, logrwsds = logrwsds)
  return(out)
}


## block rw within gibbs for poisson data model, fully correlated random effects
poisfullblockrwgibbs <- function(niter, init, datlist, cholhat, tune = TRUE,
                                 logrwsd = NULL, rwc = 0.1, rwtarget = 0.234){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  if(is.null(logrwsd)){
    logrwsd <- log(2.38/sqrt(length(nbeta + ndelta)))
  }
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  omdf <- datlist$omdf
  ominvscale <- datlist$ominvscale
  z <- dat[, 1]
  xmat <- as.matrix(dat[, 1 + 1:nbeta])
  smat <- as.matrix(dat[, 1 + nbeta + 1:ndelta])
  beta <-  as.matrix(init[1:nbeta])
  xbeta <- xmat%*%beta
  delta <- as.matrix(init[nbeta + 1:ndelta])
  sdelta <- smat%*%delta
  y <- xbeta + sdelta
  betadelta <- c(beta, delta)
  npars <- length(init)
  draws <- matrix(0, ncol = npars, nrow = niter + 1)
  ellnums <- apply(which(lower.tri(diag(ndelta), TRUE), TRUE), 1, paste, collapse = ".")
  colnames(draws) <- c(paste("beta", 1:nbeta, sep = "."),
                       paste("delta", 1:ndelta, sep = "."),
                       paste("ell", ellnums, sep = "."))
  draws[1,] <- init
  accs <- rep(0, niter)
  sighat <- tcrossprod(cholhat)
  datbetpart <- sum(z*y - exp(y)) - crossprod(beta - betamn)/(2*betavar)
  mhratios <- rep(0, niter)
  for(iter in 1:niter){
    L <- rwishchol(omdf + 1, chol2inv(chol(tcrossprod(delta) + ominvscale)))    
    ell <- L[lower.tri(L, TRUE)]
    betadeltaprop <- betadelta + cholhat%*%rnorm(nbeta + ndelta, 0, exp(logrwsd))
    betaprop <- betadeltaprop[1:nbeta]
    deltaprop <- betadeltaprop[nbeta + 1:ndelta]
    yprop <- xmat%*%betaprop + smat%*%deltaprop
    datbetpartprop <- sum(z*yprop - exp(yprop)) - crossprod(betaprop - betamn)/(2*betavar)
    lanum <-   - tcrossprod(crossprod(deltaprop, L))/2 + datbetpartprop
    ladenom <- - tcrossprod(crossprod(delta, L))/2 + datbetpart
    u <- runif(1)
    mhratios[iter] <- exp(min(0, lanum - ladenom))
    if(u < mhratios[iter]){
      accs[iter] <- 1
      betadelta <- betadeltaprop
      beta <- betaprop
      delta <- deltaprop
      y <- yprop
      datbetpart <- datbetpartprop
    }
    draws[iter + 1,] <- c(betadelta, ell)
    if(tune){
      logrwsd <- logrwsd + rwc*(mhratios[iter] - rwtarget)/2
    }    
  }
  out <- list(draws = draws, accs = accs, mhratios = mhratios, logrwsd = logrwsd)
  return(out)
}

## block rw within gibbs for lognormal data model, full random effects
lnormfullblockrwgibbs <- function(niter, init, datlist, cholhat, tune = TRUE,
                                 logrwsd = NULL, rwc = 0.1, rwtarget = 0.234){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  if(is.null(logrwsd)){
    logrwsd <- log(2.38/sqrt(length(nbeta + ndelta)))
  }
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  omdf <- datlist$omdf
  ominvscale <- datlist$ominvscale
  aphi <- datlist$aphi
  bphi <- datlist$bphi
  z <- dat[, 1]
  lz <- log(z)
  ndat <- length(z)
  xmat <- as.matrix(dat[, 1 + 1:nbeta])
  smat <- as.matrix(dat[, 1 + nbeta + 1:ndelta])
  beta <-  as.matrix(init[1:nbeta])
  delta <- as.matrix(init[nbeta + 1:ndelta])
  y <- xmat%*%beta + smat%*%delta
  betadelta <- c(beta, delta)
  npars <- length(init)
  draws <- matrix(0, ncol = npars, nrow = niter + 1)
  ellnums <- apply(which(lower.tri(diag(ndelta), TRUE), TRUE), 1, paste, collapse = ".")
  colnames(draws) <- c(paste("beta", 1:nbeta, sep = "."),
                       paste("delta", 1:ndelta, sep = "."),
                       paste("ell", ellnums, sep = "."), "phi2")
  draws[1,] <- init
  accs <- rep(0, niter)
  sighat <- tcrossprod(cholhat)
  phi2atilde <- ndat/2 + aphi
  betpart <- crossprod(beta - betamn)/(2*betavar)
  datpart <- crossprod(lz - y)/2
  mhratios <- rep(0, niter)  
  for(iter in 1:niter){
    L <- rwishchol(omdf + 1, chol2inv(chol(tcrossprod(delta) + ominvscale)))    
    ell <- L[lower.tri(L, TRUE)]
    phi2 <- 1/rgamma(1, phi2atilde, datpart + bphi)
    betadeltaprop <- betadelta + cholhat%*%rnorm(nbeta + ndelta, 0, exp(logrwsd))
    betaprop <- betadeltaprop[1:nbeta]
    deltaprop <- betadeltaprop[nbeta + 1:ndelta]
    yprop <- xmat%*%betaprop + smat%*%deltaprop
    betpartprop <- crossprod(betaprop - betamn)/(2*betavar)
    datpartprop <- crossprod(lz - yprop)/2
    lanum <-   - tcrossprod(crossprod(deltaprop, L))/2 - betpartprop - datpartprop/phi2
    ladenom <- - tcrossprod(crossprod(delta, L))/2 - betpart - datpart/phi2
    u <- runif(1)
    mhratios[iter] <- exp(min(0, lanum - ladenom))
    if(u < mhratios[iter]){
      accs[iter] <- 1
      betadelta <- betadeltaprop
      beta <- betaprop
      delta <- deltaprop
      y <- yprop
      betpart <- betpartprop
      datpart <- datpartprop
    }
    draws[iter + 1,] <- c(betadelta, ell, phi2)
    if(tune){
      logrwsd <- logrwsd + rwc*(mhratios[iter] - rwtarget)/2
    }    
  }
  out <- list(draws = draws, accs = accs, mhratios = mhratios, logrwsd = logrwsd)
  return(out)
}


## full gibbs for lognormal data model, full random effects
lnormfullgibbs <- function(niter, init, datlist){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  omdf <- datlist$omdf
  ominvscale <- datlist$ominvscale
  aphi <- datlist$aphi
  bphi <- datlist$bphi
  z <- dat[, 1]
  lz <- log(z)
  ndat <- length(z)
  xmat <- as.matrix(dat[, 1 + 1:nbeta])
  smat <- as.matrix(dat[, 1 + nbeta + 1:ndelta])
  xsmat <- dat[, 1 + 1:(nbeta + ndelta)]
  beta <-  as.matrix(init[1:nbeta])
  delta <- as.matrix(init[nbeta + 1:ndelta])
  betadelta <- c(beta, delta)
  npars <- length(init)
  draws <- matrix(0, ncol = npars, nrow = niter + 1)
  ellnums <- apply(which(lower.tri(diag(ndelta), TRUE), TRUE), 1, paste, collapse = ".")
  colnames(draws) <- c(paste("beta", 1:nbeta, sep = "."),
                       paste("delta", 1:ndelta, sep = "."),
                       paste("ell", ellnums, sep = "."), "phi2")
  draws[1,] <- init
  phi2atilde <- ndat/2 + aphi
  xsxs <- crossprod(xsmat)
  Ebeta <- diag(1/betavar, nbeta)
  cbetadelta <- c(betamn/betavar, rep(0, ndelta)) 
  xslz <- crossprod(xsmat, lz)
  for(iter in 1:niter){
    L <- rwishchol(omdf + 1, chol2inv(chol(tcrossprod(delta) + ominvscale)))    
    ell <- L[lower.tri(L, TRUE)]
    omega <- tcrossprod(L)
    y <- xsmat%*%betadelta
    phi2 <- 1/rgamma(1, phi2atilde, crossprod(lz - y)/2 + bphi)
    omegabd <- xsxs/phi2 + bdiag(Ebeta, omega)
    sigmabd <- chol2inv(chol(omegabd))
    mubd <- sigmabd%*%(xslz/phi2 + cbetadelta)
    betadelta <- drop(mubd + crossprod(chol(sigmabd), rnorm(nbeta + ndelta)))
    beta <- betadelta[1:nbeta]
    delta <- betadelta[nbeta + 1:ndelta]
    draws[iter + 1,] <- c(betadelta, ell, phi2)
  }
  out <- list(draws = draws)
  return(out)
}

## full gibbs for lognormal data model, iid random effects
lnormiidgibbs <- function(niter, init, datlist){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  sig2a <- datlist$sig2a
  sig2b <- datlist$sig2b
  aphi <- datlist$aphi
  bphi <- datlist$bphi
  z <- dat[, 1]
  lz <- log(z)
  ndat <- length(z)
  xmat <- as.matrix(dat[, 1 + 1:nbeta])
  smat <- as.matrix(dat[, 1 + nbeta + 1:ndelta])
  xsmat <- dat[, 1 + 1:(nbeta + ndelta)]
  beta <-  as.matrix(init[1:nbeta])
  delta <- as.matrix(init[nbeta + 1:ndelta])
  betadelta <- c(beta, delta)
  npars <- length(init)
  draws <- matrix(0, ncol = npars, nrow = niter + 1)
  ellnums <- apply(which(lower.tri(diag(ndelta), TRUE), TRUE), 1, paste, collapse = ".")
  colnames(draws) <- c(paste("beta", 1:nbeta, sep = "."),
                       paste("delta", 1:ndelta, sep = "."),
                       "sig2", "phi2")
  draws[1,] <- init
  phi2atilde <- ndat/2 + aphi
  sig2atilde <- ndelta/2 + sig2a
  xsxs <- crossprod(xsmat)
  Ebeta <- diag(1/betavar, nbeta)
  cbetadelta <- c(betamn/betavar, rep(0, ndelta))
  xslz <- crossprod(xsmat, lz)
  for(iter in 1:niter){
    y <- xsmat%*%betadelta
    sig2 <- 1/rgamma(1, sig2atilde, crossprod(delta)/2 + sig2b)
    phi2 <- 1/rgamma(1, phi2atilde, crossprod(lz - y)/2 + bphi)
    omegabd <- xsxs/phi2 + bdiag(Ebeta, diag(1/sig2, ndelta))
    sigmabd <- chol2inv(chol(omegabd))
    mubd <- sigmabd%*%(xslz/phi2 + cbetadelta)
    betadelta <- drop(mubd + crossprod(chol(sigmabd), rnorm(nbeta + ndelta)))
    beta <- betadelta[1:nbeta]
    delta <- betadelta[nbeta + 1:ndelta]
    draws[iter + 1,] <- c(betadelta, sig2, phi2)
  }
  out <- list(draws = draws)
  return(out)
}


## MIH w/in gibbs for poisson data model, iid random effects
## uses conditional distribution from normal approximation
poisiidindgibbs <- function(niter, init, mu, df, lpbest, datlist, tune = TRUE){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  sig2a <- datlist$sig2a
  sig2b <- datlist$sig2b
  z <- dat[, 1]
  xmat <- dat[, 1 + 1:nbeta]
  smat <- dat[, 1 + nbeta + 1:ndelta]
  betadelta <- init[1:(nbeta + ndelta)]
  beta <- init[1:nbeta]
  delta <- init[1:ndelta + nbeta]
  npars <- length(init)
  draws <- matrix(0, ncol = npars, nrow = niter + 1)
  colnames(draws) <- c(paste("beta", 1:nbeta, sep = "."),
                       paste("delta", 1:ndelta, sep = "."), "sig2")
  lpbests <- rep(0, niter)
  draws[1,] <- init
  acc <- rep(0, niter)
  sig2atilde <- ndelta/2 + sig2a
  mubd <- mu[1:(nbeta + ndelta)]
  mulsig2 <- mu[nbeta + ndelta + 1]
  hess <- poislposthessiid(mu, datlist)
  SigmaBig <- chol2inv(chol(-hess))
  Sigma1 <- SigmaBig[1:(nbeta + ndelta), 1:(nbeta + ndelta)]
  Sigma12 <- matrix(drop(SigmaBig[1:(nbeta + ndelta), nbeta + ndelta + 1]), ncol = 1)
  Sigma22 <- drop(SigmaBig[nbeta + ndelta + 1, nbeta + ndelta + 1])
  Sigma <- Sigma1 - tcrossprod(Sigma12)/Sigma22
  Rsigma <- chol(Sigma)
  Rprec <- chol(chol2inv(Rsigma))
  for(iter in 1:niter){
    ##if(iter %% 50 == 0){
    ##  print(iter)
    ##}
    sig2 <- 1/rgamma(1, sig2atilde, crossprod(delta)/2 + sig2b)
    mn <- drop(mubd + Sigma12*(log(sig2) - mulsig2)/Sigma22)
    betadeltaprop <- rmtfixed(1, mn, Sigma, df, Rsigma)
    ljprop <- dmtcholprec(betadeltaprop, mn, Rprec, df, TRUE)
    ljold <-  dmtcholprec(betadelta,     mn, Rprec, df, TRUE)
    lpprop <- poislpostiid(c(betadeltaprop, log(sig2)), datlist)
    lpold <-  poislpostiid(c(betadelta,     log(sig2)), datlist)
    la <- lpprop - lpold + ljold - ljprop
    u <- runif(1)
    if(tune &lpbest < lpprop){
      lpbest <- lpprop
      mu <- c(betadeltaprop, log(sig2))
      mbd <- mu[1:(nbeta + ndelta)]
      hess <- poislposthessiid(mu, datlist)
      SigmaBig <- chol2inv(chol(-hess))
      Sigma1 <- SigmaBig[1:(nbeta + ndelta), 1:(nbeta + ndelta)]
      Sigma12 <- matrix(drop(SigmaBig[1:(nbeta + ndelta), nbeta + ndelta + 1]), ncol = 1)
      Sigma22 <- drop(SigmaBig[nbeta + ndelta + 1, nbeta + ndelta + 1])
      Sigma <- Sigma1 - tcrossprod(Sigma12)/Sigma22
      Rsigma <- chol(Sigma)
      Rprec <- chol(chol2inv(Rsigma))
    }
    if(tune & lpbest < lpold){
      lpbest <- lpold
      mu <- c(betadelta, log(sig2))
      mbd <- mu[1:(nbeta + ndelta)]
      hess <- poislposthessiid(mu, datlist)
      SigmaBig <- chol2inv(chol(-hess))
      Sigma1 <- SigmaBig[1:(nbeta + ndelta), 1:(nbeta + ndelta)]
      Sigma12 <- matrix(drop(SigmaBig[1:(nbeta + ndelta), nbeta + ndelta + 1]), ncol = 1)
      Sigma22 <- drop(SigmaBig[nbeta + ndelta + 1, nbeta + ndelta + 1])
      Sigma <- Sigma1 - tcrossprod(Sigma12)/Sigma22
      Rsigma <- chol(Sigma)
      Rprec <- chol(chol2inv(Rsigma))
    }
    if(log(u) < la){
      betadelta <- betadeltaprop
      beta <- betadelta[1:nbeta]
      delta <- betadelta[1:ndelta + nbeta]
      acc[iter] <- 1
    }
    lpbests[iter] <- lpbest
    par <- c(betadelta, sig2)
    draws[iter + 1,] <- par
  }
  return(list(draws = draws, acc = acc, state = par, lpbest = lpbest, lpbests = lpbests,
              mu = mu, Sigma = Sigma, Rsigma = Rsigma))
}

## MIH w/in gibbs for lognormal data model, iid random effects
## uses conditional distribution from normal approximation
lnormiidindgibbs <- function(niter, init, mu, df, lpbest, datlist, tune = TRUE){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  sig2a <- datlist$sig2a
  sig2b <- datlist$sig2b
  aphi <- datlist$aphi
  bphi <- datlist$bphi
  z <- dat[, 1]
  lz <- log(z)
  ndat <- length(z)
  xmat <- as.matrix(dat[, 1 + 1:nbeta])
  smat <- as.matrix(dat[, 1 + nbeta + 1:ndelta])
  betadelta <- init[1:(nbeta + ndelta)]
  beta <- drop(init[1:nbeta])
  delta <- init[1:ndelta + nbeta]
  npars <- length(init)
  draws <- matrix(0, ncol = npars, nrow = niter + 1)
  colnames(draws) <- c(paste("beta", 1:nbeta, sep = "."),
                       paste("delta", 1:ndelta, sep = "."), "sig2", "phi2")
  lpbests <- rep(0, niter)
  draws[1,] <- init
  acc <- rep(0, niter)
  sig2atilde <- ndelta/2 + sig2a
  phi2atilde <- ndat/2 + aphi
  mubd <- mu[1:(nbeta + ndelta)]
  beta <- mu[1:nbeta]
  delta <- mu[1:ndelta + nbeta]
  mulsig2 <- mu[nbeta + ndelta + 1]
  mulphi2 <- mu[nbeta + ndelta + 2]
  hess <- lnormlposthessiid(mu, datlist)
  SigmaBig <- chol2inv(chol(-hess))
  Sigma1 <- SigmaBig[1:(nbeta + ndelta), 1:(nbeta + ndelta)]
  Sigma12 <- SigmaBig[1:(nbeta + ndelta), nbeta + ndelta + 1:2]
  Sigma22 <- SigmaBig[nbeta + ndelta + 1:2, nbeta + ndelta + 1:2]
  Sigma <- Sigma1 - Sigma12%*%tcrossprod(chol2inv(chol(Sigma22)),Sigma12)
  Rsigma <- chol(Sigma)
  Rprec <- chol(chol2inv(Rsigma))
  for(iter in 1:niter){
    ##if(iter %% 50 == 0){
    ##  print(iter)
    ##}
    sig2 <- 1/rgamma(1, sig2atilde, crossprod(delta)/2 + sig2b)
    lsig2 <- log(sig2)
    phi2 <- 1/rgamma(1, phi2atilde, crossprod(lz - xmat%*%drop(beta) - smat%*%drop(delta))/2 + bphi)
    lphi2 <- log(phi2)
    mn <- drop(mubd + Sigma12%*%chol2inv(chol(Sigma22))%*%c(lsig2 - mulsig2, lphi2 - mulphi2))
    betadeltaprop <- rmtfixed(1, mn, Sigma, df, Rsigma)
    ljprop <- dmtcholprec(betadeltaprop, mn, Rprec, df, TRUE)
    ljold <-  dmtcholprec(betadelta,     mn, Rprec, df, TRUE)
    lpprop <- lnormlpostiid(c(betadeltaprop, lsig2, lphi2), datlist)
    lpold <-  lnormlpostiid(c(betadelta,     lsig2, lphi2), datlist)
    la <- lpprop - lpold + ljold - ljprop
    u <- runif(1)
    if(tune & lpbest < lpprop){
      lpbest <- lpprop
      mu <- c(betadeltaprop, log(sig2), log(phi2))
      mbd <- mu[1:(nbeta + ndelta)]
      hess <- lnormlposthessiid(mu, datlist)
      SigmaBig <- chol2inv(chol(-hess))
      Sigma1 <- SigmaBig[1:(nbeta + ndelta), 1:(nbeta + ndelta)]
      Sigma12 <- SigmaBig[1:(nbeta + ndelta), nbeta + ndelta + 1:2]
      Sigma22 <- SigmaBig[nbeta + ndelta + 1:2, nbeta + ndelta + 1:2]
      Sigma <- Sigma1 - Sigma12%*%tcrossprod(chol2inv(chol(Sigma22)),Sigma12)
      Rsigma <- chol(Sigma)
      Rprec <- chol(chol2inv(Rsigma))
    }
    if(tune & lpbest < lpold){
      lpbest <- lpold
      mu <- c(betadelta, log(sig2), log(phi2))
      mbd <- mu[1:(nbeta + ndelta)]
      hess <- lnormlposthessiid(mu, datlist)
      SigmaBig <- chol2inv(chol(-hess))
      Sigma1 <- SigmaBig[1:(nbeta + ndelta), 1:(nbeta + ndelta)]
      Sigma12 <- SigmaBig[1:(nbeta + ndelta), nbeta + ndelta + 1:2]
      Sigma22 <- SigmaBig[nbeta + ndelta + 1:2, nbeta + ndelta + 1:2]
      Sigma <- Sigma1 - Sigma12%*%tcrossprod(chol2inv(chol(Sigma22)),Sigma12)
      Rsigma <- chol(Sigma)
      Rprec <- chol(chol2inv(Rsigma))
    }
    if(log(u) < la){
      betadelta <- betadeltaprop
      beta <- betadelta[1:nbeta]
      delta <- betadelta[1:ndelta + nbeta]
      acc[iter] <- 1
    }
    lpbests[iter] <- lpbest
    par <- c(betadelta, sig2, phi2)
    draws[iter + 1,] <- par
  }
  return(list(draws = draws, acc = acc, state = par, lpbest = lpbest, lpbests = lpbests,
              mu = mu, Sigma = Sigma, Rsigma = Rsigma))
}

## MIH w/in gibbs for poisson data model, iid random effects
## uses prior for proposal distribution
poisiidindpriorpropgibbs <- function(niter, init, datlist){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  betasd <- sqrt(betavar)
  sig2a <- datlist$sig2a
  sig2b <- datlist$sig2b
  z <- dat[, 1]
  xmat <- dat[, 1 + 1:nbeta]
  smat <- dat[, 1 + nbeta + 1:ndelta]
  betadelta <- init[1:(nbeta + ndelta)]
  beta <- init[1:nbeta]
  delta <- init[1:ndelta + nbeta]
  npars <- length(init)
  draws <- matrix(0, ncol = npars, nrow = niter + 1)
  colnames(draws) <- c(paste("beta", 1:nbeta, sep = "."),
                       paste("delta", 1:ndelta, sep = "."), "sig2")
  draws[1,] <- init
  accs <- matrix(0, nrow = niter, ncol = nbeta + ndelta)
  sig2atilde <- ndelta/2 + sig2a
  for(iter in 1:niter){
    ##if(iter %% 50 == 0){
    ##  print(iter)
    ##}
    sig2 <- 1/rgamma(1, sig2atilde, crossprod(delta)/2 + sig2b)
    betaprop <- beta
    for(i in 1:nbeta){
      betaprop[i] <- betamn + rnorm(1)*betasd
      ljprop <- dnorm(betaprop[i], betamn, betasd, TRUE)
      ljold <- dnorm(beta[i], betamn, betasd, TRUE)
      lpprop <- poislpostiid(c(betaprop, delta, log(sig2)), datlist)
      lpold <-  poislpostiid(c(beta,     delta, log(sig2)), datlist)
      la <- lpprop - lpold + ljold - ljprop
      u <- runif(1)
      if(log(u) < la){
        beta[i] <- betaprop[i]
        accs[iter,i] <- 1
      } else {
        betaprop[i] <- beta[i]
      }
    }
    deltaprop <- delta
    for(i in 1:ndelta){
      deltaprop[i] <- rnorm(1)*sqrt(sig2)
      ljprop <- dnorm(deltaprop[i], 0, sqrt(sig2), TRUE)
      ljold <- dnorm(delta[i], 0, sqrt(sig2), TRUE)
      lpprop <- poislpostiid(c(beta,     deltaprop, log(sig2)), datlist)
      lpold <-  poislpostiid(c(beta,     delta,     log(sig2)), datlist)
      la <- lpprop - lpold + ljold - ljprop
      u <- runif(1)
      if(log(u) < la){
        delta[i] <- deltaprop[i]
        accs[iter,i + nbeta] <- 1
      } else {
        deltaprop[i] <- delta[i]
      }
    }
    par <- c(beta, delta, sig2)
    draws[iter + 1,] <- par
  }
  return(list(draws = draws, accs = accs, state = par))
}

## MIH w/in gibbs for poisson data model, fully correlated random effects
## conditional normal draw
poisfullindgibbs <- function(niter, init, mu, df, lpbest, datlist, tune = TRUE){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  omdf <- datlist$omdf
  ominvscale <- datlist$ominvscale
  z <- dat[, 1]
  xmat <- dat[, 1 + 1:nbeta]
  smat <- dat[, 1 + nbeta + 1:ndelta]
  betadelta <- init[1:(nbeta + ndelta)]
  beta <- init[1:nbeta]
  delta <- init[1:ndelta + nbeta]
  npars <- length(init)
  nell <- ndelta*(ndelta + 1)/2
  draws <- matrix(0, ncol = npars, nrow = niter + 1)
  ellnums <- apply(which(lower.tri(diag(ndelta), TRUE), TRUE), 1, paste, collapse = ".")
  colnames(draws) <- c(paste("beta", 1:nbeta, sep = "."),
                       paste("delta", 1:ndelta, sep = "."),
                       paste("ell", ellnums, sep = "."))
  lpbests <- rep(0, niter)
  draws[1,] <- init
  acc <- rep(0, niter)
  mubd <- mu[1:(nbeta + ndelta)]
  muell <- mu[nbeta + ndelta + 1:nell]
  L <- matrix(0, ndelta, ndelta)
  L[lower.tri(L, TRUE)] <- muell
  omega <- tcrossprod(L)
  hess <- poislposthessfull(mu, datlist)
  SigmaBig <- chol2inv(chol(-hess))
  Sigma1 <- SigmaBig[1:(nbeta + ndelta), 1:(nbeta + ndelta)]
  Sigma12 <- SigmaBig[1:(nbeta + ndelta), nbeta + ndelta + 1:nell]
  Sigma22 <- SigmaBig[nbeta + ndelta + 1:nell, nbeta + ndelta + 1:nell]
  Sigma12Sigma22 <- Sigma12%*%chol2inv(chol(Sigma22))
  Sigma <- Sigma1 - tcrossprod(Sigma12Sigma22, Sigma12)
  Rsigma <- chol(Sigma)
  Rprec <- chol(chol2inv(Rsigma))  
  for(iter in 1:niter){
    ##if(iter %% 50 == 0){
    ##  print(iter)
    ##}
    L <- rwishchol(omdf + 1, chol2inv(chol(tcrossprod(delta) + ominvscale)))
    ell <- L[lower.tri(L, TRUE)]
    omega <- tcrossprod(L)
    mn <- drop(mubd + Sigma12Sigma22%*%(ell - muell))
    betadeltaprop <- rmtfixed(1, mn, Sigma, df, Rsigma)
    ljprop <- dmtcholprec(betadeltaprop, mn, Rprec, df, TRUE)
    ljold <-  dmtcholprec(betadelta,     mn, Rprec, df, TRUE)
    lpprop <- poislpostfullomega(betadeltaprop, omega, datlist)
    lpold <-  poislpostfullomega(betadelta,     omega, datlist)
    la <- lpprop - lpold + ljold - ljprop
    u <- runif(1)
    if(tune & lpbest < lpprop){
      lpbest <- lpprop
      mu <- c(betadeltaprop, ell)
      hess <- poislposthessfull(mu, datlist)
      SigmaBig <- chol2inv(chol(-hess))
      Sigma1 <- SigmaBig[1:(nbeta + ndelta), 1:(nbeta + ndelta)]
      Sigma12 <- SigmaBig[1:(nbeta + ndelta), nbeta + ndelta + 1:nell]
      Sigma22 <- SigmaBig[nbeta + ndelta + 1:nell, nbeta + ndelta + 1:nell]
      Sigma12Sigma22 <- Sigma12%*%chol2inv(chol(Sigma22))
      Sigma <- Sigma1 - tcrossprod(Sigma12Sigma22, Sigma12)
      Rsigma <- chol(Sigma)
      Rprec <- chol(chol2inv(Rsigma))      
      mubd <- mu[1:(nbeta + ndelta)]
      muell <- mu[nbeta + ndelta + 1:nell]
    }
    if(tune & lpbest < lpold){
      lpbest <- lpold
      mu <- c(betadelta, ell)
      hess <- poislposthessfull(mu, datlist)
      SigmaBig <- chol2inv(chol(-hess))
      Sigma1 <- SigmaBig[1:(nbeta + ndelta), 1:(nbeta + ndelta)]
      Sigma12 <- SigmaBig[1:(nbeta + ndelta), nbeta + ndelta + 1:nell]
      Sigma22 <- SigmaBig[nbeta + ndelta + 1:nell, nbeta + ndelta + 1:nell]
      Sigma12Sigma22 <- Sigma12%*%chol2inv(chol(Sigma22))
      Sigma <- Sigma1 - tcrossprod(Sigma12Sigma22, Sigma12)
      Rsigma <- chol(Sigma)
      Rprec <- chol(chol2inv(Rsigma))      
      mubd <- mu[1:(nbeta + ndelta)]
      muell <- mu[nbeta + ndelta + 1:nell]
    }
    if(log(u) < la){
      betadelta <- betadeltaprop
      beta <- betadelta[1:nbeta]
      delta <- betadelta[1:ndelta + nbeta]
      acc[iter] <- 1
    }
    lpbests[iter] <- lpbest
    par <- c(betadelta, ell)
    draws[iter + 1,] <- par
  }
  return(list(draws = draws, acc = acc, state = par, lpbest = lpbest, lpbests = lpbests,
              mu = mu, Sigma = Sigma, Rsigma = Rsigma))
}

## MIH w/in gibbs for lognormal data model, full random effects
## uses conditional distribution from normal approximation
lnormfullindgibbs <- function(niter, init, mu, df, lpbest, datlist, tune = TRUE){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  omdf <- datlist$omdf
  ominvscale <- datlist$ominvscale
  aphi <- datlist$aphi
  bphi <- datlist$bphi
  z <- dat[, 1]
  lz <- log(z)
  ndat <- length(z)
  xmat <- as.matrix(dat[, 1 + 1:nbeta])
  smat <- as.matrix(dat[, 1 + nbeta + 1:ndelta])
  betadelta <- init[1:(nbeta + ndelta)]
  beta <- drop(init[1:nbeta])
  delta <- init[1:ndelta + nbeta]
  npars <- length(init)
  nell <- ndelta*(ndelta + 1)/2
  draws <- matrix(0, ncol = npars, nrow = niter + 1)
  ellnums <- apply(which(lower.tri(diag(ndelta), TRUE), TRUE), 1, paste, collapse = ".")
  colnames(draws) <- c(paste("beta", 1:nbeta, sep = "."),
                       paste("delta", 1:ndelta, sep = "."),
                       paste("ell", ellnums, sep = "."), "phi2")
  lpbests <- rep(0, niter)
  draws[1,] <- init
  acc <- rep(0, niter)
  phi2atilde <- ndat/2 + aphi
  mubd <- mu[1:(nbeta + ndelta)]
  muell <- mu[nbeta + ndelta + 1:nell]
  mulphi2 <- mu[nbeta + ndelta + nell + 1]
  beta <- mu[1:nbeta]
  delta <- mu[1:ndelta + nbeta]
  L <- matrix(0, ndelta, ndelta)
  L[lower.tri(L, TRUE)] <- muell
  hess <- lnormlposthessfull(mu, datlist)
  SigmaBig <- chol2inv(chol(-hess))
  Sigma1 <- SigmaBig[1:(nbeta + ndelta), 1:(nbeta + ndelta)]
  Sigma12 <- SigmaBig[1:(nbeta + ndelta), nbeta + ndelta + 0:nell + 1]
  Sigma22 <- SigmaBig[nbeta + ndelta + 0:nell + 1, nbeta + ndelta + 0:nell + 1]
  Sigma12Sigma22 <- Sigma12%*%chol2inv(chol(Sigma22))
  Sigma <- Sigma1 - tcrossprod(Sigma12Sigma22, Sigma12)
  Sigma <- (Sigma + t(Sigma))/2
  Rsigma <- chol(Sigma)
  Rprec <- chol(chol2inv(Rsigma))    
  for(iter in 1:niter){
    ##if(iter %% 50 == 0){
    ##  print(iter)
    ##}
    L <- rwishchol(omdf + 1, chol2inv(chol(tcrossprod(delta) + ominvscale)))
    ell <- L[lower.tri(L, TRUE)]
    phi2 <- 1/rgamma(1, phi2atilde, crossprod(lz - xmat%*%drop(beta) - smat%*%drop(delta))/2 +
                                    bphi)
    mn <- drop(mubd + Sigma12Sigma22%*%(c(ell, log(phi2)) - c(muell, mulphi2)))
    betadeltaprop <- rmtfixed(1, mn, Sigma, df, Rsigma)
    ljprop <- dmtcholprec(betadeltaprop, mn, Rprec, df, TRUE)
    ljold <-  dmtcholprec(betadelta,     mn, Rprec, df, TRUE)
    lpprop <- lnormlpostfull(c(betadeltaprop, ell, log(phi2)), datlist)
    lpold <-  lnormlpostfull(c(betadelta,     ell, log(phi2)), datlist)
    la <- lpprop - lpold + ljold - ljprop
    u <- runif(1)
    if(tune & lpbest < lpprop){
      lpbest <- lpprop
      mu <- c(betadeltaprop, ell, log(phi2))
      mbd <- mu[1:(nbeta + ndelta)]
      hess <- lnormlposthessfull(mu, datlist)
      SigmaBig <- chol2inv(chol(-hess))
      Sigma1 <- SigmaBig[1:(nbeta + ndelta), 1:(nbeta + ndelta)]
      Sigma12 <- SigmaBig[1:(nbeta + ndelta), nbeta + ndelta + 0:nell + 1]
      Sigma22 <- SigmaBig[nbeta + ndelta + 0:nell + 1, nbeta + ndelta + 0:nell + 1]
      Sigma12Sigma22 <- Sigma12%*%chol2inv(chol(Sigma22))
      Sigma <- Sigma1 - tcrossprod(Sigma12Sigma22, Sigma12)
      Sigma <- (Sigma + t(Sigma))/2      
      Rsigma <- chol(Sigma)
      Rprec <- chol(chol2inv(Rsigma))        
    }
    if(tune & lpbest < lpold){
      lpbest <- lpold
      mu <- c(betadelta, ell, log(phi2))
      mbd <- mu[1:(nbeta + ndelta)]
      hess <- lnormlposthessfull(mu, datlist)
      SigmaBig <- chol2inv(chol(-hess))
      Sigma1 <- SigmaBig[1:(nbeta + ndelta), 1:(nbeta + ndelta)]
      Sigma12 <- SigmaBig[1:(nbeta + ndelta), nbeta + ndelta + 0:nell + 1]
      Sigma22 <- SigmaBig[nbeta + ndelta + 0:nell + 1, nbeta + ndelta + 0:nell + 1]
      Sigma12Sigma22 <- Sigma12%*%chol2inv(chol(Sigma22))
      Sigma <- Sigma1 - tcrossprod(Sigma12Sigma22, Sigma12)
      Sigma <- (Sigma + t(Sigma))/2
      Rsigma <- chol(Sigma)
      Rprec <- chol(chol2inv(Rsigma))  
    }
    if(log(u) < la){
      betadelta <- betadeltaprop
      beta <- betadelta[1:nbeta]
      delta <- betadelta[1:ndelta + nbeta]
      acc[iter] <- 1
    }
    lpbests[iter] <- lpbest
    par <- c(betadelta, ell, phi2)
    draws[iter + 1,] <- par
  }
  return(list(draws = draws, acc = acc, state = par, lpbest = lpbest, lpbests = lpbests,
              mu = mu, Sigma = Sigma, Rsigma = Rsigma))
}


## MIH w/in gibbs for poisson data model, fully correlated random effects
## proposal from prior
poisfullindpriorpropgibbs <- function(niter, init, datlist){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  betasd <- sqrt(betavar)
  omdf <- datlist$omdf
  ominvscale <- datlist$ominvscale
  z <- dat[, 1]
  xmat <- dat[, 1 + 1:nbeta]
  smat <- dat[, 1 + nbeta + 1:ndelta]
  betadelta <- init[1:(nbeta + ndelta)]
  beta <- init[1:nbeta]
  delta <- init[1:ndelta + nbeta]
  npars <- length(init)
  nell <- ndelta*(ndelta + 1)/2
  draws <- matrix(0, ncol = npars, nrow = niter + 1)
  ellnums <- apply(which(lower.tri(diag(ndelta), TRUE), TRUE), 1, paste, collapse = ".")
  colnames(draws) <- c(paste("beta", 1:nbeta, sep = "."),
                       paste("delta", 1:ndelta, sep = "."),
                       paste("ell", ellnums, sep = "."))
  draws[1,] <- init
  accs <- matrix(0, nrow = niter, ncol = nbeta + ndelta)
  for(iter in 1:niter){
    ##if(iter %% 50 == 0){
    ##  print(iter)
    ##}
    L <- rwishchol(omdf + 1, chol2inv(chol(tcrossprod(delta) + ominvscale)))
    ell <- L[lower.tri(L, TRUE)]
    omega <- tcrossprod(L)
    betaprop <- betamn + rnorm(nbeta)*betasd
    deltaprop <- drop(L%*%rnorm(ndelta))
    Sigma <- chol2inv(t(L))
    betaprop <- beta
    for(i in 1:nbeta){
      betaprop[i] <- betamn + rnorm(1)*betasd
      ljprop <- dnorm(betaprop[i], betamn, betasd, TRUE)
      ljold <- dnorm(beta[i], betamn, betasd, TRUE)
      lpprop <- poislpostfullomega(c(betaprop, delta), omega, datlist)
      lpold <-  poislpostfullomega(c(beta,     delta), omega, datlist)
      la <- lpprop - lpold + ljold - ljprop
      u <- runif(1)
      if(log(u) < la){
        beta[i] <- betaprop[i]
        accs[iter,i] <- 1
      } else {
        betaprop[i] <- beta[i]
      }
    }
    deltaprop <- delta
    for(i in 1:ndelta){
      Sigma11 <- Sigma[i,i]
      Sigma12 <- Sigma[i,-i]
      Sigma22 <- Sigma[-i,-i]
      Sigma12Sigma22inv <- Sigma12%*%chol2inv(chol(Sigma22))
      mnii <- drop(Sigma12Sigma22inv%*%delta[-i])
      varii <- drop(tcrossprod(Sigma12Sigma22inv, Sigma12))
      sdii <- sqrt(varii)
      deltaprop[i] <- mnii + rnorm(1)*sdii 
      ljprop <- dnorm(deltaprop[i], mnii, sdii, TRUE)
      ljold <- dnorm(delta[i], mnii, sdii, TRUE)
      lpprop <- poislpostfullomega(c(beta, deltaprop), omega, datlist)
      lpold <-  poislpostfullomega(c(beta, delta),     omega, datlist)
      la <- lpprop - lpold + ljold - ljprop
      u <- runif(1)
      if(log(u) < la){
        delta[i] <- deltaprop[i]
        accs[iter,i + nbeta] <- 1
      } else {
        deltaprop[i] <- delta[i]
      }
    }
    par <- c(beta, delta, ell)
    draws[iter + 1,] <- par
  }
  return(list(draws = draws, accs = accs, state = par))
}

## generic full independent metropolis
## needs log posterior (lpost) and hessian of logposterior (lposthess)
## assumes all parameters are unconstrained
indmetrop <- function(niter, lpost, lposthess, init, mu, df, lpbest, ..., tune = TRUE, w = 0, Sigma = NULL){
  npar <- length(init)
  par <- init
  if(is.null(Sigma)){
    H <- lposthess(mu, ...)
    Rprec <- chol(-H)
    Sigma <- chol2inv(Rprec)
    Rsig <- chol(Sigma)
  } else {
    Rsig <- chol(Sigma)
    Rprec <- chol(chol2inv(Rsig))        
  }
  draws <- matrix(0, ncol = npar, nrow = niter + 1)
  colnames(draws) <- paste("Par", 1:npar, sep=".")
  acc <- rep(0, niter)
  lpbests <- rep(0, niter)
  draws[1,] <- init
  lpold <- lpost(par, ...)
  ljold <- dmtcholprec(par, mu, Rprec, df, TRUE)
  for(iter in 1:niter){
    ##if(iter %% 1000 == 0){
    ##      print(iter)
    ##}
    prop <- rmtfixed(1, mu, Sigma, df, Rsig)
    lpprop <- lpost(prop, ...)
    ljprop <- dmtcholprec(prop, mu, Rprec, df, TRUE)
    la <- lpprop - lpold + ljold - ljprop
    u <- runif(1)
    if(log(u) < la){
      par <- prop
      acc[iter] <- 1
      lpold <- lpprop
      ljold <- ljprop
    }
    if(tune & lpbest < lpprop){
      lpbest <- lpprop
      mu <- prop
      H <- lposthess(mu, ...)
      Rprec <- chol(-H)
      Sigma <- chol2inv(Rprec)
      Rsig <- chol(Sigma)
    }
    lpbests[iter] <- lpbest
    draws[iter + 1,] <- par
    if(w>0){
      Sigma <- (1-w)*Sigma + w*tcrossprod(par - mu)
      Rsig <- chol(Sigma)
      Rprec <- chol(chol2inv(Rsig))        
    }
  }
  return(list(draws = draws, acc = acc, state = par, lpbests = lpbests, lpbest = lpbest,
              mu = mu, Rsig = Rsig, Sigma = Sigma))
}
