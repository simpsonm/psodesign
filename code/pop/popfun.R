library(mnormt)
library(matrixcalc)

## log posterior for poisson data model, iid random effects, 
## in terms of log(sigma^2)
poislpostiid <- function(pars, datlist){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  sig2a <- datlist$sig2a
  sig2b <- datlist$sig2b
  z <- dat[,1]
  xmat <- dat[, 1 + 1:nbeta]
  beta <- pars[1:nbeta]
  if(ndelta > 0){
    delta <- pars[nbeta + 1:ndelta]
    smat <- dat[, 1 + nbeta + 1:ndelta]
    lsig2 <- pars[nbeta + ndelta + 1]
    y <- as.matrix(xmat)%*%beta + as.matrix(smat)%*%delta
    sig2 <- exp(lsig2)
  } else {
    y <- as.matrix(xmat)%*%beta
  }
  lambda <- exp(y)
  lampart <- sum(z*y - lambda)
  if(ndelta > 0){
    latentpart <- - ndelta/2*lsig2 - crossprod(delta)/(2*sig2)
    priorpart <-  - crossprod(beta - betamn)/(2*betavar) - sig2a*lsig2 - sig2b/sig2
    out <- lampart + latentpart + priorpart
  } else {
    priorpart <-  - crossprod(beta - betamn)/(2*betavar)
    out <- lampart + priorpart
  }
  if(is.null(out) | is.na(out)){
    return(-Inf)
  }
  return(out)
}

## log posterior for poisson data model, fully correlated random effects,
## in terms of L = chol(omega), LL' = omega
poislpostfull <- function(pars, datlist){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  omdf <- datlist$omdf
  ominvscale <- datlist$ominvscale
  z <- dat[,1]
  xmat <- dat[,1 + 1:nbeta]
  beta <- pars[1:nbeta]
  if(ndelta > 0){
    nell <- ndelta*(ndelta + 1)/2
    delta <- pars[nbeta + 1:ndelta]
    smat <- dat[, 1 + nbeta + 1:ndelta]
    ell <- pars[nbeta + ndelta + 1:nell]
    L <- matrix(0, ndelta, ndelta)
    L[lower.tri(L, diag = TRUE)] <- ell
    elldiag <- diag(L)
    omega <- crossprod(L)
    y <- as.matrix(xmat)%*%beta + as.matrix(smat)%*%delta
  } else {
    y <- as.matrix(xmat)%*%beta
  }
  lambda <- exp(y)
  lampart <- sum(z*y - lambda)
  betapart <- -crossprod(beta - betamn)/(2*betavar)
  if(ndelta > 0){
    ellpart <- sum((omdf + 1 - 1:ndelta)/2*log(elldiag^2)) -
      (crossprod(crossprod(L, delta)) + sum(diag(ominvscale%*%omega)))/2
    out <- lampart + ellpart + betapart
  } else {
    out <- lampart + betapart
  }
  if(is.null(out) | is.na(out)){
    return(-Inf)
  }
  return(out)
}

## gradient log posterior for poisson data model, fully correlated random effects,
## in terms of L = chol(omega), LL' = omega
poisgradfull <- function(pars, datlist){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  omdf <- datlist$omdf
  ominvscale <- datlist$ominvscale
  z <- dat[,1]
  xmat <- dat[,1 + 1:nbeta]
  beta <- pars[1:nbeta]
  M <- datlist$M
  K <- datlist$K
  N2 <- datlist$N2
  R <- datlist$R
  MKprime <- datlist$MKprime
  dvelldvell <- datlist$dvelldvell
  if(ndelta > 0){
    nell <- ndelta*(ndelta + 1)/2
    delta <- pars[nbeta + 1:ndelta]
    smat <- dat[, 1 + nbeta + 1:ndelta]
    ell <- pars[nbeta + ndelta + 1:nell]
    L <- matrix(0, ndelta, ndelta)
    L[lower.tri(L, diag = TRUE)] <- ell
    elldiag <- diag(L)
    omega <- crossprod(L)
    y <- as.matrix(xmat)%*%beta + as.matrix(smat)%*%delta
    Indelta <- diag(ndelta)
    dvelldvell <- - (MKprime%*%tcrossprod(kronecker(tcrossprod(delta), Indelta), MKprime)) +
      dvelldvell
  } else {
    y <- as.matrix(xmat)%*%beta
  }
  lambda <- exp(y)
  out <- rep(0, length(pars))
  if(nbeta > 1){
    out[1:nbeta] <- apply(drop(z - lambda)*xmat, 2, sum) - (beta - betamn)/betavar
  } else {
    out[1:nbeta] <- sum(drop(z - lambda)*xmat) - (beta - betamn)/betavar
  }
  if(ndelta > 0){
    out[nbeta + 1:ndelta] <- apply(drop(z - lambda)*smat, 2, sum) -
      tcrossprod(crossprod(delta, L), L)
    out[nbeta + ndelta + 1:nell] <- drop(crossprod(ell, dvelldvell)) - R%*%(1/ell)
  }
  return(out)
}

## log posterior for lognormal data model, fully correlated random effects,
## in terms of L = chol(omega), LL' = omega
lnormlpostfull <- function(pars, datlist){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  omdf <- datlist$omdf
  ominvscale <- datlist$ominvscale
  aphi <- datlist$aphi
  bphi <- datlist$bphi
  z <- dat[,1]
  lz <- log(z)
  xmat <- dat[,1 + 1:nbeta]
  beta <- pars[1:nbeta]
  nell <- ndelta*(ndelta + 1)/2
  lphi2 <- pars[nbeta + ndelta + nell + 1]
  phi2 <- exp(lphi2)
  if(ndelta > 0){
    delta <- pars[nbeta + 1:ndelta]
    smat <- dat[, 1 + nbeta + 1:ndelta]
    ell <- pars[nbeta + ndelta + 1:nell]
    L <- matrix(0, ndelta, ndelta)
    L[lower.tri(L, diag = TRUE)] <- ell
    elldiag <- diag(L)
    omega <- crossprod(L)
    y <- as.matrix(xmat)%*%beta + as.matrix(smat)%*%delta
  } else {
    y <- as.matrix(xmat)%*%beta
  }
  datpart <- sum(dnorm(lz, y, sqrt(phi2), log = TRUE))
  betapart <- -crossprod(beta - betamn)/(2*betavar)
  phipart <- -lphi2*aphi - bphi/phi2
  if(ndelta > 0){
    ellpart <- sum((omdf + 1 - 1:ndelta)/2*log(elldiag^2)) -
      (crossprod(crossprod(L, delta)) + sum(diag(ominvscale%*%omega)))/2
    out <- datpart + betapart + phipart + ellpart
  } else {
    out <- datpart + betapart + phipart
  }
  if(is.null(out) | is.na(out)){
    return(-Inf)
  }
  return(out)
}

## log posterior for lognormal data model, fully correlated random effects,
## in terms of L = chol(omega), LL' = omega
lnormlpostiid <- function(pars, datlist){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  sig2a <- datlist$sig2a
  sig2b <- datlist$sig2b
  aphi <- datlist$aphi
  bphi <- datlist$bphi
  z <- dat[,1]
  lz <- log(z)
  xmat <- dat[,1 + 1:nbeta]
  beta <- pars[1:nbeta]
  if(ndelta > 0){
    delta <- pars[nbeta + 1:ndelta]
    smat <- dat[, 1 + nbeta + 1:ndelta]
    lsig2 <- pars[nbeta + ndelta + 1]
    y <- as.matrix(xmat)%*%beta + as.matrix(smat)%*%delta
    sig2 <- exp(lsig2)
    lphi2 <- pars[nbeta + ndelta + 1 + 1]
  } else {
    y <- as.matrix(xmat)%*%beta
    lphi2 <- pars[nbeta + ndelta + 1]
  }
  phi2 <- exp(lphi2)
  datpart <- sum(dnorm(lz, y, sqrt(phi2), log = TRUE))
  betapart <- -crossprod(beta - betamn)/(2*betavar)
  phipart <- -lphi2*aphi - bphi/phi2
  if(ndelta > 0){
    sig2part <- - ndelta/2*lsig2 - crossprod(delta)/(2*sig2) - sig2a*lsig2 - sig2b/sig2
    out <- datpart + betapart + phipart + sig2part
  } else {
    out <- datpart + betapart + phipart
  }
  if(is.null(out) | is.na(out)){
    return(-Inf)
  }
  return(out)
}

## log posterior for poisson data model, fully correlated random effects,
## in terms of omega
## for use in IMH w/in Gibbs along with looking for better estimates of the mode
poislpostfullomega <- function(pars, omega, datlist){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  omdf <- datlist$omdf
  ominvscale <- datlist$ominvscale
  z <- dat[,1]
  xmat <- dat[,1 + 1:nbeta]
  beta <- pars[1:nbeta]
  if(ndelta > 0){
    nell <- ndelta*(ndelta+1)/2
    smat <- dat[, 1 + nbeta + 1:ndelta]
    delta <- pars[nbeta + 1:ndelta]
    diagL <- diag(chol(omega))
    y <- as.matrix(xmat)%*%beta + as.matrix(smat)%*%delta
  } else {
    y <- as.matrix(xmat)%*%beta
  }
  lambda <- exp(y)
  lampart <- sum(z*y - lambda)
  betapart <- -crossprod(beta - betamn)/(2*betavar)
  if(ndelta > 0){
    ompart <- (omdf - ndelta)*sum(log(diagL)) -
      (crossprod(delta, omega)%*%delta + sum(diag(ominvscale%*%omega)))/2
    out <- lampart + ompart + betapart
  } else {
    out <- lampart + betapart
  }
  if(is.null(out) | is.na(out)){
    return(-Inf)
  }
  return(out)
}

## log posterior for poisson data model, fully correlated random effects,
## in terms of omega
## for use in IMH w/in Gibbs when not looking for better estimates of the mode
## (Cheaper computation than above)
poislpostfullbetadelta <- function(pars, datlist){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  omdf <- datlist$omdf
  ominvscale <- datlist$ominvscale
  z <- dat[,1]
  xmat <- dat[,1 + 1:nbeta]
  beta <- pars[1:nbeta]
  if(ndelta > 0){
    nell <- ndelta*(ndelta+1)/2
    smat <- dat[, 1 + nbeta + 1:ndelta]
    delta <- pars[nbeta + 1:ndelta]
    omegalow <- pars[nbeta + ndelta + 1:nell]
    omega <- matrix(0, ndelta, ndelta)
    omega[lower.tri(omega, diag = TRUE)] <- omegalow
    omega <- t(omega)
    omega[lower.tri(omega, diag = TRUE)] <- omegalow
    diagL <- diag(chol(omega))
    y <- as.matrix(xmat)%*%beta + as.matrix(smat)%*%delta
  } else {
    y <- as.matrix(xmat)%*%beta
  }
  lambda <- exp(y)
  lampart <- sum(z*y - lambda)
  betapart <- -crossprod(beta - betamn)/(2*betavar)
  if(ndelta > 0){
    ompart <- -(crossprod(delta, omega)%*%delta)/2
    out <- lampart + ompart + betapart
  } else {
    out <- lampart + betapart
  }
  if(is.null(out) | is.na(out)){
    return(-Inf)
  }
  return(out)
}

## hessian for poisson data model, iid random effects,
## in terms of log(sig2)
poislposthessiid <- function(pars, datlist){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  sig2b <- datlist$sig2b
  npars <- length(pars)
  ndat <- nrow(dat)
  z <- dat[,1]
  xmat <- dat[, 1 + 1:nbeta]
  beta <- pars[1:nbeta]
  if(ndelta > 0){
    smat <- dat[, 1 + nbeta + 1:ndelta]
    delta <- pars[nbeta + 1:ndelta]
    lsig2 <- pars[nbeta + ndelta + 1]
    sig2 <- exp(lsig2)
    y <- as.matrix(xmat)%*%beta + as.matrix(smat)%*%delta
  } else {
    y <- as.matrix(xmat)%*%beta
  }
  lambda <- exp(y)
  xlam <- xmat*drop(lambda)
  xxlam <- crossprod(xmat, xlam)
  dbetadbeta <- - xxlam - diag(1/betavar, nbeta)
  dbetadbeta <- (dbetadbeta + t(dbetadbeta))/2
  if(ndelta > 0){
    slam <- smat*c(lambda)
    sslam <- crossprod(smat, slam)
    xslam <- crossprod(xmat, slam)
    ddeltaddelta <- - sslam - diag(1/sig2, ndelta)
    ddeltaddelta <- (ddeltaddelta + t(ddeltaddelta))/2
    dbetaddelta <- - xslam
    ddeltadlsig2 <- drop(delta/sig2)
    dlsig2dlsig2 <- drop(-(sig2b + crossprod(delta)/2)/sig2)
    out <- matrix(0, npars, npars)
    out[1:nbeta, 1:nbeta] <- dbetadbeta
    out[nbeta + 1:ndelta, nbeta + 1:ndelta] <- ddeltaddelta
    out[1:nbeta, nbeta + 1:ndelta] <- dbetaddelta
    out[nbeta + 1:ndelta, 1:nbeta] <- t(dbetaddelta)
    out[nbeta + ndelta + 1, nbeta + ndelta + 1] <- dlsig2dlsig2
    out[nbeta + 1:ndelta, nbeta + ndelta + 1] <- ddeltadlsig2
    out[nbeta + ndelta + 1, nbeta + 1:ndelta] <- ddeltadlsig2
    if(!is.positive.definite(-out)){
      out[nbeta + 1:ndelta, nbeta + ndelta + 1] <- 0
      out[nbeta + ndelta + 1, nbeta + 1:ndelta] <- 0
    }
  } else {
    out <- dbetadbeta
  }
  return(out)
}

## hessian for lognormal data model, iid random effects,
## in terms of log(sig2)
lnormlposthessiid <- function(pars, datlist){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  sig2b <- datlist$sig2b
  bphi <- datlist$bphi
  npars <- length(pars)
  ndat <- nrow(dat)
  z <- dat[,1]
  lz <- log(z)
  xmat <- dat[, 1 + 1:nbeta]
  beta <- pars[1:nbeta]
  if(ndelta > 0){
    smat <- dat[, 1 + nbeta + 1:ndelta]
    delta <- pars[nbeta + 1:ndelta]
    lsig2 <- pars[nbeta + ndelta + 1]
    sig2 <- exp(lsig2)
    y <- as.matrix(xmat)%*%beta + as.matrix(smat)%*%delta
    lphi2 <- pars[nbeta + ndelta + 1 + 1]
  } else {
    y <- as.matrix(xmat)%*%beta
    lphi2 <- pars[nbeta + ndelta + 1]
  }
  phi2 <- exp(lphi2)
  xx <- crossprod(xmat)
  dbetadbeta <- - xx/phi2 - diag(1/betavar, nbeta)
  dbetadbeta <- (dbetadbeta + t(dbetadbeta))/2
  dlphi2dlphi2 <- - (bphi + crossprod(lz - y)/2)
  if(nbeta > 1){
    dlphi2dbeta <- - apply(drop(lz - y)*xmat, 2, sum)/phi2
  } else {
    dlphi2dbeta <- - sum(drop(lz - y)*xmat)/phi2
  }
  out <- matrix(0, npars, npars)
  out[1:nbeta, 1:nbeta] <- dbetadbeta
  out[nbeta + ndelta + 1 + 1, nbeta + ndelta + 1 + 1] <- dlphi2dlphi2
  out[nbeta + ndelta + 1 + 1, 1:nbeta] <- dlphi2dbeta
  out[1:nbeta, nbeta + ndelta + 1 + 1] <- dlphi2dbeta
  if(ndelta > 0){
    dlphi2ddelta <- - apply(drop(lz - y)*smat, 2, sum)/phi2
    ss <- crossprod(smat)
    xs <- crossprod(xmat, smat)
    ddeltaddelta <- - ss/phi2 - diag(1/sig2, ndelta)
    ddeltaddelta <- (ddeltaddelta + t(ddeltaddelta))/2
    dbetaddelta <- - xs/phi2
    ddeltadlsig2 <- drop(delta/sig2)
    dlsig2dlsig2 <- drop(-(sig2b + crossprod(delta)/2)/sig2)
    out[nbeta + 1:ndelta, nbeta + 1:ndelta] <- ddeltaddelta
    out[1:nbeta, nbeta + 1:ndelta] <- dbetaddelta
    out[nbeta + 1:ndelta, 1:nbeta] <- t(dbetaddelta)
    out[nbeta + ndelta + 1, nbeta + ndelta + 1] <- dlsig2dlsig2
    out[nbeta + 1:ndelta, nbeta + ndelta + 1] <- ddeltadlsig2
    out[nbeta + ndelta + 1, nbeta + 1:ndelta] <- ddeltadlsig2
    if(!is.positive.definite(-out)){
      out[nbeta + 1:ndelta, nbeta + ndelta + 1] <- 0
      out[nbeta + ndelta + 1, nbeta + 1:ndelta] <- 0
    }
  }
  return(out)
}

## hessian for poisson data model, iid random effects,
## in terms of sig2
## only outputs hessian for beta and delta
poislposthessiidbetadelta <- function(pars, datlist){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  npars <- length(pars)
  ndat <- nrow(dat)
  z <- dat[,1]
  xmat <- dat[, 1 + 1:nbeta]
  beta <- pars[1:nbeta]
  if(ndelta > 0){
    smat <- dat[, 1 + nbeta + 1:ndelta]
    delta <- pars[nbeta + 1:ndelta]
    lsig2 <- pars[nbeta + ndelta + 1]
    sig2 <- exp(lsig2)
    y <- as.matrix(xmat)%*%beta + as.matrix(smat)%*%delta
  } else {
    y <- as.matrix(xmat)%*%beta
  }
  lambda <- exp(y)
  xlam <- xmat*c(lambda)
  xxlam <- crossprod(xmat, xlam)
  dbetadbeta <- - xxlam - diag(1/betavar, nbeta)
  if(ndelta > 0){
    slam <- smat*c(lambda)
    sslam <- crossprod(smat, slam)
    xslam <- crossprod(xmat, slam)
    ddeltaddelta <- - sslam - diag(1/sig2, ndelta)
    dbetaddelta <- - xslam
    out <- matrix(0, nbeta + ndelta, nbeta + ndelta)
    out[1:nbeta, 1:nbeta] <- dbetadbeta
    out[nbeta + 1:ndelta, nbeta + 1:ndelta] <- ddeltaddelta
    out[1:nbeta, nbeta + 1:ndelta] <- dbetaddelta
    out[nbeta + 1:ndelta, 1:nbeta] <- t(dbetaddelta)
  } else {
    out <- dbetadbeta
  }
  return(out)
}

## hessian for poisson data model, fully correlated random effects,
## in terms of L, LL' = omega (cholesky)
poislposthessfull <- function(pars, datlist){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  omdf <- datlist$omdf
  ominvscale <- datlist$ominvscale
  M <- datlist$M
  K <- datlist$K
  N2 <- datlist$N2
  R <- datlist$R
  MKprime <- datlist$MKprime
  dvelldvell <- datlist$dvelldvell
  npars <- length(pars)
  ndat <- nrow(dat)
  z <- dat[, 1]
  xmat <- dat[, 1 + 1:nbeta]
  beta <- pars[1:nbeta]
  if(ndelta > 0){
    nell <- ndelta*(ndelta+1)/2
    smat <- dat[, 1 + nbeta + 1:ndelta]
    delta <- pars[nbeta + 1:ndelta]
    ell <- pars[nbeta + ndelta + 1:nell]
    L <- matrix(0, ndelta, ndelta)
    L[lower.tri(L, diag = TRUE)] <- ell
    elldiag <- diag(L)
    omega <- crossprod(L)
    y <- as.matrix(xmat)%*%beta + as.matrix(smat)%*%delta
  } else {
    y <- as.matrix(xmat)%*%beta
  }
  lambda <- exp(y)
  xlam <- xmat*c(lambda)
  xxlam <- crossprod(xmat, xlam)
  dbetadbeta <- - xxlam - diag(1/betavar, nbeta)
  dbetadbeta <- (dbetadbeta + t(dbetadbeta))/2
  if(ndelta > 0){
    slam <- smat*c(lambda)
    sslam <- crossprod(smat, slam)
    xslam <- crossprod(xmat, slam)
    ddeltaddelta <- - sslam - omega
    ddeltaddelta <- (ddeltaddelta + t(ddeltaddelta))/2
    dbetaddelta <- - xslam
    Indelta <- diag(ndelta)
    ddeltadvell <- -kronecker(t(delta), Indelta)%*%N2%*%tcrossprod(kronecker(L, Indelta), M)
    dvelldvell <- - (MKprime%*%tcrossprod(kronecker(tcrossprod(delta), Indelta), MKprime)) -
      R%*%diag(1/ell^2) + dvelldvell
    dvelldvell <- (dvelldvell + t(dvelldvell))/2
    out <- matrix(0, npars, npars)
    out[1:nbeta, 1:nbeta] <- dbetadbeta
    out[nbeta + 1:ndelta, nbeta + 1:ndelta] <- ddeltaddelta
    out[1:nbeta, nbeta + 1:ndelta] <- dbetaddelta
    out[nbeta + 1:ndelta, 1:nbeta] <- t(dbetaddelta)
    out[nbeta + ndelta + 1:nell, nbeta + ndelta + 1:nell] <- dvelldvell
    out[nbeta + 1:ndelta, nbeta + ndelta + 1:nell] <- ddeltadvell
    out[nbeta + ndelta + 1:nell, nbeta + 1:ndelta] <- t(ddeltadvell)
    if(!is.positive.definite(-out)){
      out[nbeta + 1:ndelta, nbeta + ndelta + 1:nell] <- 0
      out[nbeta + ndelta + 1:nell, nbeta + 1:ndelta] <- 0
    }
  } else {
    out <- dbetadbeta
  }
  return(out)
}

## hessian for poisson data model, fully correlated random effects,
## in terms of omega
## only returns hessian of beta and delta
poislposthessfullbetadelta <- function(pars, omega, datlist){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  omdf <- datlist$omdf
  ominvscale <- datlist$ominvscale
  ndat <- nrow(dat)
  z <- dat[, 1]
  xmat <- dat[, 1 + 1:nbeta]
  beta <- pars[1:nbeta]
  if(ndelta > 0){
    nell <- ndelta*(ndelta+1)/2
    smat <- dat[, 1 + nbeta + 1:ndelta]
    delta <- pars[nbeta + 1:ndelta]
    y <- as.matrix(xmat)%*%beta + as.matrix(smat)%*%delta
  } else {
    y <- as.matrix(xmat)%*%beta
  }
  lambda <- exp(y)
  xlam <- xmat*c(lambda)
  xxlam <- crossprod(xmat, xlam)
  dbetadbeta <- - xxlam - diag(1/betavar, nbeta)
  if(ndelta > 0){
    slam <- smat*c(lambda)
    sslam <- crossprod(smat, slam)
    xslam <- crossprod(xmat, slam)
    ddeltaddelta <- - sslam - omega
    dbetaddelta <- - xslam
    out <- matrix(0, nbeta + ndelta, nbeta + ndelta)
    out[1:nbeta, 1:nbeta] <- dbetadbeta
    out[nbeta + 1:ndelta, nbeta + 1:ndelta] <- ddeltaddelta
    out[1:nbeta, nbeta + 1:ndelta] <- dbetaddelta
    out[nbeta + 1:ndelta, 1:nbeta] <- t(dbetaddelta)
  } else {
    out <- dbetadbeta
  }
  return(out)
}


## hessian for lognormal data model, fully correlated random effects,
## in terms of L, LL' = omega (cholesky)
lnormlposthessfull <- function(pars, datlist){
  dat <- datlist$dat
  nbeta <- datlist$nbeta
  ndelta <- datlist$ndelta
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  omdf <- datlist$omdf
  ominvscale <- datlist$ominvscale
  M <- datlist$M
  K <- datlist$K
  N2 <- datlist$N2
  R <- datlist$R
  MKprime <- datlist$MKprime
  dvelldvell <- datlist$dvelldvell
  aphi <- datlist$aphi
  bphi <- datlist$bphi
  npars <- length(pars)
  ndat <- nrow(dat)
  z <- dat[, 1]
  lz <- log(z)
  xmat <- dat[, 1 + 1:nbeta]
  beta <- pars[1:nbeta]
  nell <- ndelta*(ndelta+1)/2
  lphi2 <- pars[nbeta + ndelta + nell + 1]
  phi2 <- exp(lphi2)
  if(ndelta > 0){
    smat <- dat[, 1 + nbeta + 1:ndelta]
    delta <- pars[nbeta + 1:ndelta]
    ell <- pars[nbeta + ndelta + 1:nell]
    L <- matrix(0, ndelta, ndelta)
    L[lower.tri(L, diag = TRUE)] <- ell
    elldiag <- diag(L)
    omega <- crossprod(L)
    y <- as.matrix(xmat)%*%beta + as.matrix(smat)%*%delta
  } else {
    y <- as.matrix(xmat)%*%beta
  }
  xx <- crossprod(xmat)
  dbetadbeta <- - xx/phi2 - diag(1/betavar, nbeta)
  dbetadbeta <- (dbetadbeta + t(dbetadbeta))/2
  dlphi2dlphi2 <- - (bphi + crossprod(lz - y)/2)
  if(nbeta > 1){
    dlphi2dbeta <- - apply(drop(lz - y)*xmat, 2, sum)/phi2
  } else {
    dlphi2dbeta <- - sum(drop(lz - y)*xmat)/phi2
  }
  out <- matrix(0, npars, npars)
  out[1:nbeta, 1:nbeta] <- dbetadbeta
  out[nbeta + ndelta + nell + 1, nbeta + ndelta + nell + 1] <- dlphi2dlphi2
  out[nbeta + ndelta + nell + 1, 1:nbeta] <- dlphi2dbeta
  out[1:nbeta, nbeta + ndelta + nell + 1] <- dlphi2dbeta
  if(ndelta > 0){
    dlphi2ddelta <- - apply(drop(lz - y)*smat, 2, sum)/phi2
    ss <- crossprod(smat)
    xs <- crossprod(xmat, smat)
    ddeltaddelta <- - ss/phi2 - omega
    ddeltaddelta <- (ddeltaddelta + t(ddeltaddelta))/2
    dbetaddelta <- - xs/phi2
    Indelta <- diag(ndelta)
    ddeltadvell <- -kronecker(t(delta), Indelta)%*%N2%*%tcrossprod(kronecker(L, Indelta), M)
    dvelldvell <- - (MKprime%*%tcrossprod(kronecker(tcrossprod(delta), Indelta), MKprime)) -
      R%*%diag(1/ell^2) + dvelldvell
    dvelldvell <- (dvelldvell + t(dvelldvell))/2
    out[nbeta + 1:ndelta, nbeta + 1:ndelta] <- ddeltaddelta
    out[1:nbeta, nbeta + 1:ndelta] <- dbetaddelta
    out[nbeta + 1:ndelta, 1:nbeta] <- t(dbetaddelta)
    out[nbeta + ndelta + 1:nell, nbeta + ndelta + 1:nell] <- dvelldvell
    out[nbeta + 1:ndelta, nbeta + ndelta + 1:nell] <- ddeltadvell
    out[nbeta + ndelta + 1:nell, nbeta + 1:ndelta] <- t(ddeltadvell)
    out[nbeta + ndelta + nell + 1, 1:ndelta + nbeta] <- dlphi2ddelta
    out[1:ndelta + nbeta, nbeta + ndelta + nell + 1] <- dlphi2ddelta
    if(!is.positive.definite(-out)){
      out[nbeta + 1:ndelta, nbeta + ndelta + 1:nell] <- 0
      out[nbeta + ndelta + 1:nell, nbeta + 1:ndelta] <- 0
    }
  }
  return(out)
}
