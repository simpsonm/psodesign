mylikefit <- function(pars, trendtype, houston, maxit = 500, reltol = .Machine$double.eps){
  u <- houston$u
  v <- houston$v
  Z <- houston$avg
  N <- length(Z)
  S <- houston$se^2
  loc <- cbind(u,v)
  X <- switch(trendtype, cte = matrix(1, ncol = 1, nrow = N),
              lin = cbind(1, u, v), quad = cbind(1, u, v, u^2, u*v, v^2))
  nbeta <- ncol(X)
  optimout <- optim(pars, krigneglike, method = "BFGS", hessian = TRUE,
                    control = list(reltol = reltol, maxit= maxit),
                    Z = Z, X = X, S = S, loc = loc, N = N,
                    nbeta = nbeta, covfun = expcov2)
  unpar <- optimout$par
  par <- c(unpar[1:nbeta], exp(unpar[nbeta + 1:2]), unpar[nbeta + 3])
  like <- -optimout$value  - N*log(2*pi)/2
  hess <- optimout$hessian
  K <- length(pars)
  ntrans <- K - nbeta - 1
  Gprime <- NULL
  covmat <- NULL
  try({
    uncovmat <- chol2inv(chol(hess))
    Gprime <- diag(c(rep(1, nbeta), par[nbeta + 1:ntrans], 1))
    covmat <- Gprime%*%uncovmat%*%Gprime
  })
  aic <- 2*K - 2*like
  bic <- K*log(N) - 2*like
  optimred <- optim(pars[-K], nospatneglike, method = "BFGS",
                    control = list(reltol = reltol, maxit= maxit),
                    Z = Z, X = X, S = S, N = N, nbeta = nbeta)
  aicred <- 2*(K - 1) + 2*optimred$value + N*log(2*pi)
  bicred <- (K - 1)*log(N) + 2*optimred$value + N*log(2*pi)
  ics <- cbind(c(aic, bic), c(aicred, bicred))
  rownames(ics) <- c("AIC", "BIC")
  colnames(ics) <- c("Full", "Reduced")
  out <- list(par = par, unpar = unpar, like = like, covmat = covmat, ics = ics,
              optimout = optimout, optimred = optimred)
  return(out)
}

## mean full uncertainty universal kriging variance at a set of target points
## (i.e. assuming only latent process and beta is unknown)
sig2fuk.mean <- function(dd, datlist){
  ## dd = variable design points
  ss <- datlist$ss ## fixed points in the design
  tt <- datlist$tt ## target points
  theta <- datlist$theta
  sig2z <- datlist$sig2z
  poly <- datlist$poly
  invCz.s <- datlist$invCz.s
  Cyy.s.t <- datlist$Cyy.s.t
  Cy.t <- datlist$Cy.t
  D.s.t <- datlist$D.s.t
  D.t <- datlist$D.t
  D.s <- datlist$D.s
  N.s <- nrow(ss)
  N.t <- nrow(tt)
  dd <- matrix(dd, ncol = 2)
  N.d <- nrow(dd)
  X.t <- cbind(1, tt)
  X.s <- cbind(1, ss)
  X.d <- cbind(1, dd)
  X <- rbind(X.s, X.d)
  ## check that all design points are in the target county
  check <- point.in.polygon(dd[,1], dd[,2], poly@coords[,1], poly@coords[,2])
  if(sum(check) == N.d){
    ## first, finish creating Cz and Czinv
    D.d <- matrix(0, N.d, N.d)
    for(i in 2:N.d){
      for(j in 1:(i-1)){
        D.d[i,j] <- sqrt(sum((dd[i,] - dd[j,])^2))
        D.d[j,i] <- D.d[i,j]
      }
    }
    D.d.s <- matrix(0, N.d, N.s)
    for(i in 1:N.d){
      for(j in 1:N.s){
        D.d.s[i,j] <- sqrt(sum((dd[i,] - ss[j,])^2))
      }
    }
    D.d.t <- matrix(0, N.d, N.t)
    for(i in 1:N.d){
      for(j in 1:N.t){
        D.d.t[i,j] <- sqrt(sum((dd[i,] - tt[j,])^2))
      }
    }
    Cz.d.s <- theta[1]*exp(-D.d.s/theta[2])
    Cz.d <- diag(sig2z, N.d) + theta[1]*exp(-D.d/theta[2])
    DD <- chol2inv(chol(Cz.d - Cz.d.s%*%tcrossprod(invCz.s, Cz.d.s)))
    AinvB <- tcrossprod(invCz.s, Cz.d.s)
    AA <- invCz.s + AinvB%*%tcrossprod(DD,AinvB)
    BB <- -AinvB%*%DD
    invCz <- rbind(cbind(AA, BB), cbind(t(BB), DD))
    ## next, finish creating Cyy
    Cyy.d.t <- theta[1]*exp(-D.d.t/theta[2])
    Cyy <- rbind(Cyy.s.t, Cyy.d.t)
    D.sd.t <- rbind(D.s.t, D.d.t)
    D.sd <- rbind(cbind(D.s, t(D.d.s)), cbind(D.d.s, D.d))
    outs1 <- apply(Cyy, 2, function(x, invCz) crossprod(x, invCz)%*%x, invCz = invCz)
    XtinvCz <- crossprod(X, invCz)
    prec <- XtinvCz%*%X
    Rprec <- chol(prec)
    invRprec <- backsolve(Rprec, diag(ncol(X)))
    delta <- t(X.t) - XtinvCz%*%Cyy
    outs2 <- apply(delta, 2, function(x, invRprec) tcrossprod(crossprod(x, invRprec)),
                   invRprec = invRprec)
    U <- invCz%*%X%*%tcrossprod(invRprec)
    V <- invCz - U%*%crossprod(X,invCz)
    uv1s <- -V%*%tcrossprod(U, X.t) - V%*%V%*%Cyy + V%*%(D.sd.t == 0)
    dtheta1 <- exp(-D.sd/theta[2])
    d2 <- invCz%*%dtheta1
    dtheta2 <- theta[1]*exp(-D.sd/theta[2])*D.sd/theta[2]^2
    d3 <- invCz%*%dtheta2
    d.sd.t.theta1 <- exp(-D.sd.t/theta[2])
    uv2s <- -V%*%dtheta1%*%tcrossprod(U, X.t) - V%*%dtheta1%*%V%*%Cyy +
      V%*%d.sd.t.theta1
    uv3s <- -V%*%dtheta2%*%tcrossprod(U, X.t) - V%*%dtheta2%*%V%*%Cyy +
      V%*%(theta[1]*d.sd.t.theta1*D.sd.t/theta[2]^2)
    UVs <- array(c(uv1s, uv2s, uv3s), c(nrow(uv1s), ncol(uv1s), 3))
    FI <- matrix(0, 3, 3)
    invCz2 <- invCz%*%invCz
    FI[1,1] <- sum(diag(invCz2))
    FI[1,2] <- sum(diag(invCz%*%d2))
    FI[2,1] <- FI[1,2]
    FI[1,3] <- sum(diag(invCz%*%d3))
    FI[3,1] <- FI[1,3]
    FI[2,2] <- sum(diag(d2%*%d2))
    FI[3,2] <- sum(diag(d2%*%d3))
    FI[2,3] <- FI[3,2]
    FI[3,3] <- sum(diag(d3%*%d3))
    FIinv <- chol2inv(chol(FI))
    parparts <- rep(0, N.t)
    for(i in 1:N.t){
       parparts[i] <- sum(diag(crossprod(UVs[,i,], invCz)%*%UVs[,i,]%*%FIinv))
    }
    out <- Cy.t + mean(outs2) - mean(outs1) + mean(parparts)
  } else {
    out <- Inf
  }
  return(out)
}

## max full uncertainty universal kriging variance at a set of target points
## (i.e. assuming only latent process and beta is unknown)
sig2fuk.max <- function(dd, datlist){
  ## dd = variable design points
  ss <- datlist$ss ## fixed points in the design
  tt <- datlist$tt ## target points
  theta <- datlist$theta
  sig2z <- datlist$sig2z
  poly <- datlist$poly
  invCz.s <- datlist$invCz.s
  Cyy.s.t <- datlist$Cyy.s.t
  Cy.t <- datlist$Cy.t
  D.s.t <- datlist$D.s.t
  D.t <- datlist$D.t
  D.s <- datlist$D.s
  N.s <- nrow(ss)
  N.t <- nrow(tt)
  dd <- matrix(dd, ncol = 2)
  N.d <- nrow(dd)
  X.t <- cbind(1, tt)
  X.s <- cbind(1, ss)
  X.d <- cbind(1, dd)
  X <- rbind(X.s, X.d)
  ## check that all design points are in the target county
  check <- point.in.polygon(dd[,1], dd[,2], poly@coords[,1], poly@coords[,2])
  if(sum(check) == N.d){
    ## first, finish creating Cz and Czinv
    D.d <- matrix(0, N.d, N.d)
    for(i in 2:N.d){
      for(j in 1:(i-1)){
        D.d[i,j] <- sqrt(sum((dd[i,] - dd[j,])^2))
        D.d[j,i] <- D.d[i,j]
      }
    }
    D.d.s <- matrix(0, N.d, N.s)
    for(i in 1:N.d){
      for(j in 1:N.s){
        D.d.s[i,j] <- sqrt(sum((dd[i,] - ss[j,])^2))
      }
    }
    D.d.t <- matrix(0, N.d, N.t)
    for(i in 1:N.d){
      for(j in 1:N.t){
        D.d.t[i,j] <- sqrt(sum((dd[i,] - tt[j,])^2))
      }
    }
    Cz.d.s <- theta[1]*exp(-D.d.s/theta[2])
    Cz.d <- diag(sig2z, N.d) + theta[1]*exp(-D.d/theta[2])
    DD <- chol2inv(chol(Cz.d - Cz.d.s%*%tcrossprod(invCz.s, Cz.d.s)))
    AinvB <- tcrossprod(invCz.s, Cz.d.s)
    AA <- invCz.s + AinvB%*%tcrossprod(DD,AinvB)
    BB <- -AinvB%*%DD
    invCz <- rbind(cbind(AA, BB), cbind(t(BB), DD))
    ## next, finish creating Cyy
    Cyy.d.t <- theta[1]*exp(-D.d.t/theta[2])
    Cyy <- rbind(Cyy.s.t, Cyy.d.t)
    D.sd.t <- rbind(D.s.t, D.d.t)
    D.sd <- rbind(cbind(D.s, t(D.d.s)), cbind(D.d.s, D.d))
    outs1 <- apply(Cyy, 2, function(x, invCz) crossprod(x, invCz)%*%x, invCz = invCz)
    XtinvCz <- crossprod(X, invCz)
    prec <- XtinvCz%*%X
    Rprec <- chol(prec)
    invRprec <- backsolve(Rprec, diag(ncol(X)))
    delta <- t(X.t) - XtinvCz%*%Cyy
    outs2 <- apply(delta, 2, function(x, invRprec) tcrossprod(crossprod(x, invRprec)),
                   invRprec = invRprec)
    U <- invCz%*%X%*%tcrossprod(invRprec)
    V <- invCz - U%*%crossprod(X,invCz)
    uv1s <- -V%*%tcrossprod(U, X.t) - V%*%V%*%Cyy + V%*%(D.sd.t == 0)
    dtheta1 <- exp(-D.sd/theta[2])
    d2 <- invCz%*%dtheta1
    dtheta2 <- theta[1]*exp(-D.sd/theta[2])*D.sd/theta[2]^2
    d3 <- invCz%*%dtheta2
    d.sd.t.theta1 <- exp(-D.sd.t/theta[2])
    uv2s <- -V%*%dtheta1%*%tcrossprod(U, X.t) - V%*%dtheta1%*%V%*%Cyy +
      V%*%d.sd.t.theta1
    uv3s <- -V%*%dtheta2%*%tcrossprod(U, X.t) - V%*%dtheta2%*%V%*%Cyy +
      V%*%(theta[1]*d.sd.t.theta1*D.sd.t/theta[2]^2)
    UVs <- array(c(uv1s, uv2s, uv3s), c(nrow(uv1s), ncol(uv1s), 3))
    FI <- matrix(0, 3, 3)
    invCz2 <- invCz%*%invCz
    FI[1,1] <- sum(diag(invCz2))
    FI[1,2] <- sum(diag(invCz%*%d2))
    FI[2,1] <- FI[1,2]
    FI[1,3] <- sum(diag(invCz%*%d3))
    FI[3,1] <- FI[1,3]
    FI[2,2] <- sum(diag(d2%*%d2))
    FI[3,2] <- sum(diag(d2%*%d3))
    FI[2,3] <- FI[3,2]
    FI[3,3] <- sum(diag(d3%*%d3))
    FIinv <- chol2inv(chol(FI))
    parparts <- rep(0, N.t)
    for(i in 1:N.t){
       parparts[i] <- sum(diag(crossprod(UVs[,i,], invCz)%*%UVs[,i,]%*%FIinv))
    }
    out <- Cy.t + max(outs2 - outs1 + parparts)
  } else {
    out <- Inf
  }
  return(out)
}

## mean negative universal kriging variance at a set of target points
## (i.e. assuming only latent process and beta is unknown)
sig2uk.mean <- function(dd, datlist){
  ## dd = variable design points
  ss <- datlist$ss ## fixed points in the design
  tt <- datlist$tt ## target points
  covfun <- datlist$covfun
  theta <- datlist$theta
  sig2z <- datlist$sig2z
  poly <- datlist$poly
  invCz.s <- datlist$invCz.s
  Cyy.s.t <- datlist$Cyy.s.t
  Cy.t <- datlist$Cy.t
  N.s <- nrow(ss)
  N.t <- nrow(tt)
  dd <- matrix(dd, ncol = 2)
  N.d <- nrow(dd)
  X.t <- cbind(1, tt)
  X.s <- cbind(1, ss)
  X.d <- cbind(1, dd)
  X <- rbind(X.s, X.d)
  ## check that all design points are in the target county
  check <- point.in.polygon(dd[,1], dd[,2], poly@coords[,1], poly@coords[,2])
  if(sum(check) == N.d){
    ## first, finish creating Cz and Czinv
    Cz.d.s <- Cyyfun(dd, ss, N.d, N.s, theta, covfun)
    Cz.d <- Czfun(dd, N.d, theta, sig2z, covfun)
    DD <- chol2inv(chol(Cz.d - Cz.d.s%*%tcrossprod(invCz.s, Cz.d.s)))
    AinvB <- tcrossprod(invCz.s, Cz.d.s)
    AA <- invCz.s + AinvB%*%tcrossprod(DD,AinvB)
    BB <- -AinvB%*%DD
    invCz <- rbind(cbind(AA, BB), cbind(t(BB), DD))
    ## next, finish creating Cyy
    Cyy.d.t <- Cyyfun(dd, tt, N.d, N.t, theta, covfun)
    Cyy <- rbind(Cyy.s.t, Cyy.d.t)
    outs1 <- apply(Cyy, 2, function(x, invCz) crossprod(x, invCz)%*%x, invCz = invCz)
    XtinvCz <- crossprod(X, invCz)
    prec <- XtinvCz%*%X
    Rprec <- chol(prec)
    invRprec <- backsolve(Rprec, diag(ncol(X)))
    delta <- t(X.t) - XtinvCz%*%Cyy
    outs2 <- apply(delta, 2, function(x, invRprec) tcrossprod(crossprod(x, invRprec)),
                   invRprec = invRprec)
    out <- Cy.t + mean(outs2) - mean(outs1)
  } else {
    out <- Inf
  }
  return(out)
}

## max universal kriging variance at a set of target points
## (i.e. assuming only latent process and beta is unknown)
sig2uk.max <- function(dd, datlist){
  ## dd = variable design points
  ss <- datlist$ss ## fixed points in the design
  tt <- datlist$tt ## target points
  covfun <- datlist$covfun
  theta <- datlist$theta
  sig2z <- datlist$sig2z
  poly <- datlist$poly
  invCz.s <- datlist$invCz.s
  Cyy.s.t <- datlist$Cyy.s.t
  Cy.t <- datlist$Cy.t
  N.s <- nrow(ss)
  N.t <- nrow(tt)
  dd <- matrix(dd, ncol = 2)
  N.d <- nrow(dd)
  X.t <- cbind(1, tt)
  X.s <- cbind(1, ss)
  X.d <- cbind(1, dd)
  X <- rbind(X.s, X.d)
  ## check that all design points are in the target county
  check <- point.in.polygon(dd[,1], dd[,2], poly@coords[,1], poly@coords[,2])
  if(sum(check) == N.d){
    ## first, finish creating Cz and Czinv
    Cz.d.s <- Cyyfun(dd, ss, N.d, N.s, theta, covfun)
    Cz.d <- Czfun(dd, N.d, theta, sig2z, covfun)
    DD <- chol2inv(chol(Cz.d - Cz.d.s%*%tcrossprod(invCz.s, Cz.d.s)))
    AinvB <- tcrossprod(invCz.s, Cz.d.s)
    AA <- invCz.s + AinvB%*%tcrossprod(DD,AinvB)
    BB <- -AinvB%*%DD
    invCz <- rbind(cbind(AA, BB), cbind(t(BB), DD))
    ## next, finish creating Cyy
    Cyy.d.t <- Cyyfun(dd, tt, N.d, N.t, theta, covfun)
    Cyy <- rbind(Cyy.s.t, Cyy.d.t)
    outs1 <- apply(Cyy, 2, function(x, invCz) crossprod(x, invCz)%*%x, invCz = invCz)
    XtinvCz <- crossprod(X, invCz)
    prec <- XtinvCz%*%X
    Rprec <- chol(prec)
    invRprec <- backsolve(Rprec, diag(ncol(X)))
    delta <- t(X.t) - XtinvCz%*%Cyy
    outs2 <- apply(delta, 2, function(x, invRprec) tcrossprod(crossprod(x, invRprec)),
                   invRprec = invRprec)
    outs <- outs2 - outs1
    out <- Cy.t + max(outs)
  } else {
    out <- Inf
  }
  return(out)
}


## mean negative universal kriging variance at a set of target points
## (i.e. assuming only latent process and beta is unknown)
negsig2uk.mean <- function(dd, datlist){
  ## dd = variable design points
  ss <- datlist$ss ## fixed points in the design
  tt <- datlist$tt ## target points
  covfun <- datlist$covfun
  theta <- datlist$theta
  sig2z <- datlist$sig2z
  poly <- datlist$poly
  invCz.s <- datlist$invCz.s
  Cyy.s.t <- datlist$Cyy.s.t
  Cy.t <- datlist$Cy.t
  N.s <- nrow(ss)
  N.t <- nrow(tt)
  dd <- matrix(dd, ncol = 2)
  N.d <- nrow(dd)
  X.t <- cbind(1, tt)
  X.s <- cbind(1, ss)
  X.d <- cbind(1, dd)
  X <- rbind(X.s, X.d)
  ## check that all design points are in the target county
  check <- point.in.polygon(dd[,1], dd[,2], poly@coords[,1], poly@coords[,2])
  if(sum(check) == N.d){
    ## first, finish creating Cz and Czinv
    Cz.d.s <- Cyyfun(dd, ss, N.d, N.s, theta, covfun)
    Cz.d <- Czfun(dd, N.d, theta, sig2z, covfun)
    DD <- chol2inv(chol(Cz.d - Cz.d.s%*%tcrossprod(invCz.s, Cz.d.s)))
    AinvB <- tcrossprod(invCz.s, Cz.d.s)
    AA <- invCz.s + AinvB%*%tcrossprod(DD,AinvB)
    BB <- -AinvB%*%DD
    invCz <- rbind(cbind(AA, BB), cbind(t(BB), DD))
    ## next, finish creating Cyy
    Cyy.d.t <- Cyyfun(dd, tt, N.d, N.t, theta, covfun)
    Cyy <- rbind(Cyy.s.t, Cyy.d.t)
    outs1 <- apply(Cyy, 2, function(x, invCz) crossprod(x, invCz)%*%x, invCz = invCz)
    XtinvCz <- crossprod(X, invCz)
    prec <- XtinvCz%*%X
    Rprec <- chol(prec)
    invRprec <- backsolve(Rprec, diag(ncol(X)))
    delta <- t(X.t) - XtinvCz%*%Cyy
    outs2 <- apply(delta, 2, function(x, invRprec) tcrossprod(crossprod(x, invRprec)),
                   invRprec = invRprec)
    out <- mean(outs1) - Cy.t - mean(outs2)
  } else {
    out <- -Inf
  }
  return(out)
}

## min negative universal kriging variance at a set of target points
## (i.e. assuming only latent process and beta is unknown)
negsig2uk.min <- function(dd, datlist){
  ## dd = variable design points
  ss <- datlist$ss ## fixed points in the design
  tt <- datlist$tt ## target points
  covfun <- datlist$covfun
  theta <- datlist$theta
  sig2z <- datlist$sig2z
  poly <- datlist$poly
  invCz.s <- datlist$invCz.s
  Cyy.s.t <- datlist$Cyy.s.t
  Cy.t <- datlist$Cy.t
  N.s <- nrow(ss)
  N.t <- nrow(tt)
  dd <- matrix(dd, ncol = 2)
  N.d <- nrow(dd)
  X.t <- cbind(1, tt)
  X.s <- cbind(1, ss)
  X.d <- cbind(1, dd)
  X <- rbind(X.s, X.d)
  ## check that all design points are in the target county
  check <- point.in.polygon(dd[,1], dd[,2], poly@coords[,1], poly@coords[,2])
  if(sum(check) == N.d){
    ## first, finish creating Cz and Czinv
    Cz.d.s <- Cyyfun(dd, ss, N.d, N.s, theta, covfun)
    Cz.d <- Czfun(dd, N.d, theta, sig2z, covfun)
    DD <- chol2inv(chol(Cz.d - Cz.d.s%*%tcrossprod(invCz.s, Cz.d.s)))
    AinvB <- tcrossprod(invCz.s, Cz.d.s)
    AA <- invCz.s + AinvB%*%tcrossprod(DD,AinvB)
    BB <- -AinvB%*%DD
    invCz <- rbind(cbind(AA, BB), cbind(t(BB), DD))
    ## next, finish creating Cyy
    Cyy.d.t <- Cyyfun(dd, tt, N.d, N.t, theta, covfun)
    Cyy <- rbind(Cyy.s.t, Cyy.d.t)
    outs1 <- apply(Cyy, 2, function(x, invCz) crossprod(x, invCz)%*%x, invCz = invCz)
    XtinvCz <- crossprod(X, invCz)
    prec <- XtinvCz%*%X
    Rprec <- chol(prec)
    invRprec <- backsolve(Rprec, diag(ncol(X)))
    delta <- t(X.t) - XtinvCz%*%Cyy
    outs2 <- apply(delta, 2, function(x, invRprec) tcrossprod(crossprod(x, invRprec)),
                   invRprec = invRprec)
    outs <- outs1 - outs2
    out <- min(outs) - Cy.t
  } else {
    out <- -Inf
  }
  return(out)
}

## mean negative simple kriging variance at a set of target points
## (i.e. assuming only latent process is unknown)
negsig2sk.mean <- function(dd, datlist){
  ## dd = variable design points
  ss <- datlist$ss ## fixed points in the design
  tt <- datlist$tt ## target points
  covfun <- datlist$covfun
  theta <- datlist$theta
  sig2z <- datlist$sig2z
  poly <- datlist$poly
  invCz.s <- datlist$invCz.s
  Cyy.s.t <- datlist$Cyy.s.t
  Cy.t <- datlist$Cy.t
  N.s <- nrow(ss)
  N.t <- nrow(tt)
  dd <- matrix(dd, ncol = 2)
  N.d <- nrow(dd)
  ## check that all design points are in the target county
  check <- point.in.polygon(dd[,1], dd[,2], poly@coords[,1], poly@coords[,2])
  if(sum(check) == N.d){
    ## first, finish creating Cz and Czinv
    Cz.d.s <- Cyyfun(dd, ss, N.d, N.s, theta, covfun)
    Cz.d <- Czfun(dd, N.d, theta, sig2z, covfun)
    DD <- chol2inv(chol(Cz.d - Cz.d.s%*%tcrossprod(invCz.s, Cz.d.s)))
    AinvB <- tcrossprod(invCz.s, Cz.d.s)
    AA <- invCz.s + AinvB%*%tcrossprod(DD,AinvB)
    BB <- -AinvB%*%DD
    invCz <- rbind(cbind(AA, BB), cbind(t(BB), DD))
    ## next, finish creating Cyy
    Cyy.d.t <- Cyyfun(dd, tt, N.d, N.t, theta, covfun)
    Cyy <- rbind(Cyy.s.t, Cyy.d.t)
    outs <- apply(Cyy, 2, function(x, invCz) crossprod(x, invCz)%*%x, invCz = invCz)
    out <- mean(outs) - Cy.t
  } else {
    out <- -Inf
  }
  return(out)
}

## min negative simple kriging variance at a set of target points
## (i.e. assuming only latent process is unknown)
negsig2sk.min <- function(dd, datlist){
  ## dd = variable design points
  ss <- datlist$ss ## fixed points in the design
  tt <- datlist$tt ## target points
  covfun <- datlist$covfun
  theta <- datlist$theta
  sig2z <- datlist$sig2z
  poly <- datlist$poly
  invCz.s <- datlist$invCz.s
  Cyy.s.t <- datlist$Cyy.s.t
  Cy.t <- datlist$Cy.t
  N.s <- nrow(ss)
  N.t <- nrow(tt)
  dd <- matrix(dd, ncol = 2)
  N.d <- nrow(dd)
  ## check that all design points are in the target county
  check <- point.in.polygon(dd[,1], dd[,2], poly@coords[,1], poly@coords[,2])
  if(sum(check) == N.d){
    ## first, finish creating Cz and Czinv
    Cz.d.s <- Cyyfun(dd, ss, N.d, N.s, theta, covfun)
    Cz.d <- Czfun(dd, N.d, theta, sig2z, covfun)
    DD <- chol2inv(chol(Cz.d - Cz.d.s%*%tcrossprod(invCz.s, Cz.d.s)))
    AinvB <- tcrossprod(invCz.s, Cz.d.s)
    AA <- invCz.s + AinvB%*%tcrossprod(DD,AinvB)
    BB <- -AinvB%*%DD
    invCz <- rbind(cbind(AA, BB), cbind(t(BB), DD))
    ## next, finish creating Cyy
    Cyy.d.t <- Cyyfun(dd, tt, N.d, N.t, theta, covfun)
    Cyy <- rbind(Cyy.s.t, Cyy.d.t)
    outs <- apply(Cyy, 2, function(x, invCz) crossprod(x, invCz)%*%x, invCz = invCz)
    out <- min(outs) - Cy.t
  } else {
    out <- -Inf
  }
  return(out)
}

logit <- function(x) return(log(x/(1-x)))

krigneglike <- function(pars, Z, X, S, loc, N, nbeta, covfun, log = TRUE){
  beta <- pars[1:nbeta]
  sig2z <- exp(pars[nbeta + 1])
  theta <- c(exp(pars[nbeta + 2]), pars[nbeta + 3])
  mu <- X%*%beta
  delta <- Z - mu
  Cz <- Czfun(loc, N, theta, sig2z, covfun) + diag(S)
  Rinv <- NULL
  try({
    R <- chol(Cz)
    Rinv <- backsolve(R, diag(N))
  })
  if(is.null(Rinv)){
    out <- Inf
  } else {
    out <- drop(tcrossprod(crossprod(delta, Rinv))/2 + sum(log(diag(R))))
  }
  return(out)
}

nospatneglike <- function(pars, Z, X, S, N, nbeta, log = TRUE){
  beta <- pars[1:nbeta]
  sig2z <- exp(pars[nbeta + 1])
  sig2y <- exp(pars[nbeta + 2])
  mu <- X%*%beta
  delta <- Z - mu
  Cz <- diag(S + sig2y)
  Rinv <- NULL
  try({
    R <- chol(Cz)
    Rinv <- backsolve(R, diag(N))
  })
  if(is.null(Rinv)){
    out <- Inf
  } else {
    out <- drop(tcrossprod(crossprod(delta, Rinv))/2 + sum(log(diag(R))))
  }
  return(out)
}


## powered exponential covariance function
powexpcov <- function(s1, s2, theta){
  theta1 <- theta[1]
  theta2 <- theta[2]
  sig0 <- theta[3]
  sig1 <- theta[4]
  h <- s1 - s2
  euc <- sqrt(sum(h^2))
  out <- sig0*(euc == 0) + sig1*exp(-(euc/theta1)^theta2)
  return(out)
}

## exponential covariance function
expcov <- function(s1, s2, theta){
  theta1 <- theta[1]
  sig0 <- theta[2]
  sig1 <- theta[3]
  h <- s1 - s2
  euc <- sqrt(sum(h^2))
  out <- sig0*(euc == 0) + exp(log(sig1) - (euc/theta1))
  return(out)
}

## exponential covariance function
## note reversed order of phi and sig2
expcov2 <- function(s1, s2, theta){
  phi <- theta[2]
  sig2 <- theta[1]
  h <- s1 - s2
  euc <- sqrt(sum(h^2))
  out <- exp(log(sig2) - (euc/phi))
  return(out)
}


## matern covariance function
materncov <- function(s1, s2, theta){
  theta1 <- theta[1]
  theta2 <- theta[2]
  sig0 <- theta[3]
  sig1 <- theta[4]
  h <- s1 - s2
  euc <- sqrt(sum(h^2))
  hthet <- euc/theta1
  if(euc > 0){
    logbesselpart <- log(besselK(hthet, theta2, TRUE)) - hthet - lgamma(theta2) +
      theta2*log(hthet) - (theta2 - 1)*log(2)
    out <- sig1*exp(logbesselpart)
  } else {
    out <- sig0 + sig1
  }
  return(out)
}

Czfun <- function(d, N, theta, sig2z, covfun){
  Cz <- matrix(0, N, N)
  for(i in 1:N){
    for(j in 1:i){
      Cz[i,j] <- covfun(d[i,], d[j,], theta)
      if(j < i){
        Cz[j,i] <- Cz[i,j]
      } else {
        Cz[i,j] <- Cz[i,j] + sig2z
      }
    }
  }
  return(Cz)
}

Cyfun <- function(s, M, theta, covfun){
  Cy <- matrix(0, M, M)
  for(i in 1:M){
    for(j in 1:i){
      Cy[i,j] <- covfun(s[i,], s[j,], theta)
      if(j < i){
        Cy[j,i] <- Cy[i,j]
      }
    }
  }
  return(Cy)
}

Cyyfun <- function(d, s, N, M, theta, covfun){
  Cyy <- matrix(0, N, M)
  for(i in 1:N){
    for(j in 1:M){
      Cyy[i,j] <- covfun(d[i,], s[j,], theta)
    }
  }
  return(Cyy)
}


## Sigma_{buk}; same as Sigma_{uk} is Sinv = 0 matrix
sigbuk <- function(d, datlist){
  N <- datlist$N
  s <- datlist$PredS
  M <- datlist$M
  Sinv <- datlist$Sinv
  theta <- datlist$theta
  sig2z <- datlist$sig2z
  covfun <- datlist$covfun
  d <- matrix(1/(1 + exp(-d)), ncol = 2)
  Xm <- cbind(1, PredS)
  Xn <- cbind(1, d)
  Cz <- Czfun(d, N, theta, sig2z, covfun)
  Cyy <- Cyyfun(d, s, N, M, theta, covfun)
  Cy <- Cyfun(s, M, theta, covfun)
  Rz <- chol(Cz)
  Rzinv <- backsolve(Rz, diag(N))
  tRzinvXn <- crossprod(Rzinv, Xn)
  tRzinvCyy <- crossprod(Rzinv, Cyy)
  Pgls <- crossprod(tRzinvXn) + Sinv
  delta <- Xm - crossprod(tRzinvCyy, tRzinvXn)
  RP <- chol(Pgls)
  RPinv <- backsolve(RP, diag(3))
  Sbuk <- Cy - crossprod(tRzinvCyy) + tcrossprod(delta%*%RPinv)
  return(Sbuk)
}

logdetbuk <- function(d, datlist){
  Sbuk <- NULL
  try(Sbuk <- sigbuk(d, datlist))
  if(is.null(Sbuk)){
    out <- -Inf
  } else{
    out <- -sum(log(diag(chol(Sbuk))))
  }
  return(out)
}

meantrbuk <- function(d, datlist){
  Sbuk <- NULL
  try(Sbuk <- sigbuk(d, datlist))
  if(is.null(Sbuk)){
    out <- -Inf
  } else{
    out <- -mean(diag(Sbuk))
  }
  return(out)
}

expentropgain <- function(d, datlist){
  sig2zsims <- datlist$sig2zsims
  thetasims <- datlist$thetasims
  betasims <- datlist$betasims
  errorsims <- datlist$errorsims
  nsim <- length(sig2zsims)
  N <- datlist$N
  s <- datlist$PredS
  M <- datlist$M
  Sinv <- datlist$Sinv
  b <- datlist$b
  covfun <- datlist$covfun
  d <- matrix(1/(1 + exp(-d)), ncol = 2)
  PredS <- datlist$PredS
  Xm <- cbind(1, PredS)
  Xn <- cbind(1, d)
  logZdens <- rep(0, nsim)
  logYZdens <- rep(0, nsim)
  Czs <- list()
  CZYs <- list()
  muZs <- list()
  muZYs <- list()
  Ysims <- matrix(0, nrow = nsim, ncol = M)
  Zsims <- matrix(0, nrow = nsim, ncol = N)
  ## construct the Y and Z sims from the error sims - depend on the design!
  ## note: don't have to redo the Y part every time though, bc fixed pred locations
  ## (can fix in future update; make sure this works first)
  out <- NULL
  try({
  for(i in 1:nsim){
    theta <- thetasims[i,]
    sig2z <- sig2zsims[i]
    beta <- betasims[i,]
    Cz <- Czfun(d, N, theta, sig2z, covfun)
    Czs[[i]] <- Cz
    Cyy <- Cyyfun(d, s, N, M, theta, covfun)
    Cy <- Cyfun(s, M, theta, covfun)
    CZY <- rbind(cbind(Cy, t(Cyy)), cbind(Cyy, Cz))
    CZYs[[i]] <- CZY
    muZ <- Xn%*%beta
    muZs[[i]] <- muZ
    muZY <- c(Xm%*%beta, muZ)
    muZYs[[i]] <- muZY
    YZsim <- muZY + crossprod(chol(CZY), drop(errorsims[i,]))
    Ysims[i,] <- YZsim[1:M]
    Zsims[i,] <- YZsim[M + 1:N]
  }
  ## now compute the expected log ratio of marginal densities
  for(i in 1:nsim){
    Z <- c(Zsims[i,])
    Y <- c(Ysims[i,])
    logYZconddens <- rep(0, nsim)
    logZconddens <- rep(0, nsim)
    for(j in 1:nsim){
      logYZconddens[j] <- dmnorm(c(Y, Z), muZYs[[j]], CZYs[[j]], TRUE)
      logZconddens[j] <- dmnorm(Z, drop(muZs[[j]]), Czs[[j]], TRUE)
    }
    Zconddens <- exp(logZconddens - max(logZconddens))
    YZconddens <- exp(logYZconddens - max(logYZconddens))
    logZdens[i] <- log(mean(Zconddens))
    logYZdens[i] <- log(mean(YZconddens))
  }
  out <-  mean(logYZdens) - mean(logZdens)})
  if(is.null(out)){
    out <- -Inf
  }
  return(out)
}


