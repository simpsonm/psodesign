logit <- function(x) return(log(x/(1-x)))

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
  out <- sig0*(euc == 0) + sig1*exp(-(euc/theta1))
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
  out <-  mean(logYZdens) - mean(logZdens)
  return(out)
}


