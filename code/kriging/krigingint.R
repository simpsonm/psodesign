library(ggplot2)
library(gridExtra)
library(fields) ## for space covering design
library(cubature)

source("../psofun.R")
source("krigingfun.R")

niter <- 100
nswarm <- 20
nrep <- 1
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496
nbhdnames <- c("ring-1", "ring-3", "global")
rates <- c(0.5)
dfs <- c(5)
ccc <- c(0.1)
psoout <- NULL
nbhd <- list()
nbhd[[1]] <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)}) ## ring-1
nbhd[[2]] <- sapply(1:nswarm, function(x){return( (x + -3:3 - 1)%%nswarm + 1)}) ## ring-3
nbhd[[3]] <- matrix(1:nswarm, ncol=nswarm, nrow=nswarm) ## global

theta1 <- 1
theta2 <- 1
sig0 <- 0.1
sig1 <- 0.9
sigz <- .1
theta <- c(theta1, sig0, sig1)

## Number of design points
N <- 10

## Number of points we're interested in
M <- 10^2

## Points where we predict - a regular grid
grid <- seq(0, 1, length.out = sqrt(M) + 2)
grid <- grid[-c(1,sqrt(M) + 2)]
PredS <- cbind(u = rep(grid, sqrt(M)), v = rep(grid, each = sqrt(M)))
S <- 100*diag(3)
Sinv <- solve(S)

datlist <- list(N=N, PredS=PredS, M=M, theta=theta, covfun=expcov, Sinv=Sinv, sig2z=sigz, S=S, b=rep(0,3), asig = 1, bsig = 1)

phi <- c(sigz, theta)
Ystar <- rnorm(M)
Z <- rnorm(N)
d <- logit(runif(2*N))

yzdenstest <- yzdens(phi, Ystar, Z, d, datlist)

system.time(inttest <- adaptIntegrate(yzdens, rep(0, 2), rep(10, 2), Ystar = Ystar, Z = Z, d = d, datlist = datlist))

yzgivenphidens <- function(Ystar, Z, Yhat, Shat, Zmu, Zprec){
  Ny <- length(Ystar)
  Nz <- length(Z)
  RY <- chol(Shat)
  RYinv <- backsolve(RY, diag(Ny))
  RZ <- chol(Zprec)
  logypart <- -tcrossprod(crossprod(Ystar - Yhat, RYinv))/2 - sum(log(diag(RY))) - Ny/2*log(2*pi)
  logzpart <- -crossprod(RZ%*%(Z - Zmu))/2 + sum(log(diag(RZ))) - Nz/2*log(2*pi)
  out <- logypart + logzpart
  return(out)
}

phidens <- function(phi, datlist){
  asig <- datlist$asig
  bsig <- datlist$bsig
  out <- sum(dgamma(1/phi, asig, bsig, log = TRUE)) - 2*sum(log(phi))
  return(out)
}

yzdens <- function(phi, Ystar, Z, d, datlist){
  N <- datlist$N
  s <- datlist$PredS
  M <- datlist$M
  Sinv <- datlist$Sinv
  S <- datlist$S
  phistar <- c(.1, 1, phi)
  theta <- phistar[-1]
  sig2z <- phistar[1]
  b <- datlist$b
  covfun <- datlist$covfun
  d <- matrix(1/(1 + exp(-d)), ncol = 2)
  Xm <- cbind(1, PredS)
  Xn <- cbind(1, d)
  ## create C_{Z}
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
  ## create C_{YY^*}
  Cyy <- matrix(0, N, M)
  for(i in 1:N){
    for(j in 1:M){
      Cyy[i,j] <- covfun(d[i,], s[j,], theta)
    }
  }
  ## create C_{Y^*}
  Cy <- matrix(0, M, M)
  for(i in 1:M){
    for(j in 1:i){
      Cy[i,j] <- covfun(s[i,], s[j,], theta)
      if(j < i){
        Cy[j,i] <- Cy[i,j]
      }
    }
  }
  R <- chol(Cz)
  Rinv <- backsolve(R, diag(N))
  tRinvXn <- crossprod(Rinv, Xn)
  tRinvCyy <- crossprod(Rinv, Cyy)
  Pgls <- crossprod(tRinvXn) + Sinv
  delta <- Xm - crossprod(tRinvCyy, tRinvXn)
  RP <- chol(Pgls)
  RPinv <- backsolve(RP, diag(3))
  Sbuk <- Cy - crossprod(tRinvCyy) + tcrossprod(delta%*%RPinv)
  betabgls <- crossprod(RPinv)%*%(Sinv%*%b +  crossprod(tRinvXn, crossprod(Rinv, Z)))
  Ybuk <- Xm%*%betabgls + crossprod(tRinvCyy, crossprod(Rinv, Z - Xn%*%betabgls))
  Zprec <- Cz + Xn%*%tcrossprod(S, Xn)
  out1 <- yzgivenphidens(Ystar, Z, Ybuk, Sbuk, Xn%*%b, Zprec)
  out2 <- phidens(phi, datlist)
  return(exp(out1 + out2))
}

