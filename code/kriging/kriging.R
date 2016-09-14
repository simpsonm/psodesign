library(ggplot2)
library(gridExtra)
library(fields) ## for space covering design

source("../psofun.R")

nrep <- 1
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496
nbhdnames <- c("ring-1", "ring-3", "global")
rates <- c(0.5)
dfs <- c(5)
ccc <- c(0.1)
alpha <- .2*niter
beta <- 1
psoout <- NULL
nbhd <- list()
nbhd[[1]] <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)}) ## ring-1
nbhd[[2]] <- sapply(1:nswarm, function(x){return( (x + -3:3 - 1)%%nswarm + 1)}) ## ring-3
nbhd[[3]] <- matrix(1:nswarm, ncol=nswarm, nrow=nswarm) ## global

theta1 <- 1
theta2 <- 1
sig0 <- 0.1
sig1 <- 0.9
theta <- c(theta1, theta2, sig0, sig1)

## powered exponential covariance function
expcov <- function(s1, s2, theta){
  theta1 <- theta[1]
  theta2 <- theta[2]
  sig0 <- theta[3]
  sig1 <- theta[4]
  h <- s1 - s2
  euc <- sqrt(sum(h^2))
  out <- sig0*(euc == 0) + sig1*exp(-(euc/theta1)^theta2)
  return(out)
}

## Number of design points
N <- 10

## Number of points we're interested in
M <- 10^2

## Points where we predict - a regular grid
grid <- seq(0, 1, length.out = sqrt(M) + 2)
grid <- grid[-c(1,sqrt(M) + 2)]
PredS <- cbind(u = rep(grid, sqrt(M)), v = rep(grid, each = sqrt(M)))

## negative MSPE, ignoring the C_Y(s,s) term
negmspe <- function(x, datlist){
  ## create Cz
  N <- datlist$N
  s <- datlist$PredS
  M <- datlist$M
  theta <- datlist$theta
  covfun <- datlist$covfun
  d <- matrix(x, ncol = 2)
  Xm <- cbind(1, PredS)
  ##  if(min(d) < 0 | max(d) > 1){
  out <- Inf
  ##  } else {
  Cz <- matrix(0, N, N)
  for(i in 1:N){
    for(j in 1:i){
      Cz[i,j] <- covfun(d[i,], d[j,], theta)
      if(j < i){
        Cz[j,i] <- Cz[i,j]
      }
    }
  }
  ## create Cy
  Cy <- matrix(0, M, N)
  for(i in 1:M){
    for(j in 1:N){
      Cy[i,j] <- covfun(s[i,], d[j,], theta)
    }
  }
  Czinv <- NULL
  try(Czinv <- chol2inv(chol(Cz)))
  if(is.null(Czinv)){
    out <- Inf
  } else {
    Xn <- cbind(1, d)
    XprimeCinv <- crossprod(Xn, Czinv)
    precmat <- crossprod(Xn, Czinv)%*%Xn
    covmat <- chol2inv(chol(precmat))
    sigyhats <- rep(0, M)
    for(i in 1:M){
      delta <- Xm[i,] - XprimeCinv%*%Cy[i,]
      sigyhats[i] <- crossprod(delta, covmat)%*%delta - crossprod(Cy[i,], Czinv)%*%Cy[i,]
    }
    out <- mean(sigyhats)
  }
  ##  }
  return(-out)
}

negmspe2 <- function(x, d, datlist){
  ## create Cz
  N <- datlist$N
  s <- datlist$PredS
  M <- datlist$M
  theta <- x
  covfun <- datlist$covfun
  d <- matrix(1/(1 + exp(-d)), ncol = 2)
  Xm <- cbind(1, PredS)
  Cz <- matrix(0, N, N)
  for(i in 1:N){
    for(j in 1:i){
      Cz[i,j] <- covfun(d[i,], d[j,], theta)
      if(j < i){
        Cz[j,i] <- Cz[i,j]
      }
    }
  }
  ## create Cy
  Cy <- matrix(0, M, N)
  for(i in 1:M){
    for(j in 1:N){
      Cy[i,j] <- covfun(s[i,], d[j,], theta)
    }
  }
  Czinv <- NULL
  try(Czinv <- chol2inv(chol(Cz)))
  if(is.null(Czinv)){
    out <- Inf
  } else {
    Xn <- cbind(1, d)
    XprimeCinv <- crossprod(Xn, Czinv)
    precmat <- crossprod(Xn, Czinv)%*%Xn
    covmat <- chol2inv(chol(precmat))
    sigyhats <- rep(0, M)
    for(i in 1:M){
      delta <- Xm[i,] - XprimeCinv%*%Cy[i,]
      sigyhats[i] <- crossprod(delta, covmat)%*%delta - crossprod(Cy[i,], Czinv)%*%Cy[i,]
    }
    out <- mean(sigyhats)
  }
return(-out)
}

exnegmspe <- function(x, datlist){
  sims <- datlist$sims
  N <- datlist$N
  s <- datlist$PredS
  M <- datlist$M
  covfun <- datlist$covfun
  nsims <- nrow(sims)
  d <- matrix(1/(1 + exp(-x)), ncol = 2)
  out <- rep(0, nsims)
  Xm <- cbind(1, PredS)
  Cz <- matrix(0, N, N)
  Cy <- matrix(0, M, N)
  for(iter in 1:nsims){
    theta <- sims[iter,]
    for(i in 1:N){
      for(j in 1:i){
        Cz[i,j] <- covfun(d[i,], d[j,], theta)
        if(j < i){
          Cz[j,i] <- Cz[i,j]
        }
      }
    }
    for(i in 1:M){
      for(j in 1:N){
        Cy[i,j] <- covfun(s[i,], d[j,], theta)
      }
    }
    Czinv <- NULL
    try(Czinv <- chol2inv(chol(Cz)))
    if(is.null(Czinv)){
      out[iter] <- Inf
    } else {
      Xn <- cbind(1, d)
      XprimeCinv <- crossprod(Xn, Czinv)
      precmat <- crossprod(Xn, Czinv)%*%Xn
      covmat <- NULL
      covmat <- chol2inv(chol(precmat))
      if(is.null(covmat)){
        out[iter] <- Inf
      } else {
        sigyhats <- rep(0, M)
        for(i in 1:M){
          delta <- Xm[i,] - XprimeCinv%*%Cy[i,]
          sigyhats[i] <- crossprod(delta, covmat)%*%delta - crossprod(Cy[i,], Czinv)%*%Cy[i,]
        }
        out[iter] <- mean(sigyhats)
      }
    }
  }
  return(-mean(out))
}

logit <- function(x) return(log(x/(1-x)))

niter <- 1000
nswarm <- 20
npar <- N*2
inits <- list()
cands <- PredS[PredS[,1]> 0 & PredS[,1]<1 & PredS[,2]> 0 & PredS[,2]<1,]
spacefill <- cover.design(cands, N)
inits[[1]] <- logit(matrix(runif(npar*nswarm, 0, 1), ncol = nswarm))
idxs <- replicate(nswarm, sample(1:nrow(cands), N))
inits[[2]] <- logit(rbind(matrix(cands[idxs,1], ncol = nswarm), matrix(cands[idxs,2], ncol = nswarm)))
inits[[2]][,1] <- logit(c(spacefill$design))
nsim <- 1
aaa <- 10
sims <- cbind(1/rgamma(nsim, aaa, (aaa-1)*1),  1/rgamma(nsim, aaa, (aaa-1)*1),
              1/rgamma(nsim, aaa, (aaa-1)*.1), 1/rgamma(nsim, aaa, (aaa-1)*.9))
##sims <- matrix(c(1, 1, .1, .9), nrow = 1)
datlist <- list(N=N, PredS=PredS, M=M, theta=theta, covfun=expcov, sims=sims)
nbhd <- list()
nbhd[[1]] <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)}) ## ring-1
nbhd[[2]] <- sapply(1:nswarm, function(x){return( (x + -3:3 - 1)%%nswarm + 1)}) ## ring-3
nbhd[[3]] <- matrix(1:nswarm, ncol=nswarm, nrow=nswarm) ## global



psoout1 <- pso(niter, nswarm, inertia, cognitive, social, inits[[1]], nbhd[[1]], exnegmspe,
               datlist = datlist)

psoout2 <- pso(niter, nswarm, inertia, cognitive, social, inits[[2]], nbhd[[1]], exnegmspe,
               datlist = datlist)

outdf <- data.frame(maxes=c(psoout1$maxes, psoout2$maxes), iter = rep(0:niter, 2), init = rep(c("A", "B"), each = niter + 1))

p0 <- qplot(iter, -maxes, color = init, data = subset(outdf, iter > 0), geom="line")
p00 <- qplot(iter, -maxes, color = init, data = subset(outdf, iter > 10), geom="line")
psodat1 <- data.frame(u = 1/(1 + exp(-psoout1$argmax[1:N])),
                      v = 1/(1 + exp(-psoout1$argmax[1:N + N])))
p1 <- qplot(u, v, data = data.frame(PredS), size = I(1), xlim=c(0,1), ylim=c(0,1)) +
  geom_point(aes(u, v), data = psodat1, color = I("red"), size = I(3))
psodat2 <- data.frame(u = 1/(1 + exp(-psoout2$argmax[1:N])),
                      v = 1/(1 + exp(-psoout2$argmax[1:N + N])))
p2 <- qplot(u, v, data = data.frame(PredS), size = I(1), xlim=c(0,1), ylim=c(0,1)) +
  geom_point(aes(u, v), data = psodat2, color = I("blue"), size = I(3))

grid.arrange(p1, p2, p0, p00, ncol = 2)


