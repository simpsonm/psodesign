theta1 <- .5
sig0 <- 1
sig1 <- 2
theta <- c(theta1, sig0, sig1)

expcov <- function(s1, s2, theta){
  sig0 <- theta[1]
  sig1 <- theta[2]
  theta1 <- theta[3]
  h <- s1 - s2
  euc <- sqrt(sum(h^2))
  out <- sig0*(euc == 0) + sig1*exp(-euc/theta1)
  return(out)
}

## num of design points
N <- 20

## num of points interested in
M <- 7*7

## X matrix for points we are interested in
grid <- seq(0, 1, length.out = sqrt(M)+2)
grid <- grid[-c(1, sqrt(M)+2)]
u <- rep(grid, sqrt(M))
v <- rep(grid, each = sqrt(M))
s <- cbind(u,v)
Xm <- cbind(1, u, v)

## pick some random design points
d <- cbind(runif(N), runif(N))
Xn <- cbind(1, d)

mspe <- function(x, datlist){
  ## create Cz
  N <- datlist$N
  s <- datlist$s
  M <- datlist$M
  theta <- datlist$theta
  covfun <- datlist$covfun
  d <- matrix(x, ncol = 2)
  if(min(d) < 0 | max(d) > 1){
    out <- Inf
  } else {
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
    Czinv <- chol2inv(chol(Cz))
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

source("../code/psofun.R")




niter <- 5000
nswarm <- 20
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496
nbhd <- sapply(1:nswarm, function(x){return( (x + -3:3 - 1)%%nswarm + 1)}) ## ring-3
init <- matrix(runif(nswarm*N*2, 0, 1), ncol = nswarm)
datlist <- list(N=N, s=s, M=M, theta=theta, covfun=expcov)

psotest <- pso(niter, nswarm, inertia, cognitive, social, init, nbhd, mspe, datlist = datlist)

psotest2 <- pso(niter, nswarm, 1, cognitive, social, init, nbhd, mspe, datlist = datlist,
               tune = TRUE)

psotest3 <- bbpso(niter, nswarm, 1, 0.3, init, nbhd, mspe, 1, TRUE, 0.5, datlist = datlist)

plot(ts(expcov(c(0,0), c(0,0), theta) - psotest$maxes))
lines(ts(expcov(c(0,0), c(0,0), theta) - psotest2$maxes), col = "red")
lines(ts(expcov(c(0,0), c(0,0), theta) - psotest3$maxes), col = "blue")


library(ggplot2)

plot1 <- qplot(u, v, data = data.frame(s), size = I(1), xlim=c(0,1), ylim=c(0,1)) +
  geom_point(aes(x=u, y=v), data = data.frame(u = psotest$argmax[1:N],
                                              v = psotest$argmax[1:N + N]),
             color = I("red"), size = I(2))

plot2 <- qplot(u, v, data = data.frame(s), size = I(1), xlim=c(0,1), ylim=c(0,1)) +
  geom_point(aes(x=u, y=v), data = data.frame(u = psotest2$argmax[1:N],
                                              v = psotest2$argmax[1:N + N]),
             color = I("blue"), size = I(2))

plot3 <- qplot(u, v, data = data.frame(s), size = I(1), xlim=c(0,1), ylim=c(0,1)) +
  geom_point(aes(x=u, y=v), data = data.frame(u = psotest3$argmax[1:N],
                                              v = psotest3$argmax[1:N + N]),
             color = I("green"), size = I(2))

library(gridExtra)

grid.arrange(plot1, plot2, plot3, ncol = 2)
