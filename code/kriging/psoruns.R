source("../psofun.R")

nrep <- 10
niter <- 1000
nswarm <- 50
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
N <- 20

## Number of points we're interested in
M <- 10^2

## Points where we predict - a regular grid
grid <- seq(0, 1, length.out = sqrt(M))
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
#  if(min(d) < 0 | max(d) > 1){
    out <- Inf
#  } else {
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
#  }
  return(-out)
}

datlist <- list(N=N, PredS=PredS, M=M, theta=theta, covfun=expcov)

npar <- N*2
inits <- list()
inits[[1]] <- matrix(runif(npar*nswarm, 0, 1), ncol = nswarm)
idxs <- replicate(nswarm, sample(1:M, N))
inits[[2]] <- rbind(matrix(PredS[idxs,1], ncol = nswarm), matrix(PredS[idxs,2], ncol = nswarm))
                    
for(m in 1:3){
  for(i in 1:2){
    for(rep in 1:nrep){
      cat("nbhd = ")
      cat(nbhdnames[m])
      cat(", init = ")
      cat(i)
      cat(", rep = ")
      cat(rep)
      cat("\n")
      print("PSO")
      psotemp <- pso(niter, nswarm, inertia, cognitive, social, inits[[i]], nbhd[[m]],
                     negmspe, datlist = datlist)
      psotempout <- data.frame(algorithm = "PSO", init = i, nbhd = nbhdnames[m], rep = rep,
                             iteration = 0:niter, maxes = psotemp$maxes)
      psoout <- rbind(psoout, psotempout)
      print("BBPSO-MC")
      psotemp <- bbpso(niter, nswarm, 0, 1, inits[[i]], nbhd[[m]], negmspe, Inf, FALSE,
                       0, datlist = datlist)
      psotempout <- data.frame(algorithm = "BBPSO-MC", init = i, nbhd = nbhdnames[m], rep = rep,
                             iteration = 0:niter, maxes = psotemp$maxes)
      psoout <- rbind(psoout, psotempout)            
      print("BBPSOxp-MC")
      psotemp <- bbpso(niter, nswarm, 0, 1, inits[[i]], nbhd[[m]], negmspe, Inf, FALSE,
                       .5, datlist = datlist)
      psotempout <- data.frame(algorithm = "BBPSOxp-MC", init = i, nbhd = nbhdnames[m], rep = rep,
                             iteration = 0:niter, maxes = psotemp$maxes)
      psoout <- rbind(psoout, psotempout)
      print("DI-PSO")
      psotemp <- pso(niter, nswarm, inertia, social, cognitive, inits[[i]], nbhd[[m]],
                     negmspe, datlist = datlist, tune = TRUE, style = "deterministic",
                     alpha = alpha, beta = beta)
      psotempout <- data.frame(algorithm = "DI-PSO", init = i, nbhd = nbhdnames[m], rep = rep,
                             iteration = 0:niter, maxes = psotemp$maxes)
      psoout <- rbind(psoout, psotempout)
      print("AT-PSO & AT-BBPSO-MC & AT-BBPSOxp-MC")
      for(rate in rates){
        psotemp <- pso(niter, nswarm, 0.9, cognitive, social, inits[[i]], nbhd[[m]],
                       negmspe, datlist = datlist, tune = TRUE, style = "adaptive",
                       rate = rate, ccc = ccc)
        psotempout <- data.frame(algorithm = paste("AT-PSO", rate, ccc, sep="-"),
                                 init = i, nbhd = nbhdnames[m], rep = rep,
                                 iteration = 0:niter, maxes = psotemp$maxes)
        psoout <- rbind(psoout, psotempout)
        for(df in dfs){
          print(paste(c(rate, df)))
          psoout <- rbind(psoout, psotempout)
          psotemp <- bbpso(niter, nswarm, 1, rate, inits[[i]], nbhd[[m]], negmspe, df,
                           TRUE, 0, datlist = datlist, ccc = ccc)
          psotempout <- data.frame(algorithm = paste("AT-BBPSO-MC", df, rate, ccc, sep="-"),
                                   init = i, nbhd = nbhdnames[m], rep = rep,
                                   iteration = 0:niter, maxes = psotemp$maxes)
          psoout <- rbind(psoout, psotempout)
          psotemp <- bbpso(niter, nswarm, 1, rate, inits[[i]], nbhd[[m]], negmspe, df,
                           TRUE, 0.5, datlist = datlist, ccc = ccc)
          psotempout <- data.frame(algorithm = paste("AT-BBPSOxp-MC", df, rate, ccc, sep="-"),
                                   init = i, nbhd = nbhdnames[m], rep = rep,
                                   iteration = 0:niter, maxes = psotemp$maxes)
          psoout <- rbind(psoout, psotempout)
        }
      }
      write.csv(psoout, file = "psoout.csv", row.names=FALSE)
    }
  }
}
