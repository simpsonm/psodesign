library(MCMCpack)
library(mnormt)
library(Matrix)

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

gelmanrwgibbs <- function(niter, init, datlist, tune = TRUE,
                          logrwsds = NULL,
                          H = 50, rwc = 0.1, rwtarget = c(0.4, 0.5)){
  nstate <- datlist$nstate
  nage <- datlist$nage
  nedu <- datlist$nedu
  nageedu <- nage*nedu
  nregion <- datlist$nregion
  if(is.null(logrwsds)){
   logrwsds <- rep(0, length(init) - 3 - 1 - nregion) ## alpha.region, beta.prev, variances
  }
  rwsds <- exp(logrwsds)
  par <- init
  polldat <- datlist$polldat
  statedat <- datlist$statedat
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  sig2a <- datlist$sig2a
  sig2b <- datlist$sig2b
  xzmat <- datlist$xzmat
  wmat <- datlist$wmat
  y <- polldat$bush
  female <- polldat$female
  black <- polldat$black
  age <- polldat$age
  edu <- polldat$edu
  ageedu <- cbind(age, edu)
  state <- polldat$state
  prev <- statedat$prev
  region <- statedat$region
  beta.0 <- par[1]
  beta.f <- par[2]
  beta.b <- par[3]
  beta.fb <- par[4]
  beta.prev <- par[5]
  betas <- par[1:4]
  alpha.state <- par[5 + 1:51]
  alpha.ageedu <- matrix(par[5 + 51 + 1:16], ncol = 4, nrow = 4)
  alpha.region <- par[5 + 51 + 16 + 1:5]
  alphas <- c(alpha.state, alpha.ageedu)
  alphabetas <- c(betas, alphas)
  lambdas <- c(beta.prev, alpha.region)
  sig2.state <- par[5 + 51 + 16 + 5 + 1]
  sig2.ageedu <- par[5 + 51 + 16 + 5 + 2]
  sig2.region <- par[5 + 51 + 16 + 5 + 3]
  sig2s <- c(sig2.state, sig2.ageedu, sig2.region)
  draws <- matrix(0, ncol = 80, nrow = niter)
  colnames(draws) <- c(paste("beta", c(0, "f", "b", "fb", "prev"), sep="."),
                       paste("alpha", "state", 1:nstate, sep = "."),
                       paste("alpha", "age", "edu", rep(1:nage, nedu),
                             rep(1:nedu, nage)[order(rep(1:nedu, nage))], sep="."),
                       paste("alpha", "region", 1:nregion, sep="."),
                       paste("sigma2", c("state", "ageedu", "region"), sep="."))
  accs <- matrix(0, ncol = length(alphabetas), nrow = niter)
  lambdabar <- c(betamn, rep(0, nregion))
  alpha.state.prop <- alpha.state
  alpha.ageedu.prop <- alpha.ageedu
  alpha.regionprop <- alpha.region
  for(iter in 1:niter){
    ## draw the variances
    sig2.state <- 1/rgamma(1, shape = sig2a + nstate/2,
                           rate = sig2b + crossprod(alpha.state - wmat%*%lambdas)/2)
    sig2.ageedu <- 1/rgamma(1, shape = sig2a + nageedu/2,
                             rate = sig2b + crossprod(c(alpha.ageedu))/2)
    sig2.region <- 1/rgamma(1, shape = sig2a + nregion/2,
                            rate = sig2b + crossprod(alpha.region)/2)
    ## draw lambda = (alpha.region, beta.prev)
    S <- diag(c(betavar, rep(sig2.region, nregion)))
    Sigma <- diag(rep(sig2.state, nstate))
    Sinv <- diag(1/diag(S))
    Sigmainv <- diag(1/diag(Sigma))
    Omegahat <- crossprod(wmat, Sigmainv)%*%wmat + Sinv
    chat <- crossprod(wmat, Sigmainv)%*%alpha.state + Sinv%*%lambdabar
    Sigmahat <- chol2inv(chol(Omegahat))
    muhat <- Sigmahat%*%chat
    cholhat <- chol(Sigmahat)
    lambdas <- drop(muhat + crossprod(cholhat, rnorm(nregion + 1)))
    alpha.region <- lambdas[1:nregion+1]
    beta.prev <- lambdas[1]
    ## random walk steps for all other alphas and betas ##
    lpold <- gelmanlpost(c(beta.0, beta.f, beta.b, beta.fb, alpha.state, alpha.ageedu, beta.prev,
                           alpha.region, log(sig2.state), log(sig2.ageedu), log(sig2.region)),
                         datlist)
    ## beta.0 rw step
    beta.0.prop <- beta.0 + rwsds[1]*rnorm(1)
    lpprop <- gelmanlpost(c(beta.0.prop, beta.f, beta.b, beta.fb, alpha.state, alpha.ageedu,
                           beta.prev,
                           alpha.region, log(sig2.state), log(sig2.ageedu), log(sig2.region)),
                         datlist)
    la <- lpprop - lpold
    u <- runif(1)
    if(la > log(u)){
      beta.0 <- beta.0.prop
      accs[iter,1] <- 1
      lpold <- lpprop
    }
    ## beta.female rw step
    beta.f.prop <- beta.f + rwsds[2]*rnorm(1)
    lpold <- gelmanlpost(c(beta.0, beta.f, beta.b, beta.fb, alpha.state, alpha.ageedu, beta.prev,
                           alpha.region, log(sig2.state), log(sig2.ageedu), log(sig2.region)),
                         datlist)
    lpprop <- gelmanlpost(c(beta.0, beta.f.prop, beta.b, beta.fb, alpha.state, alpha.ageedu,
                            beta.prev,
                            alpha.region, log(sig2.state), log(sig2.ageedu), log(sig2.region)),
                          datlist)
    la <- lpprop - lpold
    u <- runif(1)
    if(la > log(u)){
      beta.f <- beta.f.prop
      lpold <- lpprop
      accs[iter,2] <- 1
    }
    ## beta.black rw step
    beta.b.prop <- beta.b + rwsds[3]*rnorm(1)
    lpold <- gelmanlpost(c(beta.0, beta.f, beta.b, beta.fb, alpha.state, alpha.ageedu, beta.prev,
                           alpha.region, log(sig2.state), log(sig2.ageedu), log(sig2.region)),
                         datlist)    
    lpprop <- gelmanlpost(c(beta.0, beta.f, beta.b.prop, beta.fb, alpha.state, alpha.ageedu,
                            beta.prev,
                            alpha.region, log(sig2.state), log(sig2.ageedu), log(sig2.region)),
                          datlist)
    la <- lpprop - lpold
    u <- runif(1)
    if(la > log(u)){
      beta.b <- beta.b.prop
      lpold <- lpprop
      accs[iter,3] <- 1
    }
    ## beta.female.black rw step
    beta.fb.prop <- beta.fb + rwsds[4]*rnorm(1)
    lpold <- gelmanlpost(c(beta.0, beta.f, beta.b, beta.fb, alpha.state, alpha.ageedu, beta.prev,
                           alpha.region, log(sig2.state), log(sig2.ageedu), log(sig2.region)),
                         datlist)    
    lpprop <- gelmanlpost(c(beta.0, beta.f, beta.b, beta.fb.prop, alpha.state, alpha.ageedu,
                            beta.prev,
                            alpha.region, log(sig2.state), log(sig2.ageedu), log(sig2.region)),
                          datlist)
    la <- lpprop - lpold
    u <- runif(1)
    if(la > log(u)){
      beta.fb <- beta.fb.prop
      lpold <- lpprop
      accs[iter,4] <- 1
    }
    betas <- c(beta.0, beta.f, beta.b, beta.fb)
    ## alpha.state rw steps
    alpha.state.prop <- alpha.state
    for(i in 1:nstate){
      alpha.state.i.prop <- alpha.state[i] + rwsds[4 + i]*rnorm(1)
      alpha.state.prop[i] <- alpha.state.i.prop
      lpold <- gelmanlpost(c(beta.0, beta.f, beta.b, beta.fb, alpha.state, alpha.ageedu,
                             beta.prev,
                             alpha.region, log(sig2.state), log(sig2.ageedu), log(sig2.region)),
                           datlist)      
      lpprop <- gelmanlpost(c(beta.0, beta.f, beta.b, beta.fb, alpha.state.prop, alpha.ageedu,
                            beta.prev,
                            alpha.region, log(sig2.state), log(sig2.ageedu), log(sig2.region)),
                          datlist)
      la <- lpprop - lpold
      u <- runif(1)
      if(la > log(u)){
        alpha.state[i] <- alpha.state.i.prop
        lpold <- lpprop
        accs[iter,4 + i] <- 1
      } else {
        alpha.state.prop[i] <- alpha.state[i]
      }
    }
    ## alpha.ageedu rw steps
    for(i in 1:nage){
      for(j in 1:nedu){
        alpha.ageedu.ij.prop <- alpha.ageedu[i,j] +
          rwsds[4 + nstate + (j-1)*nage + i]*rnorm(1)
        alpha.ageedu.prop[i,j] <- alpha.ageedu.ij.prop
        lpold <- gelmanlpost(c(beta.0, beta.f, beta.b, beta.fb, alpha.state, alpha.ageedu,
                               beta.prev,
                               alpha.region, log(sig2.state), log(sig2.ageedu), log(sig2.region)),
                             datlist)        
        lpprop <- gelmanlpost(c(beta.0, beta.f, beta.b, beta.fb, alpha.state, alpha.ageedu.prop,
                                beta.prev, alpha.region,
                                log(sig2.state), log(sig2.ageedu), log(sig2.region)),
                              datlist)
        la <- lpprop - lpold
        u <- runif(1)
        if(la > log(u)){
          alpha.ageedu[i,j] <- alpha.ageedu.prop[i,j]
          lpold <- lpprop
          accs[iter,4 + nstate + (j-1)*nage + i] <- 1
        } else {
          alpha.ageedu.prop[i,j] <- alpha.ageedu[i,j]
        }
      }
    }
    draws[iter,] <- c(beta.0, beta.f, beta.b, beta.fb, beta.prev, alpha.state,
                      c(alpha.ageedu), alpha.region, sig2.state, sig2.ageedu, sig2.region)
    if(tune && iter %% H == 0){
      ##print(c(iter, "tuning!"))
      accrates <- apply(accs[(iter-H+1):iter,], 2, mean)
      logrwsds <- logrwsds + ((accrates > rwtarget[2]) - (accrates < rwtarget[1]))*rwc
    }
  }
  out <- list(draws = draws, accs = accs, logrwsds = logrwsds)
  return(out)
}

gelmanblockrwgibbs <- function(niter, init, datlist, tune = TRUE,
                               logrwsd = 0,
                               H = 50, rwc = 0.1, rwtarget = c(0.2, 0.3)){
  rwsd <- exp(logrwsd)
  par <- init
  polldat <- datlist$polldat
  statedat <- datlist$statedat
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  sig2a <- datlist$sig2a
  sig2b <- datlist$sig2b
  xzmat <- datlist$xzmat
  wmat <- datlist$wmat
  y <- polldat$bush
  female <- polldat$female
  black <- polldat$black
  age <- polldat$age
  edu <- polldat$edu
  ageedu <- cbind(age, edu)
  state <- polldat$state
  prev <- statedat$prev
  region <- statedat$region
  nstate <- datlist$nstate
  nage <- datlist$nage
  nedu <- datlist$nedu
  nageedu <- nage*nedu
  nregion <- datlist$nregion
  beta.0 <- par[1]
  beta.female <- par[2]
  beta.black <- par[3]
  beta.female.black <- par[4]
  beta.prev <- par[5]
  betas <- par[1:4]
  alpha.state <- par[5 + 1:51]
  alpha.ageedu <- matrix(par[5 + 51 + 1:16], ncol = 4, nrow = 4)
  alpha.region <- par[5 + 51 + 16 + 1:5]
  alphas <- c(alpha.state, alpha.ageedu)
  alphabetas <- c(betas, alphas)
  lambdas <- c(beta.prev, alpha.region)
  sig2.state <- par[5 + 51 + 16 + 5 + 1]
  sig2.ageedu <- par[5 + 51 + 16 + 5 + 2]
  sig2.region <- par[5 + 51 + 16 + 5 + 3]
  sig2s <- c(sig2.state, sig2.ageedu, sig2.region)
  draws <- matrix(0, ncol = 80, nrow = niter)
  colnames(draws) <- c(paste("beta", c(0, "f", "b", "fb", "prev"), sep="."),
                       paste("alpha", "state", 1:nstate, sep = "."),
                       paste("alpha", "age", "edu", rep(1:nage, nedu),
                             rep(1:nedu, nage)[order(rep(1:nedu, nage))], sep="."),
                       paste("alpha", "region", 1:nregion),
                       paste("sigma2", c("state", "ageedu", "region"), sep="."))
  accs <- rep(0, niter)
  lambdabar <- c(betamn, rep(0, nregion))
  muold <- xzmat%*%alphabetas
  lpdatold <- sum(y*muold - log(1 + exp(muold)))
  for(iter in 1:niter){
    ## draw the variances
    sig2.state <- 1/rgamma(1, shape = sig2a + nstate/2,
                           rate = sig2b + crossprod(alpha.state - wmat%*%lambdas)/2)
    sig2.ageedu <- 1/rgamma(1, shape = sig2a + nageedu/2,
                             rate = sig2b + crossprod(c(alpha.ageedu))/2)
    sig2.region <- 1/rgamma(1, shape = sig2a + nregion/2,
                            rate = sig2b + crossprod(alpha.region)/2)
    ## draw lambda = (alpha.region, beta.prev)
    S <- diag(c(rep(sig2.region, nregion), betavar))
    Sigma <- diag(rep(sig2.state, nstate))
    Sinv <- diag(1/diag(S))
    Sigmainv <- diag(1/diag(Sigma))
    Omegahat <- crossprod(wmat, Sigmainv)%*%wmat + Sinv
    chat <- crossprod(wmat, Sigmainv)%*%alpha.state + Sinv%*%lambdabar
    Sigmahat <- chol2inv(chol(Omegahat))
    muhat <- Sigmahat%*%chat
    cholhat <- chol(Sigmahat)
    lambdas <- drop(muhat + crossprod(cholhat, rnorm(nregion + 1)))
    alpha.region <- lambdas[1:nregion + 1]
    beta.prev <- lambdas[1]
    ## draw alphabeta
    alphabetasprop <- alphabetas + cholR%*%rnorm(nab, sd = rwsd)
    muprop <- xzmat%*%alphabetasprop
    lpdatprop <- sum(y*muprop - log(1 + exp(muprop)))
    lppriorold <- - sum((betas - betamn)^2)/(2*betavar) -
      sum((alpha.state - wmat*lambdas)^2)/(2*sig2.state) -
      sum(alpha.ageedu^2)/(2*sig2.ageedu)
    lppriorprop <- - sum((alphabetasprop[1:4] - betamn)^2)/(2*betavar) -
      sum((alphabetas[4 + 1:nstate] - wmat*lambdas)^2)/(2*sig2.state) -
      sum(alphabetas[4 + nstate + 1:nageedu]^2)/(2*sig2.ageedu)
    la <- lpdatprop + lppriorprop - lpdatold - lppriorold
    u <- runif(1)
    if(la > log(u)){
      accs[iter] <- 1
      alphabetas <- alphabetasprop
      betas <- alphabetas[1:4]
      alpha.state <- alphabetas[4 + 1:nstate]
      alpha.ageedu <- alphabetas[4 + nstate + 1:nageedu]
      lpdatold <- lpdatprop
    }
    draws[iter,] <- c(betas, beta.prev, alpha.state,
                    c(alpha.ageedu), alpha.region, sig2.state,
                    sig2.ageedu, sig2.region)
    if(tune && iter %% H == 0){
      ##print(c(iter, "tuning!"))
      accrate <- mean(accs[(iter-H+1):iter])
      logrwsd <- logrwsd + ((accrate > rwtarget[2]) - (accrate < rwtarget[1]))*rwc
    }
  }
  out <- list(draws = draws, accs = accs, logrwsd = logrwsd)
  return(out)
}


gelmanindwithingibbs <- function(niter, init, mu, df, datlist){
  nstate <- datlist$nstate
  nage <- datlist$nage
  nedu <- datlist$nedu
  nageedu <- nage*nedu
  nregion <- datlist$nregion
  nbeta <- datlist$nbeta
  par <- init
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  sig2a <- datlist$sig2a
  sig2b <- datlist$sig2b
  y <- datlist$y
  xmat <- datlist$xmat
  statemat <- datlist$statemat
  ageedumat <- datlist$agedumat
  regionmat <- datlist$regionmat
  prev <- datlist$prev
  beta.y <- par[1:nbeta]
  alpha.state <- par[nbeta + 1:nstate]
  alpha.ageedu <- par[nbeta + nstate + 1:nageedu]
  beta.prev <- par[nbeta + nstate + nageedu + 1]
  alpha.region <- par[nbeta + nstate + nageedu + 1 + 1:nregion]
  old <- c(beta.y, alpha.state, alpha.ageedu)
  lambdas <- c(beta.prev, alpha.region)
  wmat <- cbind(prev, regionmat)
  sig2.state <- par[nbeta + nstate + nageedu + 1 + nregion + 1]
  sig2.ageedu <- par[nbeta + nstate + nageedu + 1 + nregion + 2]
  sig2.region <- par[nbeta + nstate + nageedu + 1 + nregion + 3]
  sig2s <- c(sig2.state, sig2.ageedu, sig2.region)
  draws <- matrix(0, ncol = nbeta + nstate + nageedu + 1 + nregion + 3, nrow = niter)
  colnames(draws) <- c(paste("beta", 1:nbeta, sep="."),
                       paste("alpha", "state", 1:nstate, sep = "."),
                       paste("alpha", "age", "edu", rep(1:nage, nedu),
                             rep(1:nedu, nage)[order(rep(1:nedu, nage))], sep="."),
                       "beta.prev",
                       paste("alpha", "region", 1:nregion, sep="."),
                       paste("sigma2", c("state", "ageedu", "region"), sep="."))
  acc <- rep(0, niter)
  lambdabar <- c(betamn, rep(0, nregion))
  hess <- gelmanlposthess(mu, datlist)
  SigmaBig <- chol2inv(chol(-hess))
  Sigma1 <- SigmaBig[1:(nbeta + nstate + nageedu), 1:(nbeta + nstate + nageedu)]
  Sigma12 <- SigmaBig[1:(nbeta + nstate + nageedu), nbeta + nstate + nageedu + 1:(1 + nregion + 3)]
  Sigma22 <- SigmaBig[nbeta + nstate + nageedu + 1:(1 + nregion + 3),
                      nbeta + nstate + nageedu + 1:(1 + nregion + 3)]
  Sigma22inv <- chol2inv(chol(Sigma22))
  Sigma <- Sigma1 - tcrossprod(Sigma12%*%Sigma22inv, Sigma12)
  Rsigma <- chol(Sigma)
  Rprec <- chol(chol2inv(Rsigma))
  mu1 <- mu[1:(nbeta + nstate + nageedu)]
  mu2 <- mu[-c(1:(nbeta + nstate + nageedu))]
  for(iter in 1:niter){
    ## draw the variances
    sig2.state <- 1/rgamma(1, shape = sig2a + nstate/2,
                           rate = sig2b + crossprod(alpha.state - wmat%*%lambdas)/2)
    sig2.ageedu <- 1/rgamma(1, shape = sig2a + nageedu/2,
                             rate = sig2b + crossprod(c(alpha.ageedu))/2)
    sig2.region <- 1/rgamma(1, shape = sig2a + nregion/2,
                            rate = sig2b + crossprod(alpha.region)/2)
    ## draw lambda = (alpha.region, beta.prev)
    S <- diag(c(betavar, rep(sig2.region, nregion)))
    Sigma <- diag(rep(sig2.state, nstate))
    Sinv <- diag(1/diag(S))
    Sigmainv <- diag(1/diag(Sigma))
    Omegahat <- crossprod(wmat, Sigmainv)%*%wmat + Sinv
    chat <- crossprod(wmat, Sigmainv)%*%alpha.state + Sinv%*%lambdabar
    Sigmahat <- chol2inv(chol(Omegahat))
    muhat <- Sigmahat%*%chat
    cholhat <- chol(Sigmahat)
    lambdas <- drop(muhat + crossprod(cholhat, rnorm(nregion + 1)))
    alpha.region <- lambdas[1:nregion+1]
    beta.prev <- lambdas[1]
    otherdraw <- c(beta.prev, alpha.region, log(sig2.state), log(sig2.ageedu),
                   log(sig2.region))
    mn <- drop(mu1 + Sigma12%*%Sigma22inv%*%(otherdraw - mu2))
    prop <- rmtfixed(1, mn, Sigma, df, Rsigma)
    ljprop <- dmtcholprec(prop, mn, Rprec, df, TRUE)
    ljold <-  dmtcholprec(old,  mn, Rprec, df, TRUE)
    lpprop <- gelmanlpost(c(prop, otherdraw), datlist)
    lpold <- gelmanlpost(c(old, otherdraw), datlist)
    la <- lpprop - lpold + ljold - ljprop
    u <- runif(1)
    if(log(u) < la){
      old <- prop
      alpha.state <- old[nbeta + 1:nstate]
      alpha.ageedu <- old[nbeta + nstate + 1:nageedu]
      acc[iter] <- 1
    }
    draws[iter,] <- c(old, otherdraw)
  }    
  return(list(draws = draws, acc = acc, state = c(old, otherdraw)))
}

gelmanplusindwithingibbs <- function(niter, init, mu, df, datlist){
  nstate <- datlist$nstate
  nage <- datlist$nage
  nedu <- datlist$nedu
  nageedu <- nage*nedu
  nregion <- datlist$nregion
  npoll <- datlist$npoll
  nbeta <- datlist$nbeta
  par <- init
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  sig2a <- datlist$sig2a
  sig2b <- datlist$sig2b
  y <- datlist$y
  xmat <- datlist$xmat
  statemat <- datlist$statemat
  ageedumat <- datlist$ageedumat
  pollmat <- datlist$pollmat
  prev <- datlist$prev
  regionmat <- datlist$regionmat
  beta.y <- par[1:nbeta]
  alpha.state <- par[nbeta + 1:nstate]
  alpha.ageedu <- par[nbeta + nstate + 1:nageedu]
  alpha.poll <- par[nbeta + nstate + nageedu + 1:npoll]
  beta.prev <- par[nbeta + nstate + nageedu + npoll + 1]
  alpha.region <- par[nbeta + nstate + nageedu + npoll + 1 + 1:nregion]
  old <- c(beta.y, alpha.state, alpha.ageedu, alpha.poll)
  lambdas <- c(beta.prev, alpha.region)
  wmat <- cbind(prev, regionmat)
  sig2.state <-   par[nbeta + nstate + nageedu + npoll + 1 + nregion + 1]
  sig2.ageedu <- par[nbeta + nstate + nageedu + npoll + 1 + nregion + 2]
  sig2.poll <-    par[nbeta + nstate + nageedu + npoll + 1 + nregion + 3]
  sig2.region <-  par[nbeta + nstate + nageedu + npoll + 1 + nregion + 4]
  sig2s <- c(sig2.state, sig2.ageedu, sig2.poll, sig2.region)
  draws <- matrix(0, ncol = nbeta + nstate + nageedu + npoll + 1 + nregion + 4, nrow = niter)
  colnames(draws) <- c(paste("beta", 1:nbeta, sep="."),
                       paste("alpha", "state", 1:nstate, sep = "."),
                       paste("alpha", "age", "edu", rep(1:nage, nedu),
                             rep(1:nedu, nage)[order(rep(1:nedu, nage))], sep="."),
                       paste("alpha", "poll", 1:npoll, sep = "."),
                       "beta.prev",
                       paste("alpha", "region", 1:nregion, sep="."),
                       paste("sigma2", c("state", "ageedu", "poll", "region"), sep="."))
  acc <- rep(0, niter)
  lambdabar <- c(betamn, rep(0, nregion))
  hess <- gelmanpluslposthess(mu, datlist)
  SigmaBig <- chol2inv(chol(-hess))
  Sigma1 <- SigmaBig[1:(nbeta + nstate + nageedu + npoll), 1:(nbeta + nstate + nageedu + npoll)]
  Sigma12 <- SigmaBig[1:(nbeta + nstate + nageedu + npoll),
                      nbeta + nstate + nageedu + npoll + 1:(1 + nregion + 4)]
  Sigma22 <- SigmaBig[nbeta + nstate + nageedu + npoll + 1:(1 + nregion + 4),
                      nbeta + nstate + nageedu + npoll + 1:(1 + nregion + 4)]
  Sigma22inv <- chol2inv(chol(Sigma22))
  Sigma <- Sigma1 - tcrossprod(Sigma12%*%Sigma22inv, Sigma12)
  Rsigma <- chol(Sigma)
  Rprec <- chol(chol2inv(Rsigma))
  mu1 <- mu[1:(nbeta + nstate + nageedu + npoll)]
  mu2 <- mu[-c(1:(nbeta + nstate + nageedu + npoll))]
  for(iter in 1:niter){
    ## draw the variances
    sig2.state <- 1/rgamma(1, shape = sig2a + nstate/2,
                           rate = sig2b + crossprod(alpha.state - wmat%*%lambdas)/2)
    sig2.ageedu <- 1/rgamma(1, shape = sig2a + nageedu/2,
                             rate = sig2b + crossprod(c(alpha.ageedu))/2)
    sig2.poll <- 1/rgamma(1, shape = sig2a + npoll/2,
                             rate = sig2b + crossprod(c(alpha.poll))/2)
    sig2.region <- 1/rgamma(1, shape = sig2a + nregion/2,
                            rate = sig2b + crossprod(alpha.region)/2)
    ## draw lambda = (alpha.region, beta.prev)
    S <- diag(c(betavar, rep(sig2.region, nregion)))
    Sigma <- diag(rep(sig2.state, nstate))
    Sinv <- diag(1/diag(S))
    Sigmainv <- diag(1/diag(Sigma))
    Omegahat <- crossprod(wmat, Sigmainv)%*%wmat + Sinv
    chat <- crossprod(wmat, Sigmainv)%*%alpha.state + Sinv%*%lambdabar
    Sigmahat <- chol2inv(chol(Omegahat))
    muhat <- Sigmahat%*%chat
    cholhat <- chol(Sigmahat)
    lambdas <- drop(muhat + crossprod(cholhat, rnorm(nregion + 1)))
    alpha.region <- lambdas[1:nregion+1]
    beta.prev <- lambdas[1]
    otherdraw <- c(beta.prev, alpha.region, log(sig2.state), log(sig2.ageedu),
                   log(sig2.poll), log(sig2.region))
    mn <- drop(mu1 + Sigma12%*%Sigma22inv%*%(otherdraw - mu2))
    prop <- rmtfixed(1, mn, Sigma, df, Rsigma)
    ljprop <- dmtcholprec(prop, mn, Rprec, df, TRUE)
    ljold <-  dmtcholprec(old,  mn, Rprec, df, TRUE)
    lpprop <- gelmanpluslpost(c(prop, otherdraw), datlist)
    lpold <- gelmanpluslpost(c(old, otherdraw), datlist)
    la <- lpprop - lpold + ljold - ljprop
    u <- runif(1)
    if(log(u) < la){
      old <- prop
      alpha.state <- old[nbeta + 1:nstate]
      alpha.ageedu <- old[nbeta + nstate + 1:nageedu]
      alpha.poll <- old[nbeta + nstate + nageedu + 1:npoll]
      acc[iter] <- 1
    }
    draws[iter,] <- c(old, otherdraw)
  }    
  return(list(draws = draws, acc = acc, state = c(old, otherdraw)))
}
