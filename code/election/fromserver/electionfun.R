gelmanlpost <- function(par, datlist){
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  sig2a <- datlist$sig2a
  sig2b <- datlist$sig2b
  y <- datlist$y
  statemat <- datlist$statemat
  ageedumat <- datlist$ageedumat
  regionmat <- datlist$regionmat
  prev <- datlist$prev
  xmat <- datlist$xmat
  nstate <- datlist$nstate
  nage <- datlist$nage
  nedu <- datlist$nedu
  nageedu <- nage*nedu
  nregion <- datlist$nregion
  nbeta <- datlist$nbeta
  beta.y <- par[1:nbeta]
  alpha.state <-  par[nbeta + 1:nstate]
  alpha.ageedu <- par[nbeta + nstate + 1:nageedu]
  beta.prev <-    par[nbeta + nstate + nageedu + 1]
  betas <- c(par[1:nbeta], beta.prev)
  alpha.region <- par[nbeta + nstate + nageedu + 1 + 1:nregion]
  lsig2.state <-  par[nbeta + nstate + nageedu + 1 + nregion + 1]
  lsig2.ageedu <- par[nbeta + nstate + nageedu + 1 + nregion + 2]
  lsig2.region <- par[nbeta + nstate + nageedu + 1 + nregion + 3]
  sig2.state <- exp(lsig2.state)
  sig2.ageedu <- exp(lsig2.ageedu)
  sig2.region <- exp(lsig2.region)
  sig2s <- c(sig2.state, sig2.ageedu, sig2.region)
  lsig2s <- c(lsig2.state, lsig2.ageedu, lsig2.region)
  mu.y <- xmat%*%beta.y + ageedumat%*%alpha.ageedu + statemat%*%alpha.state
  log1minusp <- -log(1 + exp(mu.y))
  mu.state <- regionmat%*%alpha.region + prev*beta.prev
  datterm <- sum(y*mu.y) + sum(log1minusp)
  stateterm <- - sum((alpha.state - mu.state)^2)/(2*sig2.state) - nstate/2*lsig2.state
  randterms <- - sum(alpha.ageedu^2)/(2*sig2.ageedu) - nageedu/2*lsig2.ageedu +
      - sum(alpha.region^2)/(2*sig2.region) - nregion/2*lsig2.region
  priorandjacobianterms <- -sum((betas - betamn)^2)/(2*betavar) +
    sum(- sig2b/sig2s - sig2a*lsig2s)
  out <- datterm + stateterm + randterms + priorandjacobianterms
  return(out)
}

gelmanpluslpost <- function(par, datlist){
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  sig2a <- datlist$sig2a
  sig2b <- datlist$sig2b
  y <- datlist$y
  statemat <- datlist$statemat
  pollmat <- datlist$pollmat
  ageedumat <- datlist$ageedumat
  regionmat <- datlist$regionmat
  prev <- datlist$prev
  xmat <- datlist$xmat
  nstate <- datlist$nstate
  nage <- datlist$nage
  nedu <- datlist$nedu
  nageedu <- nage*nedu
  nregion <- datlist$nregion
  nbeta <- datlist$nbeta
  npoll <- datlist$npoll
  beta.y <- par[1:nbeta]
  alpha.state <-  par[nbeta + 1:nstate]
  alpha.ageedu <- par[nbeta + nstate + 1:nageedu]
  alpha.poll <-   par[nbeta + nstate + nageedu + 1:npoll]
  beta.prev <-    par[nbeta + nstate + nageedu + npoll + 1]
  betas <- c(par[1:4], beta.prev)
  alpha.region <- par[nbeta + nstate + nageedu + npoll + 1 + 1:nregion]
  lsig2.state <-  par[nbeta + nstate + nageedu + npoll + 1 + nregion + 1]
  lsig2.ageedu <- par[nbeta + nstate + nageedu + npoll + 1 + nregion + 2]
  lsig2.poll <-   par[nbeta + nstate + nageedu + npoll + 1 + nregion + 3]
  lsig2.region <- par[nbeta + nstate + nageedu + npoll + 1 + nregion + 4]
  sig2.state <- exp(lsig2.state)
  sig2.ageedu <- exp(lsig2.ageedu)
  sig2.poll <- exp(lsig2.poll)
  sig2.region <- exp(lsig2.region)
  sig2s <- c(sig2.state, sig2.ageedu, sig2.poll, sig2.region)
  lsig2s <- c(lsig2.state, lsig2.ageedu, lsig2.poll, lsig2.region)
  mu.y <- xmat%*%beta.y + ageedumat%*%alpha.ageedu + statemat%*%alpha.state +
    pollmat%*%alpha.poll
  log1minusp <- -log(1 + exp(mu.y))
  mu.state <- regionmat%*%alpha.region + prev*beta.prev
  datterm <- sum(y*mu.y) + sum(log1minusp)
  stateterm <- - sum((alpha.state - mu.state)^2)/(2*sig2.state) - nstate/2*lsig2.state
  randterms <- - sum(alpha.ageedu^2)/(2*sig2.ageedu) - nageedu/2*lsig2.ageedu +
    - sum(alpha.region^2)/(2*sig2.region) - nregion/2*lsig2.region +
    - sum(alpha.poll^2)/(2*sig2.poll) - npoll/2*lsig2.poll
  priorandjacobianterms <- -sum((betas - betamn)^2)/(2*betavar) +
    sum(- sig2b/sig2s - sig2a*lsig2s)
  out <- datterm + stateterm + randterms + priorandjacobianterms
  return(out)
}


gelmanlposthess <- function(par, datlist){
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  sig2a <- datlist$sig2a
  sig2b <- datlist$sig2b
  y <- datlist$y
  prev <- datlist$prev
  statemat <- datlist$statemat
  ageedumat <- datlist$ageedumat
  regionmat <- datlist$regionmat
  xmat <- datlist$xmat
  nstate <- datlist$nstate
  nage <- datlist$nage
  nedu <- datlist$nedu
  nageedu <- nage*nedu
  nregion <- datlist$nregion
  nbeta <- datlist$nbeta
  beta.y <- par[1:nbeta]
  alpha.state <-  par[nbeta + 1:nstate]
  alpha.ageedu <- par[nbeta + nstate + 1:nageedu]
  beta.prev <-    par[nbeta + nstate + nageedu + 1]
  alpha.region <- par[nbeta + nstate + nageedu + 1 + 1:nregion]
  lsig2.state <-  par[nbeta + nstate + nageedu + 1 + nregion + 1]
  lsig2.ageedu <- par[nbeta + nstate + nageedu + 1 + nregion + 2]
  lsig2.region <- par[nbeta + nstate + nageedu + 1 + nregion + 3]
  sig2.state <- exp(lsig2.state)
  sig2.ageedu <- exp(lsig2.ageedu)
  sig2.region <- exp(lsig2.region)
  sig2s <- c(sig2.state, sig2.ageedu, sig2.region)
  lsig2s <- c(lsig2.state, lsig2.ageedu, lsig2.region)
  mu.y <- xmat%*%beta.y + ageedumat%*%alpha.ageedu + statemat%*%alpha.state
  mu.state <- regionmat%*%alpha.region + prev*beta.prev
  dpdmu.y <- c(exp(mu.y)/(1 + exp(mu.y))^2)
  out <- matrix(0, nbeta + nstate + nageedu + 1 + nregion + 3,
                   nbeta + nstate + nageedu + 1 + nregion + 3)
  ## beta beta
  out[1:nbeta, 1:nbeta] <- - crossprod(xmat, dpdmu.y*xmat) - diag(1/betavar, nbeta)
  ## alpha.state alpha.state
  out[nbeta + 1:nstate, nbeta + 1:nstate] <-
    - crossprod(statemat, dpdmu.y*statemat) - diag(1/sig2.state, nstate)
  ## alpha.ageedu alpha.ageedu
  out[nbeta + nstate + 1:nageedu, nbeta + nstate + 1:nageedu] <-
    - crossprod(ageedumat, dpdmu.y*ageedumat) - diag(1/sig2.ageedu, nageedu)
  ## beta alpha.state
  out[1:nbeta, nbeta + 1:nstate] <- - crossprod(xmat, dpdmu.y*statemat)
  out[nbeta + 1:nstate, 1:nbeta] <- t(out[1:nbeta, nbeta + 1:nstate])
  ## beta alpha.ageedu
  out[1:nbeta, nbeta + nstate + 1:nageedu] <- - crossprod(xmat, dpdmu.y*ageedumat)
  out[nbeta + nstate + 1:nageedu, 1:nbeta] <- t(out[1:nbeta, nbeta + nstate + 1:nageedu])
  ## alpha.state alpha.ageedu
  out[nbeta + 1:nstate, nbeta + nstate + 1:nageedu] <- - crossprod(statemat, dpdmu.y*ageedumat)
  out[nbeta + nstate + 1:nageedu, nbeta + 1:nstate] <-
    t(out[nbeta + 1:nstate, nbeta + nstate + 1:nageedu])
  ## beta.prev
  out[nbeta + nstate + nageedu + 1, nbeta + nstate + nageedu + 1] <-
    - crossprod(prev)/sig2.state - 1/betavar
  ## alpha.region
  out[nbeta + nstate + nageedu + 1 + 1:nregion, nbeta + nstate + nageedu + 1 + 1:nregion] <-
    - crossprod(regionmat)/sig2.state - 1/sig2.region
  ## alpha.state alpha.region
  out[nbeta + 1:nstate, nbeta + nstate + nageedu + 1 + 1:nregion] <- regionmat/sig2.state
  out[nbeta + nstate + nageedu + 1 + 1:nregion, nbeta + 1:nstate] <- t(regionmat)/sig2.state
  ## beta.prev alpha.region
  out[nbeta + nstate + nageedu + 1, nbeta + nstate + nageedu + 1 + 1:nregion] <-
    - crossprod(prev, regionmat)/sig2.state
  out[nbeta + nstate + nageedu + 1 + 1:nregion, nbeta + nstate + nageedu + 1] <-
    - crossprod(regionmat, prev)/sig2.state
  ## alpha.state beta.prev
  out[nbeta + 1:nstate, nbeta + nstate + nageedu + 1] <- prev/sig2.state
  out[nbeta + nstate + nageedu + 1, nbeta + 1:nstate] <- prev/sig2.state
  ## tau.state tau.state
  out[nbeta + nstate + nageedu + 1 + nregion + 1, nbeta + nstate + nageedu + 1 + nregion + 1] <-
    - (sig2b + crossprod(alpha.state - regionmat%*%alpha.region - prev*beta.prev)/2)/sig2.state
  ## alpha.state tau.state
  out[nbeta + 1:nstate, nbeta + nstate + nageedu + 1 + nregion + 1] <-
    c(alpha.state - regionmat%*%alpha.region - prev*beta.prev)/sig2.state
  out[nbeta + nstate + nageedu + 1 + nregion + 1, nbeta + 1:nstate] <-
    c(alpha.state - regionmat%*%alpha.region - prev*beta.prev)/sig2.state
  ## alpha.region tau.state
  out[nbeta + nstate + nageedu + 1 + 1:nregion, nbeta + nstate + nageedu + 1 + nregion + 1] <-
    crossprod(regionmat, alpha.state - regionmat%*%alpha.region - prev*beta.prev)/sig2.state
  out[nbeta + nstate + nageedu + 1 + nregion + 1, nbeta + nstate + nageedu + 1 + 1:nregion] <-
    crossprod(alpha.state - regionmat%*%alpha.region - prev*beta.prev, regionmat)/sig2.state
  ## beta.prev tau.state
  out[nbeta + nstate + nageedu + 1 + nregion + 1, nbeta + nstate + nageedu + 1] <-
    crossprod(prev, alpha.state - regionmat%*%alpha.region - prev*beta.prev)/sig2.state
  out[nbeta + nstate + nageedu + 1, nbeta + nstate + nageedu + 1 + nregion + 1] <-
    crossprod(alpha.state - regionmat%*%alpha.region - prev*beta.prev, prev)/sig2.state
  ## tau.ageedu tau.ageedu
  out[nbeta + nstate + nageedu + 1 + nregion + 2, nbeta + nstate + nageedu + 1 + nregion + 2] <-
    - (sig2b + crossprod(alpha.ageedu)/2)/sig2.ageedu
  ## alpha.ageedu tau.age
  out[nbeta + nstate + 1:nageedu, nbeta + nstate + nageedu + 1 + nregion + 2] <-
    1/sig2.ageedu*alpha.ageedu
  out[nbeta + nstate + nageedu + 1 + nregion + 2, nbeta + nstate + 1:nageedu] <-
    1/sig2.ageedu*alpha.ageedu
  ## tau.region tau.region
  out[nbeta + nstate + nageedu + 1 + nregion + 3, nbeta + nstate + nageedu + 1 + nregion + 3] <-
    - (sig2b + crossprod(alpha.region)/2)/sig2.region
  ## alpha.region tau.region
  out[nbeta + nstate + nageedu + 1 + 1:nregion, nbeta + nstate + nageedu + 1 + nregion + 3] <-
    1/sig2.region*alpha.region
  out[nbeta + nstate + nageedu + 1 + nregion + 3, nbeta + nstate + nageedu + 1 + 1:nregion] <-
    1/sig2.region*alpha.region
  return(out)
}

gelmanpluslposthess <- function(par, datlist){
  betamn <- datlist$betamn
  betavar <- datlist$betavar
  sig2a <- datlist$sig2a
  sig2b <- datlist$sig2b
  y <- datlist$y
  prev <- datlist$prev
  statemat <- datlist$statemat
  ageedumat <- datlist$ageedumat
  pollmat <- datlist$pollmat
  regionmat <- datlist$regionmat
  xmat <- datlist$xmat
  nstate <- datlist$nstate
  nage <- datlist$nage
  nedu <- datlist$nedu
  nageedu <- nage*nedu
  nregion <- datlist$nregion
  nbeta <- datlist$nbeta
  npoll <- datlist$npoll
  beta.y <- par[1:nbeta]
  alpha.state <-  par[nbeta + 1:nstate]
  alpha.ageedu <- par[nbeta + nstate + 1:nageedu]
  alpha.poll <-   par[nbeta + nstate + nageedu + 1:npoll]
  beta.prev <-    par[nbeta + nstate + nageedu + npoll + 1]
  alpha.region <- par[nbeta + nstate + nageedu + npoll + 1 + 1:nregion]
  lsig2.state <-  par[nbeta + nstate + nageedu + npoll + 1 + nregion + 1]
  lsig2.ageedu <- par[nbeta + nstate + nageedu + npoll + 1 + nregion + 2]
  lsig2.poll <-   par[nbeta + nstate + nageedu + npoll + 1 + nregion + 3]
  lsig2.region <- par[nbeta + nstate + nageedu + npoll + 1 + nregion + 4]
  sig2.state <- exp(lsig2.state)
  sig2.ageedu <- exp(lsig2.ageedu)
  sig2.poll <- exp(lsig2.poll)
  sig2.region <- exp(lsig2.region)
  sig2s <- c(sig2.state, sig2.ageedu, sig2.poll, sig2.region)
  lsig2s <- c(lsig2.state, lsig2.ageedu, lsig2.poll, lsig2.region)
  mu.y <- xmat%*%beta.y + ageedumat%*%alpha.ageedu + statemat%*%alpha.state +
    pollmat%*%alpha.poll
  mu.state <- regionmat%*%alpha.region + prev*beta.prev
  dpdmu.y <- c(exp(mu.y)/(1 + exp(mu.y))^2)
  out <- matrix(0, nbeta + nstate + nageedu + npoll + 1 + nregion + 4,
                   nbeta + nstate + nageedu + npoll + 1 + nregion + 4)
  ## beta beta
  out[1:nbeta, 1:nbeta] <- - crossprod(xmat, dpdmu.y*xmat) - diag(1/betavar, nbeta)
  ## alpha.state alpha.state
  out[nbeta + 1:nstate, nbeta + 1:nstate] <-
    - crossprod(statemat, dpdmu.y*statemat) - diag(1/sig2.state, nstate)
  ## alpha.ageedu alpha.ageedu
  out[nbeta + nstate + 1:nageedu, nbeta + nstate + 1:nageedu] <-
    - crossprod(ageedumat, dpdmu.y*ageedumat) - diag(1/sig2.ageedu, nageedu)
  ## alpha.poll alpha.poll
  out[nbeta + nstate + nageedu + 1:npoll, nbeta + nstate + nageedu + 1:npoll] <-
    - crossprod(pollmat, dpdmu.y*pollmat) - diag(1/sig2.poll, npoll)
  ## beta alpha.state
  out[1:nbeta, nbeta + 1:nstate] <- - crossprod(xmat, dpdmu.y*statemat)
  out[nbeta + 1:nstate, 1:nbeta] <- t(out[1:nbeta, nbeta + 1:nstate])
  ## beta alpha.ageedu
  out[1:nbeta, nbeta + nstate + 1:nageedu] <- - crossprod(xmat, dpdmu.y*ageedumat)
  out[nbeta + nstate + 1:nageedu, 1:nbeta] <- t(out[1:nbeta, nbeta + nstate + 1:nageedu])
  ## alpha.state alpha.ageedu
  out[nbeta + 1:nstate, nbeta + nstate + 1:nageedu] <- - crossprod(statemat, dpdmu.y*ageedumat)
  out[nbeta + nstate + 1:nageedu, nbeta + 1:nstate] <-
    t(out[nbeta + 1:nstate, nbeta + nstate + 1:nageedu])
  ## beta alpha.poll
  out[1:nbeta, nbeta + nstate + nageedu + 1:npoll] <- - crossprod(xmat, dpdmu.y*pollmat)
  out[nbeta + nstate + nageedu + 1:npoll, 1:nbeta] <-
    t(out[1:nbeta, nbeta + nstate + nageedu + 1:npoll])
  ## alpha.state alpha.poll
  out[nbeta + 1:nstate, nbeta + nstate + nageedu + 1:npoll] <-
    - crossprod(statemat, dpdmu.y*pollmat)
  out[nbeta + nstate + nageedu + 1:npoll, nbeta + 1:nstate] <-
    t(out[nbeta + 1:nstate, nbeta + nstate + nageedu + 1:npoll])
  ## alpha.ageedu alpha.poll
  out[nbeta + nstate + 1:nageedu, nbeta + nstate + nageedu + 1:npoll] <-
    - crossprod(ageedumat, dpdmu.y*pollmat)
  out[nbeta + nstate + nageedu + 1:npoll, nbeta + nstate + 1:nageedu] <-
    t(out[nbeta + nstate + 1:nageedu, nbeta + nstate + nageedu + 1:npoll])
  ## beta.prev
  out[nbeta + nstate + nageedu + npoll + 1, nbeta + nstate + nageedu + npoll + 1] <-
    - crossprod(prev)/sig2.state - 1/betavar
  ## alpha.region
  out[nbeta + nstate + nageedu + npoll + 1 + 1:nregion,
      nbeta + nstate + nageedu + npoll + 1 + 1:nregion] <-
    - crossprod(regionmat)/sig2.state - diag(1/sig2.region, nregion)
  ## alpha.state alpha.region
  out[nbeta + 1:nstate, nbeta + nstate + nageedu + npoll + 1 + 1:nregion] <-
    regionmat/sig2.state
  out[nbeta + nstate + nageedu + npoll + 1 + 1:nregion, nbeta + 1:nstate] <-
    t(regionmat)/sig2.state
  ## beta.prev alpha.region
  out[nbeta + nstate + nageedu + npoll + 1, nbeta + nstate + nageedu + npoll + 1 + 1:nregion] <-
    - crossprod(prev, regionmat)/sig2.state
  out[nbeta + nstate + nageedu + npoll + 1 + 1:nregion, nbeta + nstate + nageedu + npoll + 1] <-
    - crossprod(regionmat, prev)/sig2.state
  ## alpha.state beta.prev
  out[nbeta + 1:nstate,            nbeta + nstate + nageedu + npoll + 1] <- prev/sig2.state
  out[nbeta + nstate + nageedu + npoll + 1, nbeta + 1:nstate           ] <- prev/sig2.state
  ## tau.state tau.state
  out[nbeta + nstate + nageedu + npoll + 1 + nregion + 1,
      nbeta + nstate + nageedu + npoll + 1 + nregion + 1] <-
    - (sig2b + crossprod(alpha.state - regionmat%*%alpha.region - prev*beta.prev)/2)/sig2.state
  ## alpha.state tau.state
  out[nbeta + 1:nstate, nbeta + nstate + nageedu + npoll + 1 + nregion + 1] <-
    c(alpha.state - regionmat%*%alpha.region - prev*beta.prev)/sig2.state
  out[nbeta + nstate + nageedu + npoll + 1 + nregion + 1, nbeta + 1:nstate] <-
    t(out[nbeta + 1:nstate, nbeta + nstate + nageedu + npoll + 1 + nregion + 1])
  ## alpha.region tau.state
  out[nbeta + nstate + nageedu + npoll + 1 + 1:nregion,
      nbeta + nstate + nageedu + npoll + 1 + nregion + 1] <-
    crossprod(regionmat, alpha.state - regionmat%*%alpha.region - prev*beta.prev)/sig2.state
  out[nbeta + nstate + nageedu + npoll + 1 + nregion + 1,
      nbeta + nstate + nageedu + npoll + 1 + 1:nregion  ] <-
    crossprod(alpha.state - regionmat%*%alpha.region - prev*beta.prev, regionmat)/sig2.state
  ## beta.prev tau.state
  out[nbeta + nstate + nageedu + npoll + 1 + nregion + 1,
      nbeta + nstate + nageedu + npoll + 1        ] <-
    crossprod(prev, alpha.state - regionmat%*%alpha.region - prev*beta.prev)/sig2.state
  out[nbeta + nstate + nageedu + npoll + 1, nbeta + nstate + nageedu + npoll + 1 + nregion + 1] <-
    crossprod(alpha.state - regionmat%*%alpha.region - prev*beta.prev, prev)/sig2.state
  ## tau.ageedu tau.ageedu
  out[nbeta + nstate + nageedu + npoll + 1 + nregion + 2,
      nbeta + nstate + nageedu + npoll + 1 + nregion + 2] <-
    - (sig2b + crossprod(alpha.ageedu)/2)/sig2.ageedu
  ## tau.poll tau.poll
  out[nbeta + nstate + nageedu + npoll + 1 + nregion + 3,
      nbeta + nstate + nageedu + npoll + 1 + nregion + 3] <-
    - (sig2b + crossprod(alpha.poll)/2)/sig2.poll
  ## alpha.poll tau.poll
  out[nbeta + nstate + 1:nageedu, nbeta + nstate + nageedu + npoll + 1 + nregion + 3] <-
    1/sig2.poll*alpha.poll
  out[nbeta + nstate + nageedu + npoll + 1 + nregion + 3, nbeta + nstate + 1:nageedu] <-
    1/sig2.poll*alpha.poll
  ## alpha.ageedu tau.age
  out[nbeta + nstate + 1:nageedu, nbeta + nstate + nageedu + npoll + 1 + nregion + 2] <-
    1/sig2.ageedu*alpha.ageedu
  out[nbeta + nstate + nageedu + npoll + 1 + nregion + 2, nbeta + nstate + 1:nageedu] <-
    1/sig2.ageedu*alpha.ageedu
  ## tau.region tau.region
  out[nbeta + nstate + nageedu + npoll + 1 + nregion + 4,
      nbeta + nstate + nageedu + npoll + 1 + nregion + 4] <-
    - (sig2b + crossprod(alpha.region)/2)/sig2.region
  ## alpha.region tau.region
  out[nbeta + nstate + nageedu + npoll + 1 + 1:nregion,
      nbeta + nstate + nageedu + npoll + 1 + nregion + 4] <- 1/sig2.region*alpha.region
  out[nbeta + nstate + nageedu + npoll + 1 + nregion + 4,
      nbeta + nstate + nageedu + npoll + 1 + 1:nregion  ] <- 1/sig2.region*alpha.region
  return(out)
}
