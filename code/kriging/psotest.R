library(parallel)
library(doParallel)
library(foreach)
library(sp)
library(maptools)
library(rgeos)
source("../psofun.R")
source("krigingfun.R")
load("datlist.Rdata")
datlist$sppoly <- SpatialPolygons(list(b=Polygons(list(a=datlist$poly), "a")))

ncores <- 4
registerDoParallel(ncores)    #stopImplicitCluster()

nswarm <- 40
niter <- 10
nrep <- 4
inertias <- c(0.7298, 1/(log(2)*2))
cognitives <- c(1.496, log(2) + 1/2)
socials <- c(1.496, log(2) + 1/2)
nnbors <- c(1, 3, 40)
alpha <- 0.2*niter
beta <- 2
rates <- c(0.3, 0.5)
ccc <- 0.1
df <- 1
pcuts <- c(0, 0.5)
sig0 <- 1
inertia0 <- 1.2

ndesign <- 21
min.d <- apply(datlist$poly@coords, 2, min)
max.d <- apply(datlist$poly@coords, 2, max)
lower <- rep(min.d, each = ndesign)
upper <- rep(max.d, each = ndesign)

parsets <- 1:2
CFs <- c("CF", "nCF")
objnames <- c("sig2fuk.mean", "sig2fuk.max")

specs <- c(outer(c(outer(c(outer(paste(c("CI", "DI", "AT1", "AT2"), "PSO", sep="-"), parsets, paste, sep = "-")),
                        c(outer(CFs, nnbors, paste, sep="-")), paste, sep="-"),
                  c(outer(c(outer(paste(c("AT1", "AT2"), "BBPSO", sep="-"), parsets, paste, sep = "-")),
                          c(outer(CFs, nnbors, paste, sep="-")), paste, sep="-"))), objnames, paste, sep="-"))

seeds <- rep(3123, length(specs))

test <- psowrap(1, datlist, specs, seeds)

test <- foreach(i=c(1:40), .packages = c("sp", "maptools", "rgeos", "mnormt")) %dopar%
  psowrap(i, datlist, specs, seeds)

psowrap <- function(i, datlist, specs, seeds){
  spec <- specs[i]
  set.seed(seeds[i])
  splt <- strsplit(spec, "-")[[1]]
  style <- splt[1]
  alg <- splt[2]
  parset <- as.numeric(splt[3])
  CF <- splt[4]=="CF"
  nnbor <- as.numeric(splt[5])
  objname <- splt[6]
  obj <- switch(objname, sig2fuk.mean = sig2fuk.mean, sig2fuk.max = sig2fuk.max)
  repl <- 1
  if(alg == "BBPSO"){
    rate <- ifelse(style=="AT1", rates[1], rates[2])
    pcut <- pcuts[parset]
    init <- replicate(nswarm, c(spsample(datlist$poly, ndesign, "random")@coords))
    temp <- sbbpso(niter, nswarm, nnbor, sig0, obj, lower, upper,
                   pcut = pcut, CF = CF, AT = TRUE, rate = rate, df = df, ccc = 0.1,
                   init = init, boundaryfun = movetoboundary, datlist = datlist)
    algid <- paste("BBPSO", parset, ifelse(CF, "CF", "notCF"), style, sep = "-")
    tempdat <- data.frame(obj = objname, logpost = temp[["values"]],
                          time = time, algid = algid,
                          type = "BBPSO", parset = parset, CF = CF,
                          style = style, nbhd = nnbor, rep = repl,
                          inertias = temp$sigs)
    temppar <- data.frame(obj = objname, logpost = temp[["value"]],
                          algid = algid, type = "BBPSO", parset = parset, CF = CF,
                          style = style, nbhd = nnbor, rep = repl,
                          parid = 1:(ndesign*2), par = temp[["par"]])
  } else {
    rate <- ifelse(style=="AT1", rates[1], rates[2])
    c.in <- ifelse(style=="CI", inertias[parset], inertia0)
    c.co <- cognitives[parset]
    c.so <- socials[parset]
    init <- replicate(nswarm, c(spsample(datlist$poly, ndesign, "random")@coords))
    temp <- spso(niter, nswarm, nnbor, c.in, c.co, c.so, obj, lower, upper,
                 style = substr(style, 1, 2), CF = CF, alpha = alpha, beta = beta,
                 rate = rate, ccc = ccc,
                 init = init, boundaryfun = movetoboundary, datlist = datlist)
    algid <- paste("PSO", parset, ifelse(CF, "CF", "notCF"), style, sep = "-")
    tempdat <- data.frame(obj = objname, logpost = temp[["values"]],
                          time = time, algid = algid,
                          type = "PSO", parset = parset, CF = CF,
                          style = style, nbhd = nnbor, rep = repl,
                          inertias = temp$inertias)
    temppar <- data.frame(obj = objname, logpost = temp[["value"]],
                          algid = algid, type = "PSO", parset = parset, CF = CF,
                          style = style, nbhd = nnbor, rep = repl,
                          parid = 1:(ndesign*2), par = temp[["par"]])
  }
  out <- list(values = tempdat, pars = temppar)
  return(out)
}


time <- 0:niter

set.seed(3453)

nnbor <- 3
CF <- FALSE
parset <- 2
obj <- sig2fuk.mean
style <- "AT1"
rate <- ifelse(style=="AT1", rates[1], rates[2])
pcut <- pcuts[parset]


test <- foreach(nnbor = nnbors, .packages = c("sp", "maptools", "rgeos", "mnormt")) %dopar% {
  init <- replicate(nswarm, c(spsample(datlist$poly, ndesign, "random")@coords))
  sbbpso(niter, nswarm, nnbor, sig0, obj, lower, upper,
         pcut = pcut, CF = CF, AT = TRUE, rate = rate, df = df, ccc = 0.1,
         init = init, boundaryfun = movetoboundary, datlist = datlist)
}

style <- "CI"
parset <- 1
c.in <- ifelse(style=="CI", inertias[parset], inertia0)
c.co <- cognitives[parset]
c.so <- socials[parset]

init <- replicate(nswarm, c(spsample(datlist$poly, ndesign, "random")@coords))
cipso1 <- spso(niter, nswarm, nnbor, c.in, c.co, c.so, obj, lower, upper,
               style = substr(style, 1, 2), CF = CF, alpha = alpha, beta = beta,
               rate = rate, ccc = ccc,
               init = init, boundaryfun = movetoboundary, datlist = datlist)

init <- replicate(nswarm, c(spsample(datlist$poly, ndesign, "random")@coords))
atpso1 <- spso(niter, nswarm, nnbor, inertia0, c.co, c.so, obj, lower, upper,
               style = "AT", CF = CF, alpha = alpha, beta = beta,
               rate = rate, ccc = ccc,
               init = init, boundaryfun = movetoboundary, datlist = datlist)


c(at3bbpso$value, atpso1$value, cipso1$value)

plot(ts(at3bbpso$values), ylim = c(min(c(at3bbpso$values, atpso1$values, cipso1$values)),
                                   max(c(at3bbpso$values, atpso1$values, cipso1$values))))
lines(ts(atpso1$values), col = "red")
lines(ts(cipso1$values), col = "blue")

library(ggplot2)

testdat <- rbind(lower[c(1,ndesign + 1)],upper[c(1,ndesign + 1)])
testdat <- data.frame(cbind(rep(testdat[,1], each = 2), rep(testdat[,2], 2)))
testdat <- data.frame(rbind(testdat[1:2,], testdat[4,], testdat[3,]))
names(testdat) <- c("X", "Y")

ggplot(data = testdat, aes(x=X, y=Y)) + geom_line() +
  geom_point(data = data.frame(matrix(at3bbpso$par, ncol = 2)), aes(X1, X2), color = "red", size = I(4)) +
  geom_point(data = data.frame(matrix(atpso1$par, ncol = 2)), aes(X1, X2), color = "blue", size = I(4)) +
  geom_point(data = data.frame(matrix(cipso1$par, ncol = 2)), aes(X1, X2), color = "black", size = I(4))



dd <- matrix(at3bbpso, ncol = 2)
sig2fuk.mean(c(dd), datlist)

plot(datlist$poly@coords, type = "l")
points(dd)

dddist <- matrix(0, nrow(dd), nrow(dd))
for(i in 2:nrow(dd)){
  for(j in 1:(i-1)){
    dddist[i,j] <- sqrt(euc(dd[i,], dd[j,]))
    dddist[j,i] <- dddist[i,j]
  }
}

idx0 <- which(dddist == 0, TRUE)

dd[unique(idx0[which(idx0[,1] != idx0[,2]),][,1]),]

tt <- datlist$tt

vals <- NULL
idxs <- NULL
val <- 0
while(val < Inf){
  idx <- sample(1:1229, 30)
  val <- sig2fuk.mean(tt[idx,], datlist)
  idxs <- cbind(idxs, idx)
  vals <- c(vals, val)
}



ncand <- 2000
temp <- exch(ncand, obj, datlist$poly@coords, 10, 22, datlist = datlist)

dd <- matrix(at3bbpso, ncol = 2)
sig2fuk.mean(temp, datlist)


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
  DDDtemp <- Cz.d.s%*%tcrossprod(invCz.s, Cz.d.s)
  DDtemp <- Cz.d - (DDDtemp + t(DDDtemp))/2
  DD <- solve((DDtemp + t(DDtemp))/2)
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
  Rprec <- NULL
  try(Rprec <- chol(prec))
  if(is.null(Rprec))
    return(Inf)
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





## see if bbpso does better after lots of repetitions!!!

          for(style in c("CI", "DI", "AT1", "AT2")){
            rate <- ifelse(style=="AT1", rates[1], rates[2])
            algid <- paste("PSO", parset, ifelse(CF, "CF", "notCF"), style, sep = "-")
            tempdat <- data.frame(obj = objname, logpost = temp[["values"]],
                                  time = time, algid = algid,
                                  type = "PSO", parset = parset, CF = CF,
                                  style = style, nbhd = nnbor, rep = repl,
                                  inertias = temp$inertias)
            psoout <- rbind(psoout, tempdat)
            temppar <- data.frame(obj = objname, logpost = temp[["value"]],
                                  algid = algid, type = "PSO", parset = parset, CF = CF,
                                  style = style, nbhd = nnbor, rep = repl,
                                  parid = 1:(ndesign*2), par = temp[["par"]])
            parout <- rbind(parout, temppar)
          }
          write.csv(psoout, file = "psosimsout.csv", row.names=FALSE)
          write.csv(parout, file = "parsimsout.csv", row.names=FALSE)
        }
      }
    }
  }
}


set.seed(324280)

nbatches <- c(1,2)
nchrome <- 2
nrun <- ndesign
mutvars <- c(1,2)
mutrates <- c(1/100, 1/10)
nexnbors <- c(5, 10)
ncand <- 2000
parout2 <- NULL
gaout <- NULL
exout <- NULL
for(repl in 1:nrep){
  print("rep")
  print(repl)
  for(objnum in 1:2){
    print("obj")
    print(objnum)
    if(objnum == 1){
      obj <- sig2fuk.mean
      objname <- "sig2fuk.mean"
    } else {
      obj <- sig2fuk.max
      objname <- "sig2fuk.max"
    }
    for(nbatch in nbatches){
      for(mutvar in mutvars){
        for(mutrate in mutrates){
          print("GA")
          print(c(nbatch, mutvar, mutrate))
          temp <- ga(niter/nbatch, nbatch, floor(nswarm/2), nchrome, nrun, mutvar, mutrate,
                     lower, upper, obj, datlist=datlist)
          algid <- paste("GA", nbatch, mutrate, mutvar, sep="-")
          tempdat <- data.frame(obj = objname, logpost = temp[["values"]],
                                time = time, algid = algid, type = "GA",
                                nbatch = nbatch, mutrate = mutrate, mutvar = mutvar,
                                rep = repl)
          gaout <- rbind(gaout, tempdat)
          temppar <- data.frame(obj = objname, logpost = temp[["value"]],
                                algid = algid, type = "GA", parset = NA, CF = NA,
                                style = NA, nbhd = NA, rep = repl,
                                parid = 1:(ndesign*2), par = temp[["par"]])
          parout2 <- rbind(parout2, temppar)
        }
      }
    }
    for(nexnbor in nexnbors){
      print("EX")
      ### this needs to be rewritten because of exch not having a fixed niter
      temp <- exch(ncand, obj, datlist$poly@coords, nexnbor, ndesign, datlist = datlist)
      algid <- paste("EX", nexnbor, sep="-")
      tempdat <- data.frame(obj = objname, logpost = temp[["values"]],
                            objcount = temp$objcount,
                            time = 1:length(temp$values), algid = algid,
                            type = "EX", nnbor = nexnbor, rep = repl)
      exout <- rbind(exout, tempdat)
      temppar <- data.frame(obj = objname, logpost = temp[["value"]],
                            algid = algid, type = "EX", parset = NA, CF = NA,
                            style = NA, nbhd = NA, rep = repl,
                            parid = 1:(ndesign*2), par = temp[["par"]])
      parout2 <- rbind(parout2, temppar)
    }
    write.csv(exout, file = "exsimsout.csv", row.names=FALSE)
    write.csv(gaout, file = "gasimsout.csv", row.names=FALSE)
    write.csv(parout2, file = "parsimsout2.csv", row.names=FALSE)
  }
}








datlist2 <- datlist

dd0 <- spsample(datlist$sppoly, 10, "random")@coords
ss <- datlist$ss
dd2 <- c(rbind(dd0, ss))
dd1 <- c(dd0)
sig2fuk.mean(dd1, datlist)
sig2fuk.new.mean(dd2, datlist)
