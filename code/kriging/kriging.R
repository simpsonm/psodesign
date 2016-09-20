library(ggplot2)
library(gridExtra)
library(fields) ## for space covering design

source("../psofun.R")
source("krigingfun.R")

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

theta1 <- .0010
theta2 <- .0010
sig0 <- 0.1
sig1 <- 0.9
sigz <- .1
theta <- c(theta1, theta2, sig0, sig1)

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

nsim <- 100
thetasims <- matrix(1/rgamma(nsim*3, 10, 10), ncol = 3)
sig2zsims <- 1/rgamma(nsim, 10, 10)
betasims <- matrix(rnorm(3*nsim, 0, 10), ncol=3)
errorsims <- matrix(rnorm(nsim*(N + M)), ncol = N + M)

datlist <- list(N=N, PredS=PredS, M=M, theta=theta, covfun=expcov, Sinv=Sinv, sig2z=sigz, S=S,
                b=rep(0,3), sig2zsims=sig2zsims, thetasims=thetasims, betasims=betasims,
                errorsims=errorsims)

niter <- 100
nswarm <- 20
nbhd <- list()
nbhd[[1]] <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)}) ## ring-1
nbhd[[2]] <- sapply(1:nswarm, function(x){return( (x + -3:3 - 1)%%nswarm + 1)}) ## ring-3
nbhd[[3]] <- matrix(1:nswarm, ncol=nswarm, nrow=nswarm) ## global

npar <- N*2
inits <- list()
cands <- PredS[PredS[,1]> 0 & PredS[,1]<1 & PredS[,2]> 0 & PredS[,2]<1,]
spacefill <- cover.design(cands, N)
inits[[1]] <- logit(matrix(runif(npar*nswarm, 0, 1), ncol = nswarm))
idxs <- replicate(nswarm, sample(1:nrow(cands), N))
inits[[2]] <- logit(rbind(matrix(cands[idxs,1], ncol = nswarm), matrix(cands[idxs,2], ncol = nswarm)))
inits[[2]][,1] <- logit(c(spacefill$design))



psoout1 <- pso(niter, nswarm, inertia, cognitive, social, inits[[1]], nbhd[[1]], expentropgain,
               datlist = datlist)

psoout2 <- pso(niter, nswarm, inertia, cognitive, social, inits[[2]], nbhd[[1]], expentropgain,
               datlist = datlist)

outdf <- data.frame(maxes=c(psoout1$maxes, psoout2$maxes), iter = rep(0:niter, 2), init = rep(c("A", "B"), each = niter + 1))

p0 <- qplot(iter, maxes, color = init, data = subset(outdf, init %in% c("A", "B")), geom="line")
psodat1 <- data.frame(u = 1/(1 + exp(-psoout1$argmax[1:N])),
                      v = 1/(1 + exp(-psoout1$argmax[1:N + N])))
p1 <- qplot(u, v, data = data.frame(PredS), size = I(1), xlim=c(0,1), ylim=c(0,1)) +
  geom_point(aes(u, v), data = psodat1, color = I("red"), size = I(3))
psodat2 <- data.frame(u = 1/(1 + exp(-psoout2$argmax[1:N])),
                      v = 1/(1 + exp(-psoout2$argmax[1:N + N])))
p2 <- qplot(u, v, data = data.frame(PredS), size = I(1), xlim=c(0,1), ylim=c(0,1)) +
  geom_point(aes(u, v), data = psodat2, color = I("blue"), size = I(3))

grid.arrange(p1, p2, p0, ncol = 2)






NITER <- 1000
curmax <- -Inf

for(i in 1:NITER){
  newd <- PredS[sample(1:M, N),]
  newmax <- logdetbuk(c(logit(newd)), datlist)
  if(newmax > curmax){
    curd <- newd
    curmax <- newmax
    print(c(i, curmax))
  }
}


NITER <- 10000
for(i in 1:NITER){
  newd <- curd
  newd[sample(1:N, 2),] <- PredS[sample(1:M, 2),]
  newmax <- logdetbuk(c(logit(newd)), datlist)
  if(newmax > curmax){
    curd <- newd
    curmax <- newmax
    print(c(i, curmax))
  }
}

testd <- rbind(c(1/11, 1/11), c(1/11, 10/11), c(10/11, 1/11), c(10/11, 10/11),
               c(1/11, 2/11), c(2/11, 1/11),  c(1/11, 10/11), c(10/11, 1/11),
               c(9/11, 10/11), c(10/11, 9/11))
colnames(testd) <- c("u", "v")

psodatNEW <- data.frame(u =curd[1:N], v = curd[1:N + N])
newplot <- qplot(u, v, data = data.frame(PredS), size = I(1), xlim=c(0,1), ylim=c(0,1)) +
  geom_point(aes(u, v), data = psodatNEW, color = I("red"), size = I(3))

psodatTEST <- data.frame(testd)
testplot <- qplot(u, v, data = data.frame(PredS), size = I(1), xlim=c(0,1), ylim=c(0,1)) +
  geom_point(aes(u, v), data = psodatTEST, color = I("blue"), size = I(3))


grid.arrange(newplot, testplot, ncol = 2)

Scur <- sigbuk(c(logit(curd)), datlist)
Stest <- sigbuk(c(logit(testd)), datlist)

-log(det(Scur))/2
-log(det(Stest))/2


oldcurd <- curd
oldmax <- logdetbuk(c(logit(oldcurd)), datlist)

curd <- testd
curmax <- logdetbuk(c(logit(curd)), datlist)
NITER <- 10000
for(i in 1:NITER){
  newd <- curd
  newd[sample(1:N, 1),] <- PredS[sample(1:M, 1),]
  newmax <- logdetbuk(c(logit(newd)), datlist)
  if(newmax > curmax){
    curd <- newd
    curmax <- newmax
    print(c(i, curmax, oldmax))
  }
}


psodatNEWNEW <- data.frame(u =curd[1:N], v = curd[1:N + N])
newnewplot <- qplot(u, v, data = data.frame(PredS), size = I(1), xlim=c(0,1), ylim=c(0,1)) +
  geom_point(aes(u, v), data = psodatNEWNEW, color = I("red"), size = I(3))

grid.arrange(newplot, testplot, newnewplot, ncol = 2)
