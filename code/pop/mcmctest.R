source("popfun.R")
source("mcmcfun.R")
source("psofun.R")
load("popdat/popdat.RData")
library(rstan)
library(MCMCpack)




nswarms <- c(100, 1000)
ndeltas <- c(50)
psoiter <- 10000
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496
mcmciter <- 50000
df <- 100

out <- NULL
psoarglist <- list()
bbpso1arglist <- list()
bbpso2arglist <- list()

for(nswarm in nswarms){
  i <- which(nswarm == nswarms)
  psoarglist[[i]] <- list()
  bbpso1arglist[[i]] <- list()
  bbpso2arglist[[i]] <- list()
  for(ndelta in ndeltas){
    j <- which(ndelta == ndeltas)
    print(c(i,j))
    nbhd <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)})
    poispsoiidinit <- matrix(runif((ndelta + nbeta + 1)*nswarm, -1, 1), ncol = nswarm)
    datlistiid <- list(dat = popdat, nbeta = nbeta, ndelta = ndelta,
                    betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                    aphi = 1, bphi = 1)
    psopoisiid <- pso(psoiter, nswarm, inertia, cognitive, social, poispsoiidinit, nbhd,
                      poislpostiid, datlistiid)
    poisiidmu <- psopoisiid$argmax
    poisiidinit <- psopoisiid$argmax
    poisiidlpbest <- psopoisiid$max
    psopoisiid2 <- bbpso(psoiter, nswarm, 1, 0.1, poispsoiidinit, nbhd,
                   poislpostiid, 1, TRUE, c(0.5, 0.5), datlistiid)
    poisiidmu2 <- psopoisiid2$argmax
    poisiidinit2 <- psopoisiid2$argmax
    poisiidlpbest2 <- psopoisiid2$max
    psopoisiid5 <- bbpso(psoiter, nswarm, 1, 0.3, poispsoiidinit, nbhd,
                         poislpostiid, 5, TRUE, c(0.5, 0.5), datlistiid)
    poisiidmu5 <- psopoisiid5$argmax
    poisiidinit5 <- psopoisiid5$argmax
    poisiidlpbest5 <- psopoisiid5$max
    poisiidind <- indmetrop(mcmciter, poislpostiid, poislposthessiid,
                         poisiidmu, poisiidmu, df,
                         poisiidlpbest, datlistiid)
    poisiidind2 <- indmetrop(mcmciter, poislpostiid, poislposthessiid,
                             poisiidmu2, poisiidmu2, df,
                             poisiidlpbest2, datlistiid)
    poisiidind5 <- indmetrop(mcmciter, poislpostiid, poislposthessiid,
                             poisiidmu5, poisiidmu5, df,
                             poisiidlpbest5, datlistiid)
    out <- rbind(out, data.frame(ndelta = ndelta, pso = "pso", nswarm = nswarm,
                                 max = poisiidlpbest,
                                 acc = mean(poisiidind$acc[-c(1:10000)])))
    out <- rbind(out, data.frame(ndelta = ndelta, pso = "bbpso1", nswarm = nswarm,
                                 max = poisiidlpbest2,
                                 acc = mean(poisiidind2$acc[-c(1:10000)])))
    out <- rbind(out, data.frame(ndelta = ndelta, pso = "bbpso2", nswarm = nswarm,
                                 max = poisiidlpbest5,
                                 acc = mean(poisiidind5$acc[-c(1:10000)])))
    psoarglist[[i]][[j]] <- poisiidmu
    bbpso1arglist[[i]][[j]] <- poisiidmu2
    bbpso2arglist[[i]][[j]] <- poisiidmu5
  }
}


out


   ndelta    pso nswarm        max      acc
1      10    pso    100 3313269857 0.887800
2      10 bbpso1    100 3313269857 0.885325
3      10 bbpso2    100 3313269857 0.884325
4      30    pso    100 3345208711 0.133775
5      30 bbpso1    100 3345208651 0.004300
6      30 bbpso2    100 3345208048 0.010000
7      10    pso   1000 3313269857 0.883575
8      10 bbpso1   1000 3313269857 0.888575
9      10 bbpso2   1000 3313269857 0.889375
10     30    pso   1000 3345208712 0.383975
11     30 bbpso1   1000 3345208713 0.860650
12     30 bbpso2   1000 3345208713 0.865650


par(mfrow=c(3,1))
plot(ts(psopoisiid$maxes[-c(1:100)]))
plot(ts(psopoisiid2$maxes[-c(1:100)]))
plot(ts(psopoisiid5$maxes[-c(1:100)]))

k1 <- 1
k2 <- 1500
k3 <- 2000
plot(ts(psopoisiid$maxes[k1:k3]), ylim = c(min(c(psopoisiid$maxes[k2:k3], psopoisiid2$maxes[k2:k3], psopoisiid5$maxes[k2:k3])), max(c(psopoisiid$maxes[k2:k3], psopoisiid2$maxes[k2:k3], psopoisiid5$maxes[k2:k3]))))
lines(k1:k3, psopoisiid2$maxes[k1:k3], col = "red")
lines(k1:k3, psopoisiid5$maxes[k1:k3], col = "blue")

nbeta <- 1
ndeltasiid <- c(10, 20, 30)
ndeltasfull <- c(5, 7, 9)

psoiter <- 1000
nswarm <- 100
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496

nbhd <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)})

nsubs <- floor(nrow(popdat)/c(4, 2))

mcmciter <- 10000
dfs <- c(1, 4, 10, 100)


ndelta <- 30

poispsoiidinit <- matrix(runif((ndelta + nbeta + 1)*nswarm, -1, 1), ncol = nswarm)
lnormpsoiidinit <- matrix(runif((ndelta + nbeta + 1 + 1)*nswarm, -1, 1), ncol = nswarm)
datlistiid <- list(dat = popdat, nbeta = nbeta, ndelta = ndelta,
                    betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                    aphi = 1, bphi = 1)

psopoisiid <- bbpso(psoiter, nswarm, 1, 0.3, poispsoiidinit, nbhd,
                    poislpostiid, 3, TRUE, c(0.5, 0.5), datlistiid, ccc = 0.1)
poisiidmu <- psopoisiid$argmax
poisiidlpbest <- psopoisiid$max

psopoisiid2 <- pso(psoiter, nswarm, inertia, cognitive, social, poispsoiidinit,
                   nbhd, poislpostiid, datlistiid)
poisiidmu2 <- psopoisiid2$argmax
poisiidlpbest2 <- psopoisiid2$max


mcmciter <- 50000
df <- 100

poisiidind <- indmetrop(mcmciter, poislpostiid, poislposthessiid,
                         poisiidmu, poisiidmu, df,
                         poisiidlpbest, datlistiid)

poisiidind2 <- indmetrop(mcmciter, poislpostiid, poislposthessiid,
                         poisiidmu2, poisiidmu2, df,
                         poisiidlpbest2, datlistiid)



plot(ts(poisiidind$lpbests[-c(1:10000)]), ylim=c(poisiidind$lpbests[10000+1], poisiidind2$lpbests[1]))
lines(1:(length(poisiidind2$lpbests)-10000), poisiidind2$lpbests[-c(1:10000)], col = "red")

mean(poisiidind$acc[-c(1:10000)])
mean(poisiidind2$acc[-c(1:10000)])

mean(poisiidind$acc[c(1:10000)])
mean(poisiidind2$acc[c(1:10000)])


c(max(poisiidind$lpbests), poisiidlpbest)
c(max(poisiidind2$lpbests), poisiidlpbest2)

k3 <- 300
k4 <- psoiter
plot(ts(psopoisiid6$maxes[k3:k4]))
lines(1:(k4 + 1 - k3), psopoisiid$maxes[k3:k4], col="red")
lines(1:(k4 + 1 - k3), psopoisiid2$maxes[k3:k4], col="blue")
lines(1:(k4 + 1 - k3), psopoisiid3$maxes[k3:k4], col="green")
lines(1:(k4 + 1 - k3), psopoisiid4$maxes[k3:k4], col="orange")
lines(1:(k4 + 1 - k3), psopoisiid5$maxes[k3:k4], col="brown")



psopoisiid3 <- bbpso(psoiter, nswarm, 1, 0.3, poispsoiidinit, nbhd,
                   poislpostiid, 1, TRUE, c(0.5, 0.5), datlistiid)
poisiidmu3 <- psopoisiid3$argmax
poisiidinit3 <- psopoisiid3$argmax
poisiidlpbest3 <- psopoisiid3$max

psopoisiid4 <- bbpso(psoiter, nswarm, 1, 0.1, poispsoiidinit, nbhd,
                   poislpostiid, 5, TRUE, c(0.5, 0.5), datlistiid)
poisiidmu4 <- psopoisiid4$argmax
poisiidinit4 <- psopoisiid4$argmax
poisiidlpbest4 <- psopoisiid4$max

psopoisiid5 <- bbpso(psoiter, nswarm, 1, 0.3, poispsoiidinit, nbhd,
                   poislpostiid, 5, TRUE, c(0.5, 0.5), datlistiid)
poisiidmu5 <- psopoisiid4$argmax
poisiidinit5 <- psopoisiid4$argmax
poisiidlpbest5 <- psopoisiid4$max


c(poisiidlpbest, poisiidlpbest2, poisiidlpbest3, poisiidlpbest4, poisiidlpbest5)

c(poisiidlpbest, poisiidlpbest2, poisiidlpbest3, poisiidlpbest4, poisiidlpbest5)) - 3345208713


k <- 900
plot(ts(psopoisiid$maxes[-c(1:k)]), ylim = c(min(c(psopoisiid4$maxes[-c(1:k)],
                                                   psopoisiid$maxes[-c(1:k)],
                                                   psopoisiid2$maxes[-c(1:k)],
                                                    psopoisiid3$maxes[-c(1:k)])),
                                             max(c(psopoisiid$maxes[-c(1:k)],
                                                   psopoisiid4$maxes[-c(1:k)],
                                                   psopoisiid2$maxes[-c(1:k)],
                                                    psopoisiid3$maxes[-c(1:k)]))))
lines(1:length(psopoisiid$maxes[-c(1:k)]), psopoisiid2$maxes[-c(1:k)], col = "red")
lines(1:length(psopoisiid$maxes[-c(1:k)]), psopoisiid3$maxes[-c(1:k)], col = "blue")
lines(1:length(psopoisiid$maxes[-c(1:k)]), psopoisiid4$maxes[-c(1:k)], col = "green")

k <- 1000
k2 <- 4000
plot(ts(psopoisiid$maxes[-c(1:k)]), ylim = c(min(c(psopoisiid4$maxes[-c(1:k2)],
                                                   psopoisiid$maxes[-c(1:k2)],
                                                   psopoisiid2$maxes[-c(1:k2)],
                                                   psopoisiid3$maxes[-c(1:k2)],
                                                   psopoisiid5$maxes[-c(1:k2)])),
                                             max(c(psopoisiid$maxes[-c(1:k2)],
                                                   psopoisiid4$maxes[-c(1:k2)],
                                                   psopoisiid$maxes[-c(1:k2)],
                                                   psopoisiid3$maxes[-c(1:k2)]))))
lines(1:length(psopoisiid$maxes[-c(1:k)]), psopoisiid2$maxes[-c(1:k)], col = "red")
lines(1:length(psopoisiid$maxes[-c(1:k)]), psopoisiid3$maxes[-c(1:k)], col = "blue")
lines(1:length(psopoisiid$maxes[-c(1:k)]), psopoisiid4$maxes[-c(1:k)], col = "green")
lines(1:length(psopoisiid$maxes[-c(1:k)]), psopoisiid5$maxes[-c(1:k)], col = "orange")


mcmciter <- 10000

poisiidind <- indmetrop(mcmciter, poislpostiid, poislposthessiid,
                         poisiidmu, poisiidmu, df,
                         poisiidlpbest, datlistiid)

poisiidind2 <- indmetrop(mcmciter, poislpostiid, poislposthessiid,
                         poisiidmu2, poisiidmu2, df,
                         poisiidlpbest2, datlistiid)

poisiidind3 <- indmetrop(mcmciter, poislpostiid, poislposthessiid,
                         poisiidmu3, poisiidmu3, df,
                         poisiidlpbest3, datlistiid)


c(mean(poisiidind$acc[-c(1:1000)]),
  mean(poisiidind2$acc[-c(1:1000)]))


poisiidindgibbsout <- poisiidindgibbs(mcmciter, poisiidmu, poisiidmu, df,
                                          poisiidlpbest, datlistiid)


poisiidrwgibbsout <- poisiidrwgibbs(mcmciter,
                                        rep(0, nbeta + ndelta + 1),
                                        datlistiid)

sighatiid <- cov(poisiidrwgibbsout$draws[-c(1:10000),1:(nbeta + ndelta)])
cholhatiid <- t(chol(sighatiid))
blockinit <- poisiidrwgibbsout$draws[mcmciter + 1,]
poisiidblockrwgibbsout <- poisiidblockrwgibbs(mcmciter, blockinit, datlistiid, cholhatiid)

c(mean(poisiidind$acc),
  mean(poisiidindgibbsout$acc),
  mean(poisiidrwgibbsout$acc),
  mean(poisiidblockrwgibbsout$acc))


par(mfrow=c(4,4))
for(i in 1:4){
  plot(ts((poisiidind$draws[-c(1:50000),i])))
  plot(ts(poisiidindgibbsout$draws[-c(1:50000),i]))
  plot(ts(poisiidrwgibbsout$draws[-c(1:50000),i]))
  plot(ts(poisiidblockrwgibbsout$draws[-c(1:50000),i]))
}

cbind(summary(mcmc((poisiidind$draws[-c(1:10000),])))[[1]][,1],
      summary(mcmc(poisiidindgibbsout$draws[-c(1:10000),]))[[1]][,1],
      summary(mcmc(poisiidrwgibbsout$draws[-c(1:10000),]))[[1]][,1],
      summary(mcmc(poisiidblockrwgibbsout$draws[-c(1:10000),]))[[1]][,1])



cbind(summary(mcmc(exp(poisiidind$draws[-c(1:10000),])))[[1]][4,1],
      summary(mcmc(poisiidindgibbsout$draws[-c(1:10000),]))[[1]][4,1],
      summary(mcmc(poisiidrwgibbsout$draws[-c(1:10000),]))[[1]][4,1],
      summary(mcmc(poisiidblockrwgibbsout$draws[-c(1:10000),]))[[1]][4,1])






nell <- ndelta*(ndelta + 1)/2
poispsofullinit <- matrix(runif((ndelta + nbeta + nell)*nswarm, -1, 1), ncol = nswarm)
lnormpsofullinit <- matrix(runif((ndelta + nbeta + nell + 1)*nswarm, -1, 1), ncol = nswarm)
ominvscale <- diag(ndelta)
K <- commutation.matrix(ndelta)
M <- elimination.matrix(ndelta)
N2 <- 2*N.matrix(ndelta)
R <- matrix(0, nell, nell)
MKprime <- tcrossprod(M, K)
diag(R)[(ndelta+1)*(1:ndelta) - (1:ndelta)*(1 + 1:ndelta)/2 + 1:ndelta - ndelta] <-
  (ndelta + 2 - 1:ndelta)
dvelldvell <- - M%*%tcrossprod(kronecker(diag(ndelta), ominvscale), M)
datlistfull <- list(dat = popdat, nbeta = nbeta, ndelta = ndelta,
                    betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                    aphi = 1, bphi = 1, omdf = ndelta + 1, ominvscale = diag(ndelta),
                    K = K, M = M, N2 = N2, R = R, MKprime = MKprime, dvelldvell = dvelldvell)

psopoisfull <- pso(psoiter, nswarm, inertia, cognitive, social, poispsofullinit, nbhd,
                   poislpostfull, datlistfull)
poisfullmu <- psopoisfull$argmax
poisfullinit <- psopoisfull$argmax
poisfulllpbest <- psopoisfull$max
poisfulllpbest

poisfullindgibbsout <- poisfullindgibbs(mcmciter, poisfullmu, poisfullmu, df,
                                        poisfulllpbest, datlistfull)
mean(poisfullindgibbsout$acc[-c(1:1000)])


psolnormfull <- pso(psoiter, nswarm, inertia, cognitive, social, lnormpsofullinit, nbhd,
                    lnormlpostfull, datlistfull)
lnormfullmu <- psolnormfull$argmax
lnormfullinit <- psolnormfull$argmax
lnormfulllpbest <- psolnormfull$max
lnormfulllpbest


psopoisfull2 <- bbpso(psoiter, nswarm, 1, 0.3, poispsofullinit, nbhd,
                      poislpostfull, 5, TRUE, c(0.5, 0.5), datlistfull)
poisfullmu2 <- psopoisfull2$argmax
poisfullinit2 <- psopoisfull2$argmax
poisfulllpbest2 <- psopoisfull2$max
poisfulllpbest2

psolnormfull2 <- bbpso(psoiter, nswarm, 1, 0.3, lnormpsofullinit, nbhd,
                   lnormlpostfull, 5, TRUE, c(0.5, 0.5), datlistfull)
lnormfullmu2 <- psolnormfull2$argmax
lnormfullinit2 <- psolnormfull2$argmax
lnormfulllpbest2 <- psolnormfull2$max
lnormfulllpbest2

df <- 100
mcmciter <- 10000

poisfullind <- indmetrop(mcmciter, poislpostfull, poislposthessfull,
                         poisfullmu, poisfullmu, df,
                         poisfulllpbest, datlistfull)

mean(poisfullind$acc [-c(1:1000)])


lnormfullind <- indmetrop(mcmciter, lnormlpostfull, lnormlposthessfull,
                          lnormfullmu, lnormfullmu, df,
                          lnormfulllpbest, datlistfull)

poisfullind2 <- indmetrop(mcmciter, poislpostfull, poislposthessfull,
                         poisfullmu2, poisfullmu2, df,
                         poisfulllpbest2, datlistfull)

lnormfullind2 <- indmetrop(mcmciter, lnormlpostfull, lnormlposthessfull,
                          lnormfullmu2, lnormfullmu2, df,
                          lnormfulllpbest2, datlistfull)


c(mean(poisfullind$acc [-c(1:1000)]), mean(poisfullind2$acc [-c(1:1000)]))
c(mean(lnormfullind$acc[-c(1:1000)]), mean(lnormfullind2$acc[-c(1:1000)]))


lnormfullindgibbsout <- lnormfullindgibbs(mcmciter, lnormfullmu, lnormfullmu, df,
                                          lnormfulllpbest, datlistfull)

lnormfullindgibbsout2 <- lnormfullindgibbs(mcmciter, lnormfullmu2, lnormfullmu2, df,
                                           lnormfulllpbest2, datlistfull)

poisfullindgibbsout <- poisfullindgibbs(mcmciter, poisfullmu, poisfullmu, df,
                                        poisfulllpbest, datlistfull)

poisfullindgibbsout2 <- poisfullindgibbs(mcmciter, poisfullmu2, poisfullmu2, df,
                                         poisfulllpbest2, datlistfull)

c(mean(poisfullindgibbsout$acc [-c(1:1000)]), mean(poisfullindgibbsout2$acc [-c(1:1000)]))
c(mean(lnormfullindgibbsout$acc[-c(1:1000)]), mean(lnormfullindgibbsout2$acc[-c(1:1000)]))



lnormfullrwgibbsout <- lnormfullrwgibbs(mcmciter,
                                        rep(0, nbeta + ndelta + ndelta*(ndelta + 1)/2 + 1),
                                        datlistfull)

sighatfull <- cov(lnormfullrwgibbsout$draws[-c(1:10000),1:(nbeta + ndelta)])
cholhatfull <- t(chol(sighatfull))
blockinit <- lnormfullrwgibbsout$draws[mcmciter + 1,]
lnormfullblockrwgibbsout <- lnormfullblockrwgibbs(mcmciter, blockinit, datlistfull, cholhatfull)

mcmciter <- 10000

lnormfullgibbsout <- lnormfullgibbs(mcmciter, rep(0, nbeta + ndelta + ndelta*(ndelta + 1)/2 + 1),
                                    datlistfull)

c(mean(lnormfullind$acc),
  mean(lnormfullindgibbsout$acc),
  mean(lnormfullrwgibbsout$acc),
  mean(lnormfullblockrwgibbsout$acc))


par(mfrow=c(4,5))
for(i in 1:4){
  plot(ts((lnormfullind$draws[-c(1:50000),i])))
  plot(ts(lnormfullindgibbsout$draws[-c(1:50000),i]))
  plot(ts(lnormfullrwgibbsout$draws[-c(1:50000),i]))
  plot(ts(lnormfullblockrwgibbsout$draws[-c(1:50000),i]))
  plot(ts(lnormfullgibbsout$draws[-c(1:1000),i]))
}

cbind(summary(mcmc((lnormfullind$draws[-c(1:10000),])))[[1]][,1],
      summary(mcmc(lnormfullindgibbsout$draws[-c(1:10000),]))[[1]][,1],
      summary(mcmc(lnormfullrwgibbsout$draws[-c(1:10000),]))[[1]][,1],
      summary(mcmc(lnormfullblockrwgibbsout$draws[-c(1:10000),]))[[1]][,1],
      summary(mcmc(lnormfullgibbsout$draws[-c(1:1000),]))[[1]][,1])


cbind(summary(mcmc(exp(lnormfullind$draws[-c(1:10000),])))[[1]][7,1],
      summary(mcmc(lnormfullindgibbsout$draws[-c(1:10000),]))[[1]][7,1],
      summary(mcmc(lnormfullrwgibbsout$draws[-c(1:10000),]))[[1]][7,1],
      summary(mcmc(lnormfullblockrwgibbsout$draws[-c(1:10000),]))[[1]][7,1],
      summary(mcmc(lnormfullgibbsout$draws[-c(1:1000),]))[[1]][7,1])




mcmciter <- 100000
df <- 1

poisfullind <- indmetrop(mcmciter, poislpostfull, poislposthessfull,
                          poisfullmu, poisfullmu, df,
                          poisfulllpbest, datlistfull)

poisfullindgibbsout <- poisfullindgibbs(mcmciter, poisfullmu, poisfullmu, df,
                                          poisfulllpbest, datlistfull)


poisfullrwgibbsout <- poisfullrwgibbs(mcmciter,
                                        rep(0, nbeta + ndelta + ndelta*(ndelta + 1)/2),
                                        datlistfull)

sighatfull <- cov(poisfullrwgibbsout$draws[-c(1:10000),1:(nbeta + ndelta)])
cholhatfull <- t(chol(sighatfull))
blockinit <- poisfullrwgibbsout$draws[mcmciter + 1,]
poisfullblockrwgibbsout <- poisfullblockrwgibbs(mcmciter, blockinit, datlistfull, cholhatfull)

c(mean(poisfullind$acc),
  mean(poisfullindgibbsout$acc),
  mean(poisfullrwgibbsout$acc),
  mean(poisfullblockrwgibbsout$acc))

par(mfrow=c(4,4))
for(i in 1:4){
  plot(ts((poisfullind$draws[-c(1:50000),i+2])))
  plot(ts(poisfullindgibbsout$draws[-c(1:50000),i+2]))
  plot(ts(poisfullrwgibbsout$draws[-c(1:50000),i+2]))
  plot(ts(poisfullblockrwgibbsout$draws[-c(1:50000),i+2]))
}

cbind(summary(mcmc((poisfullind$draws[-c(1:10000),])))[[1]][,1],
      summary(mcmc(poisfullindgibbsout$draws[-c(1:10000),]))[[1]][,1],
      summary(mcmc(poisfullrwgibbsout$draws[-c(1:10000),]))[[1]][,1],
      summary(mcmc(poisfullblockrwgibbsout$draws[-c(1:10000),]))[[1]][,1])


cbind(summary(mcmc(abs(poisfullind$draws[-c(1:10000),])))[[1]][,1],
      summary(mcmc(poisfullindgibbsout$draws[-c(1:10000),]))[[1]][,1],
      summary(mcmc(poisfullrwgibbsout$draws[-c(1:10000),]))[[1]][,1],
      summary(mcmc(poisfullblockrwgibbsout$draws[-c(1:10000),]))[[1]][,1])





poispsoiidinit <- matrix(runif((ndelta + nbeta + 1)*nswarm, -1, 1), ncol = nswarm)
lnormpsoiidinit <- matrix(runif((ndelta + nbeta + 1 + 1)*nswarm, -1, 1), ncol = nswarm)
datlistiid <- list(dat = popdat, nbeta = nbeta, ndelta = ndelta,
                    betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                    aphi = 1, bphi = 1)

psopoisiid <- pso(psoiter, nswarm, inertia, cognitive, social, poispsoiidinit, nbhd,
                   poislpostiid, datlistiid)
poisiidmu <- psopoisiid$argmax
poisiidinit <- psopoisiid$argmax
poisiidlpbest <- psopoisiid$max

psolnormiid <- pso(psoiter, nswarm, inertia, cognitive, social, lnormpsoiidinit, nbhd,
                    lnormlpostiid, datlistiid)
lnormiidmu <- psolnormiid$argmax
lnormiidinit <- psolnormiid$argmax
lnormiidlpbest <- psolnormiid$max

mcmciter <- 100000

lnormiidind <- indmetrop(mcmciter, lnormlpostiid, lnormlposthessiid,
                          lnormiidmu, lnormiidmu, df,
                          lnormiidlpbest, datlistiid)

lnormiidindgibbsout <- lnormiidindgibbs(mcmciter, lnormiidmu, lnormiidmu, df,
                                          lnormiidlpbest, datlistiid)


lnormiidrwgibbsout <- lnormiidrwgibbs(mcmciter,
                                        rep(0, nbeta + ndelta + 1 + 1),
                                        datlistiid)

sighatiid <- cov(lnormiidrwgibbsout$draws[-c(1:10000),1:(nbeta + ndelta)])
cholhatiid <- t(chol(sighatiid))
blockinit <- lnormiidrwgibbsout$draws[mcmciter + 1,]
lnormiidblockrwgibbsout <- lnormiidblockrwgibbs(mcmciter, blockinit, datlistiid, cholhatiid)

mcmciter <- 10000

lnormiidgibbsout <- lnormiidgibbs(mcmciter, rep(0, nbeta + ndelta + 1 + 1),
                                    datlistiid)

c(mean(lnormiidind$acc),
  mean(lnormiidindgibbsout$acc),
  mean(lnormiidrwgibbsout$acc),
  mean(lnormiidblockrwgibbsout$acc))


par(mfrow=c(4,5))
for(i in 1:4){
  plot(ts((lnormiidind$draws[-c(1:50000),i])))
  plot(ts(lnormiidindgibbsout$draws[-c(1:50000),i]))
  plot(ts(lnormiidrwgibbsout$draws[-c(1:50000),i]))
  plot(ts(lnormiidblockrwgibbsout$draws[-c(1:50000),i]))
  plot(ts(lnormiidgibbsout$draws[-c(1:1000),i]))
}

cbind(summary(mcmc((lnormiidind$draws[-c(1:10000),])))[[1]][,1],
      summary(mcmc(lnormiidindgibbsout$draws[-c(1:10000),]))[[1]][,1],
      summary(mcmc(lnormiidrwgibbsout$draws[-c(1:10000),]))[[1]][,1],
      summary(mcmc(lnormiidblockrwgibbsout$draws[-c(1:10000),]))[[1]][,1],
      summary(mcmc(lnormiidgibbsout$draws[-c(1:1000),]))[[1]][,1])


cbind(summary(mcmc(exp(lnormiidind$draws[-c(1:10000),])))[[1]][4,1],
      summary(mcmc(lnormiidindgibbsout$draws[-c(1:10000),]))[[1]][4,1],
      summary(mcmc(lnormiidrwgibbsout$draws[-c(1:10000),]))[[1]][4,1],
      summary(mcmc(lnormiidblockrwgibbsout$draws[-c(1:10000),]))[[1]][4,1],
      summary(mcmc(lnormiidgibbsout$draws[-c(1:1000),]))[[1]][4,1])

cbind(summary(mcmc(exp(lnormiidind$draws[-c(1:10000),])))[[1]][5,1],
      summary(mcmc(lnormiidindgibbsout$draws[-c(1:10000),]))[[1]][5,1],
      summary(mcmc(lnormiidrwgibbsout$draws[-c(1:10000),]))[[1]][5,1],
      summary(mcmc(lnormiidblockrwgibbsout$draws[-c(1:10000),]))[[1]][5,1],
      summary(mcmc(lnormiidgibbsout$draws[-c(1:1000),]))[[1]][5,1])





mcmciter <- 100000
nbeta <- 1
iidinit <- c(rep(0, ndelta + nbeta), exp(0))
fullinit <- c(rep(0, ndelta + nbeta), rep(exp(0), ndelta*(ndelta + 1)/2))

fullstandat <- list(nobs = length(z), ndelta = ndelta, nbeta = ncol(xmat), zobs = z,
                    lzobs = log(z), S = smat[,1:ndelta], X = xmat,
                    a_phi = aphi, b_phi = bphi, a_sig = asig, b_sig = bsig,
                    mu_beta = mubeta, sig_beta = sigbeta,
                    d_omega = d, Einv_omega = Einv)

### to compare against
poisfulltest <- stan(file = "poppoisfull.stan", data = fullstandat, chains = 1, iter = 1)
poisfullfit <-  stan(fit = poisfulltest,        data = fullstandat, chains = 1, iter = 50000)
poisfullfit <-  stan(fit = poisfullfit,         data = fullstandat, chains = 1, iter = 100000)

poisiidtest <-  stan(file = "poppoisiid.stan", data = fullstandat, chains = 1, iter = 1)
poisiidfit <-   stan(fit = poisiidtest,        data = fullstandat, chains = 1, iter = 10000)
poisiidfit <-   stan(fit = poisiidfit,         data = fullstandat, chains = 1, iter = 50000)












poisfullrwgibbsout <- poisfullrwgibbs(mcmciter,
                                        rep(0, nbeta + ndelta + ndelta*(ndelta + 1)/2),
                                        datlistfull)

sighatfull <- cov(poisfullrwgibbsout$draws[-c(1:10000),1:(nbeta + ndelta)])
cholhatfull <- t(chol(sighatfull))
blockinit <- poisfullrwgibbsout$draws[mcmciter + 1,]
poisfullblockrwgibbsout <- poisfullblockrwgibbs(mcmciter, blockinit, datlistfull, cholhatfull)


cbind(summary(mcmc(poisfullrwgibbsout$draws[-c(1:10000),]))[[1]][,1],
      summary(mcmc(poisfullblockrwgibbsout$draws[-c(1:10000),]))[[1]][,1])






ndelta <- 30
mcmciter <- 10000

poispsoiidinit <- matrix(runif((ndelta + nbeta + 1)*nswarm, -1, 1), ncol = nswarm)
lnormpsoiidinit <- matrix(runif((ndelta + nbeta + 1 + 1)*nswarm, -1, 1), ncol = nswarm)
datlistiid <- list(dat = popdat, nbeta = nbeta, ndelta = ndelta,
                   betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                   aphi = 1, bphi = 1)
psopoisiid <- pso(psoiter, nswarm, inertia, cognitive, social, poispsoiidinit, nbhd,
                  poislpostiid, datlistiid)
poisiidmu <- psopoisiid$argmax
poisiidinit <- psopoisiid$argmax
poisiidlpbest <- psopoisiid$max
psolnormiid <- pso(psoiter, nswarm, inertia, cognitive, social, lnormpsoiidinit, nbhd,
                   lnormlpostiid, datlistiid)
lnormiidmu <- psolnormiid$argmax
lnormiidinit <- psolnormiid$argmax
lnormiidlpbest <- psolnormiid$max

poisiidind <- indmetrop(mcmciter, poislpostiid, poislposthessiid,
                        poisiidmu, poisiidmu, df,
                        poisiidlpbest, datlistiid)

lnormiidind <- indmetrop(mcmciter, lnormlpostiid, lnormlposthessiid,
                         lnormiidmu, lnormiidmu, df,
                         lnormiidlpbest, datlistiid)

lnormiidindgibbsout <- lnormiidindgibbs(mcmciter, lnormiidinit, lnormiidmu, df, lnormiidlpbest, datlistiid)



lnormiidrwgibbsout <- lnormiidrwgibbs(mcmciter*2, rep(0, nbeta + ndelta + 2), datlistiid) 

sighatiid <- cov(lnormiidrwgibbsout$draws[-c(1:2000),1:(nbeta + ndelta)])
cholhatiid <- t(chol(sighatiid))
blockinit <- lnormiidrwgibbsout$draws[mcmciter + 1,]
lnormiidblockrwgibbsout <- lnormiidblockrwgibbs(mcmciter*5, blockinit, datlistiid, cholhatiid, logrwsd = - 6)


mean(lnormiidrwgibbsout$acc)

mean(lnormiidblockrwgibbsout$acc)

par(mfrow=c(4,4))
for(i in 1:4){
  plot(ts(lnormiidrwgibbsout$draws[-c(1:2000),i]))
  plot(ts(lnormiidblockrwgibbsout$draws[-c(1:5000),i]))
  plot(ts(lnormiidind$draws[,i]))
  plot(ts(lnormiidindgibbsout$draws[,i]))
}




cbind(summary(mcmc(lnormiidrwgibbsout$draws[-c(1:2000),]))[[1]][,1],
      summary(mcmc(lnormiidblockrwgibbsout$draws[-c(1:5000),]))[[1]][,1],
      summary(mcmc((lnormiidind$draws)))[[1]][,1],
      summary(mcmc(lnormiidindgibbsout$draws))[[1]][,1])



mubd <- lnormiidmu[1:(nbeta + ndelta)]
sig2 <- exp(lnormiidmu[nbeta + ndelta + 1])
phi2 <- exp(lnormiidmu[nbeta + ndelta + 1 + 1])
mulsig2 <- lnormiidmu[nbeta + ndelta + 1]
mulphi2 <- lnormiidmu[nbeta + ndelta + 2]

hess <- lnormlposthessiid(lnormiidmu, datlistiid)
SigmaBig <- chol2inv(chol(-hess))
Sigma1 <- SigmaBig[1:(nbeta + ndelta), 1:(nbeta + ndelta)]
Sigma12 <- SigmaBig[1:(nbeta + ndelta), nbeta + ndelta + 1:2]
Sigma22 <- SigmaBig[nbeta + ndelta + 1:2, nbeta + ndelta + 1:2]

mn <- mubd + Sigma12%*%chol2inv(chol(Sigma22))%*%(c(log(sig2), log(phi2)) - c(mulsig2, mulphi2))





z <- datlistiid$dat[, 1]
lz <- log(z)
ndat <- length(z)
xmat <- datlistiid$dat[, 1 + 1:nbeta]
smat <- datlistiid$dat[, 1 + nbeta + 1:ndelta]


hess <- lnormlposthessiid(lnormiidmu, datlistiid)
SigmaBig <- chol2inv(chol(-hess))
Sigma1 <- SigmaBig[1:(nbeta + ndelta), 1:(nbeta + ndelta)]
Sigma12 <- SigmaBig[1:(nbeta + ndelta), nbeta + ndelta + 1:2]
Sigma22 <- SigmaBig[nbeta + ndelta + 1:2, nbeta + ndelta + 1:2]

Sigma <- Sigma1 - Sigma12%*%tcrossprod(chol2inv(chol(Sigma22)),Sigma12)
Rsigma <- chol(Sigma)


mean(poisiidind$acc)
mean(lnormiidind$acc)
mean(lnormiidindgibbsout$acc)

par(mfrow=c(2,2))
for(i in 1:2){
  plot(ts(exp(lnormiidind$draws[,ndelta + nbeta + i])))
  plot(ts(lnormiidindgibbsout$draws[,ndelta + nbeta + i]))
}

cbind(summary(mcmc(exp(lnormiidind$draws[,nbeta + ndelta + 1:2])))[[1]][,1],
      summary(mcmc(lnormiidindgibbsout$draws[,nbeta + ndelta + 1:2]))[[1]][,1])

cbind(summary(mcmc((lnormiidind$draws[,nbeta + ndelta + 1:2])))[[1]][,1],
      summary(mcmc(log(lnormiidindgibbsout$draws[,nbeta + ndelta + 1:2])))[[1]][,1])





library(microbenchmark)

xmat <- crossprod(matrix(rnorm(100*3), ncol = 3))

microbenchmark(
  test1 <- pd.solve(xmat),
  test2 <- chol2inv(chol(xmat)))


sigma <- crossprod(matrix(rnorm(27), ncol = 3))
omega <- pd.solve(sigma)
cholsig <- chol(sigma)

rmtcholtest <- rmtchol(10000, rep(0, 3), cholsig, 5)
rmtfixedtest <- rmtfixed(10000, rep(0, 3), 5, )

Sigma <- crossprod(matrix(rnorm(30), ncol = 3))
Prec <- chol2inv(chol(Sigma))
Rprec <- chol(Prec)
mu <- rnorm(3)
par <- rnorm(3) + mu
prop <- rnorm(3) + mu
                         
c(dmt(par, mu, Sigma, df, TRUE) - dmt(prop, mu, Sigma, df, TRUE)) - c(dmtcholprec(par, mu, Rprec, df, TRUE) - dmtcholprec(prop, mu, Rprec, df, TRUE))



                         
H <- -crossprod(matrix(rnorm(30), ncol = 3))
Rprec <- chol(-H)
Sigma <- chol2inv(Rprec)
Rsig <- chol(Sigma)
Rprec2 <- chol(chol2inv(Rsig))

Rprec
Rprec2







nswarm <- 100
psoiter <- 1000
ndelta <- 10
nbeta <- 1
nbhd <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)})
poispsoiidinit <- matrix(runif((ndelta + nbeta + 1)*nswarm, -1, 1), ncol = nswarm)
datlistiid <- list(dat = popdat, nbeta = nbeta, ndelta = ndelta,
                   betamn = 0, betavar = 100, sig2a = 1, sig2b = 1,
                   aphi = 1, bphi = 1)

psoout <- pso(psoiter, nswarm, inertia, cognitive, social, poispsoiidinit, nbhd,
              poislpostiid, datlistiid)

psotuneout <- pso(psoiter, nswarm, inertia, cognitive, social, poispsoiidinit, nbhd,
                  poislpostiid, datlistiid, tune = TRUE, rate = 0.5, ccc = 0.1)

psoinertout <- psoinertfun(psoiter, nswarm, 1, psoiter*0.4, cognitive, social, poispsoiidinit, nbhd,
                    poislpostiid, datlistiid)

psotuneout2 <- pso(psoiter, nswarm, inertia, cognitive, social, poispsoiidinit, nbhd,
                  poislpostiid, datlistiid, tune = TRUE, rate = 0.5, ccc = 0.01)



k1 <- 150
k2 <- 300
plot(ts(psoout$maxes[k1:k2]))
lines(1:(k2-k1 + 1), psoinertout$maxes[k1:k2], col = "red")
lines(1:(k2-k1 + 1), psotuneout$maxes[k1:k2], col = "blue")

lines(1:(k2-k1 + 1), psotuneout2$maxes[k1:k2], col = "green")

k1 <- 200
k2 <- 300
plot(ts(psotuneout$maxes[k1:k2]), col = "blue")
lines(1:(k2-k1 + 1), psotuneout2$maxes[k1:k2], col = "green")
lines(1:(k2-k1 + 1), psoinertout$maxes[k1:k2], col = "red")

psoiter <- 1000
slope <- psoiter*0.4
speed <- 1
plot(ts(1/(1 + (1:psoiter/slope)^speed)), ylim=c(0,1))
slope <- psoiter*0.1
lines(1:psoiter, 1/(1 + (1:psoiter/slope)^speed), col="red")


inertia <- 1/(1 + (iter/slope)^speed)
