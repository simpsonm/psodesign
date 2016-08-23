library(xtable)
library(plyr)
library(ggplot2)
library(reshape2)

psoout <- read.csv("psoout.csv")

sumlast <- subset(psoout, iteration == max(psoout$iteration))


psosum <- ddply(sumlast, .(model, ranef, ndelta, algorithm, nbhd), summarise,
                mean = mean(maxes), median = median(maxes), sd = sd(maxes),
                min=min(maxes), q10=quantile(maxes, 0.1), q25=quantile(maxes, 0.25),
                q75=quantile(maxes, 0.75), q90=quantile(maxes, 0.9), max=max(maxes))

psosum <- psosum[rev(1:nrow(psosum)),]

algset <- c("PSO", "DI-PSO", "BBPSO-MC", "BBPSOxp-MC", "AT-PSO-0.5-0.1", "AT-BBPSOxp-MC-5-0.5-0.1",
            "AT-BBPSO-MC-5-0.5-0.1")

psosum <- subset(psosum, algorithm %in% algset)

psosum$algorithm <-
  mapvalues(psosum$algorithm,
            c("AT-PSO-0.5-0.1", "AT-BBPSO-MC-5-0.5-0.1", "AT-BBPSOxp-MC-5-0.5-0.1"),
            c("AT-PSO", "AT-BBPSO-MC", "AT-BBPSOxp-MC"))
         


poisiid <- subset(psosum, model == "pois" & ranef == "iid" & ndelta == 30)[,c(4,5,6,8)]
poisiid$mean <- round(poisiid$mean - poisiid$mean[3],2)
poisiid$sd[poisiid$mean < -100000] <- NA
poisiid$mean[poisiid$mean < -100000] <- NA
poisiid <- rbind(subset(poisiid, nbhd == "global")[c(1,4,3,2,5,7,6),],
                 subset(poisiid, nbhd == "ring-3")[c(1,4,3,2,5,7,6),],
                 subset(poisiid, nbhd == "ring-1")[c(1,4,3,2,5,7,6),])
poisiid

poisfull <- subset(psosum, model == "pois" & ranef == "full" & ndelta == 15)[,c(4,5,6,8)]
poisfull$mean <- round(poisfull$mean - poisfull$mean[1],2)
poisfull$sd[poisfull$mean < -100000] <- NA
poisfull$mean[poisfull$mean < -100000] <- NA
poisfull <- rbind(subset(poisfull, nbhd == "global")[c(1,4,3,2,5,7,6),],
                  subset(poisfull, nbhd == "ring-3")[c(1,4,3,2,5,7,6),],
                  subset(poisfull, nbhd == "ring-1")[c(1,4,3,2,5,7,6),])
poisfull

lnormiid <- subset(psosum, model == "lnorm" & ranef == "iid" & ndelta == 30)[,c(4,5,6,8)]
lnormiid$mean <- round(lnormiid$mean - lnormiid$mean[1],2)
lnormiid$sd[lnormiid$mean < -100000] <- NA
lnormiid$mean[lnormiid$mean < -100000] <- NA
lnormiid <- rbind(subset(lnormiid, nbhd == "global")[c(1,4,3,2,5,7,6),],
                  subset(lnormiid, nbhd == "ring-3")[c(1,4,3,2,5,7,6),],
                  subset(lnormiid, nbhd == "ring-1")[c(1,4,3,2,5,7,6),])
lnormiid

lnormfull <- subset(psosum, model == "lnorm" & ranef == "full" & ndelta == 15)[,c(4,5,6,8)]
lnormfull$mean <- round(lnormfull$mean - lnormfull$mean[3],2)
lnormfull$sd[lnormfull$mean < -100000] <- NA
lnormfull$mean[lnormfull$mean < -100000] <- NA
lnormfull <- rbind(subset(lnormfull, nbhd == "global")[c(1,4,3,2,5,7,6),],
                   subset(lnormfull, nbhd == "ring-3")[c(1,4,3,2,5,7,6),],
                   subset(lnormfull, nbhd == "ring-1")[c(1,4,3,2,5,7,6),])
lnormfull



library(xtable)
print(xtable(lnormfull), include.rownames=FALSE)







which.min(poisiid$mean)

subset(psosum, model == "pois" & ranef == "iid" & ndelta == 30 & nbhd == "ring-3")
subset(psosum, model == "pois" & ranef == "iid" & ndelta == 30 & nbhd == "global")

subset(psosum, model == "pois" & ranef == "iid" & ndelta == 30 & nbhd == "global")[14,]
subset(psosum, model == "pois" & ranef == "iid" & ndelta == 30 & nbhd == "ring-3")[10,]

subset(psosum, model == "pois" & ranef == "iid" & ndelta == 30 & nbhd == "ring-3")[10,]

load("accout.RData")

accout$id <- as.factor(paste(accout$model, accout$ranef, accout$ndelta, accout$niter, accout$pso, accout$mcmc, sep=""))
accout2 <- accout[1,]

for(i in 1:length(levels(accout$id))){
  accout2 <- rbind(accout2, subset(accout, id == levels(accout$id)[i])[1,])
}
accout2 <- accout2[-1,]



longmelt <- melt(accout2[,1:8],
                 id.vars = c("model", "ranef", "ndelta", "nswarm", "niter", "pso", "mcmc"))
longmelt$pso <- factor(longmelt$pso, levels = rev(levels(accout2$pso)))
##longmelt$pso <- mapvalues(longmelt$pso, from=c("AT-BBPSOxp-MC-1", "AT-BBPSOxp-MC-2"),
##                          to=c("AT-1-0.1", "AT-5-0.3"))
longmelt$ranef <- factor(longmelt$ranef, levels = c("iid", "full"))

imhcast <- dcast(longmelt, niter + pso ~ ranef + model + ndelta,
                 subset = .(mcmc == "IMH"))

imhcast$niter[-c(1,5,9, 13)] <- ""
print(xtable(imhcast[,c(2,1,3:14)], digits=c(0,0,0,rep(2,12))), include.rownames=FALSE)

imhwgcast <- dcast(longmelt, niter + pso ~ ranef + model + ndelta,
                   subset = .(mcmc == "IMHwG"))
imhwgcast$niter[-c(1,5,9,13)] <- ""
print(xtable(imhwgcast[,c(2,1,3:14)], digits=c(0,0,0,rep(2,12))), include.rownames=FALSE)





print(xtable(dcast(longmelt, niter + ndelta ~ pso,
                   subset = .(model == "pois"  & ranef == "iid" & niter > 0 &
                              mcmc == "IMH" & !(pso %in% c("AT-BBPSO-2", "AT-PSO-2")))),
           , digits=c(0,0,0,2,2,2,2,2)),
      include.rownames=FALSE)

print(xtable(dcast(longmelt, niter + ndelta ~ pso,
                   subset = .(model == "pois"  & ranef == "full" & niter > 0 &
                              mcmc == "IMH" & !(pso %in% c("AT-BBPSO-2", "AT-PSO-2")))),
           , digits=c(0,0,0,2,2,2,2,2)),
      include.rownames=FALSE)

print(xtable(dcast(longmelt, niter + ndelta ~ pso,
                   subset = .(model == "lnorm"  & ranef == "iid" & niter > 0 &
                              mcmc == "IMH" & !(pso %in% c("AT-BBPSO-2", "AT-PSO-2")))),
           , digits=c(0,0,0,2,2,2,2,2)),
      include.rownames=FALSE)

print(xtable(dcast(longmelt, niter + ndelta ~ pso,
                   subset = .(model == "lnorm"  & ranef == "full" & niter > 0 &
                              mcmc == "IMH" & !(pso %in% c("AT-BBPSO-2", "AT-PSO-2")))),
           , digits=c(0,0,0,2,2,2,2,2)),
      include.rownames=FALSE)





print(xtable(dcast(longmelt, niter + ndelta ~ pso,
                   subset = .(model == "pois"  & ranef == "iid" & niter > 0 &
                              mcmc == "IMHwG" & !(pso %in% c("AT-BBPSO-2", "AT-PSO-2")))),
           , digits=c(0,0,0,2,2,2,2,2)),
      include.rownames=FALSE)

print(xtable(dcast(longmelt, niter + ndelta ~ pso,
                   subset = .(model == "pois"  & ranef == "full" & niter > 0 &
                              mcmc == "IMHwG" & !(pso %in% c("AT-BBPSO-2", "AT-PSO-2")))),
           , digits=c(0,0,0,2,2,2,2,2)),
      include.rownames=FALSE)

print(xtable(dcast(longmelt, niter + ndelta ~ pso,
                   subset = .(model == "lnorm"  & ranef == "iid" & niter > 0 &
                              mcmc == "IMHwG" & !(pso %in% c("AT-BBPSO-2", "AT-PSO-2")))),
           , digits=c(0,0,0,2,2,2,2,2)),
      include.rownames=FALSE)

print(xtable(dcast(longmelt, niter + ndelta ~ pso,
                   subset = .(model == "lnorm"  & ranef == "full" & niter > 0 &
                              mcmc == "IMHwG" & !(pso %in% c("AT-BBPSO-2", "AT-PSO-2")))),
           , digits=c(0,0,0,2,2,2,2,2)),
      include.rownames=FALSE)




load("mcmcout.RData")
rownames(mcmcout) <- NULL
head(mcmcout)
mcmcout$nefftime <- mcmcout$itertime/mcmcout$neff*10000
mcmcout$alg <- mapvalues(mcmcout$alg, c("IMHwGibbs", "blockRWwG"), c("IMHwG", "B-RWwG"))
mcmcout$alg <- factor(mcmcout$alg, levels(mcmcout$alg)[c(6,5,4,3,2,1)])
mcmcout$model <- mapvalues(mcmcout$model, c("lnorm", "pois"), c("lognormal", "Poisson"))

mcmcoutmelt <- melt(mcmcout, id.vars = c("model", "ranef", "ndelta", "alg"))

poisiidneff <- dcast(mcmcoutmelt, ndelta ~ alg,
                     subset = .(model == "Poisson" & ranef == "iid" & variable == "neff"))
poisiidnefftime <- dcast(mcmcoutmelt, ndelta ~ alg,
                         subset = .(model == "Poisson" & ranef == "iid" & variable == "nefftime"))
poisfullneff <- dcast(mcmcoutmelt, ndelta ~ alg,
                      subset = .(model == "Poisson" & ranef == "full" & variable == "neff"))
poisfullnefftime <- dcast(mcmcoutmelt, ndelta ~ alg,
                          subset = .(model == "Poisson" & ranef == "full" & variable == "nefftime"))
lnormiidneff <- dcast(mcmcoutmelt, ndelta ~ alg,
                      subset = .(model == "lognormal" & ranef == "iid" & variable == "neff"))
lnormiidnefftime <- dcast(mcmcoutmelt, ndelta ~ alg,
                          subset = .(model == "lognormal" & ranef == "iid" & variable == "nefftime"))
lnormfullneff <- dcast(mcmcoutmelt, ndelta ~ alg,
                       subset = .(model == "lognormal" & ranef == "full" & variable == "neff"))
lnormfullnefftime <- dcast(mcmcoutmelt, ndelta ~ alg,
                           subset = .(model == "lognormal" & ranef == "full" & variable == "nefftime"))

print(xtable(rbind(poisiidneff, poisfullneff), digits = 0), include.rownames = FALSE)

print(xtable(rbind(lnormiidneff, lnormfullneff), digits = 0), include.rownames = FALSE)

print(xtable(rbind(poisiidnefftime, poisfullnefftime), digits = c(0,0,rep(0,5))),
      include.rownames = FALSE)

print(xtable(rbind(lnormiidnefftime, lnormfullnefftime), digits = c(0,0,0,rep(0,5))),
      include.rownames = FALSE)





load("psolist.RData")

k1 <- 9001
k2 <- 10001
outdat <- data.frame(logpost = psolist[["pois"]][["full"]][[1]][[3]]$maxes[k1:k2], iter = k1:k2)

lnormfullplot <- ggplot(outdat) + geom_line(aes(iter, logpost)) +
  theme(axis.text = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 11)) +
  xlab("iteration") + ylab("log posterior")

options(scipen = -2)
ht <- 5
wd <- 5
ggsave("poisfullplot.png", lnormfullplot, width = wd, height = ht)


k1 <- 8000
k2 <- 10001
plot(ts(psolist[["lnorm"]][["full"]][[1]][[1]]$maxes[k1:k2]))


load("psolist.RData")
load("bbpsoxpmclist.RData")
load("atbbpsoxpmc1list.RData")
load("atbbpsoxpmc2list.RData")

k1 <- 200
k2 <- 1000
maxiter <- rbind(data.frame(logpost = psolist[["pois"]][["full"]][[1]][[1]]$maxes[k1:k2],
                            iter = k1:k2, alg = "pso"),
                 data.frame(logpost = bbpsoxpmclist[["pois"]][["full"]][[1]][[1]]$maxes[k1:k2],
                            iter = k1:k2, alg = "bbpsoxpmc"),
                 data.frame(logpost = atbbpsoxpmc1list[["pois"]][["full"]][[1]][[1]]$maxes[k1:k2],
                            iter = k1:k2, alg = "atbbpsoxpmc1"),
                 data.frame(logpost = atbbpsoxpmc2list[["pois"]][["full"]][[1]][[1]]$maxes[k1:k2],
                            iter = k1:k2, alg = "atbbpsoxpmc2"))

ggplot(maxiter) + geom_line(aes(iter, logpost, color = alg))


k1 <- 200
k2 <- 1000
maxiter <- rbind(data.frame(logpost = psolist[["lnorm"]][["full"]][[1]][[2]]$maxes[k1:k2],
                            iter = k1:k2, alg = "pso"),
                 data.frame(logpost = bbpsoxpmclist[["lnorm"]][["full"]][[1]][[2]]$maxes[k1:k2],
                            iter = k1:k2, alg = "bbpsoxpmc"),
                 data.frame(logpost = atbbpsoxpmc1list[["lnorm"]][["full"]][[1]][[2]]$maxes[k1:k2],
                            iter = k1:k2, alg = "atbbpsoxpmc1"),
                 data.frame(logpost = atbbpsoxpmc2list[["lnorm"]][["full"]][[1]][[2]]$maxes[k1:k2],
                            iter = k1:k2, alg = "atbbpsoxpmc2"))

ggplot(maxiter) + geom_line(aes(iter, logpost, color = alg))



load("psoout.RData")


psosum <- ddply(subset(psoout, iteration == 1000),
                .(obj, ndelta, algorithm, nbhd, nsubs), summarise,
                mean=mean(maxes), median=median(maxes), sd=sd(maxes),
                min=min(maxes), q025=quantile(maxes, 0.1),
                q25=quantile(maxes, 0.25), q75=quantile(maxes, 0.75),
                q975=quantile(maxes, 0.9), max=max(maxes))
psosum$nbhd <- as.character(psosum$nbhd)
psosum$algorithm <- mapvalues(psosum$algorithm,
          from = c("pso", "bbpso-mc", "bbpsoxp-mc", "at-bbpso-mc-3", "at-bbpsoxp-mc-3",
                   "at-bbpso-mc-5", "at-bbpsoxp-mc-5", "at-bbpso-mc-Inf", "at-bbpsoxp-mc-Inf"),
          to = c("PSO", "BBPSO-MC", "BBPSOxp-MC", "AT-BBPSO-MC-3", "AT-BBPSOxp-MC-3",
                 "AT-BBPSO-MC-5", "AT-BBPSOxp-MC-5", "AT-BBPSO-MC-Inf", "AT-BBPSOxp-MC-Inf"))
psosum$obj <- as.character(psosum$obj)
psosum$algorithm <- factor(psosum$algorithm,
                           levels = c("PSO", "BBPSO-MC", "BBPSOxp-MC",
                                      "AT-BBPSO-MC-3", "AT-BBPSOxp-MC-3",
                                      "AT-BBPSO-MC-5", "AT-BBPSOxp-MC-5",
                                      "AT-BBPSO-MC-Inf", "AT-BBPSOxp-MC-Inf"))

psosum$q025[psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "PSO" & psosum$nbhd == "global"] <- NA
psosum$q25[psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "PSO" & psosum$nbhd == "global"] <- NA
psosum$q75[psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "PSO" & psosum$nbhd == "global"] <- NA
psosum$q975[psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "PSO" & psosum$nbhd == "global"] <- NA
psosum$median[psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "PSO" & psosum$nbhd == "global"] <- NA
psosum$q025  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "BBPSO-MC" &
              psosum$nbhd == "global" & psosum$obj == "iid"] <- NA
psosum$q25   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "BBPSO-MC" &
              psosum$nbhd == "global" & psosum$obj == "iid"] <- NA
psosum$q75   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "BBPSO-MC" &
              psosum$nbhd == "global" & psosum$obj == "iid"] <- NA
psosum$q975  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "BBPSO-MC" &
              psosum$nbhd == "global" & psosum$obj == "iid"] <- NA
psosum$median[psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "BBPSO-MC" &
              psosum$nbhd == "global" & psosum$obj == "iid"] <- NA
psosum$q025  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-3" &
              psosum$nbhd == "global" & psosum$obj == "iid"] <- NA
psosum$q25   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-3" &
              psosum$nbhd == "global" & psosum$obj == "iid"] <- NA
psosum$q75   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-3" &
              psosum$nbhd == "global" & psosum$obj == "iid"] <- NA
psosum$q975  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-3" &
              psosum$nbhd == "global" & psosum$obj == "iid"] <- NA
psosum$median[psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-3" &
              psosum$nbhd == "global" & psosum$obj == "iid"] <- NA
psosum$q025  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-5" &
              psosum$nbhd == "global" & psosum$obj == "iid"] <- NA
psosum$q25   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-5" &
              psosum$nbhd == "global" & psosum$obj == "iid"] <- NA
psosum$q75   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-5" &
              psosum$nbhd == "global" & psosum$obj == "iid"] <- NA
psosum$q975  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-5" &
              psosum$nbhd == "global" & psosum$obj == "iid"] <- NA
psosum$median[psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-5" &
              psosum$nbhd == "global" & psosum$obj == "iid"] <- NA
psosum$q025  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-Inf" &
              psosum$nbhd == "global" & psosum$obj == "iid"] <- NA
psosum$q25   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-Inf" &
              psosum$nbhd == "global" & psosum$obj == "iid"] <- NA
psosum$q75   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-Inf" &
              psosum$nbhd == "global" & psosum$obj == "iid"] <- NA
psosum$q975  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-Inf" &
              psosum$nbhd == "global" & psosum$obj == "iid"] <- NA
psosum$median[psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-Inf" &
              psosum$nbhd == "global" & psosum$obj == "iid"] <- NA

r <- 30
nsubs <- unique(psosum$nsubs)
objs <- c("full", "iid")

## paper plots
maxlist <- list()
for(j in 1:3){
  maxlist[[j]] <- list()
  for(k in 1:2){
    temp <- subset(psosum, obj == objs[k] & nsubs == nsubs[j] & ndelta == r)  
    maxplot <- ggplot(temp) +
      geom_boxplot(aes(x=algorithm, lower = q25, upper = q75, middle = median, ymin = q025,
                       ymax = q975), stat="identity") +
      facet_wrap(nbhd~nsubs, ncol = 1,
                 labeller = label_bquote(paste(.(nbhd), " , M = ", .(nsubs)))) +
      theme(axis.text = element_text(angle = 45, hjust = 1),
            strip.text = element_text(size = 11)) +
      scale_x_discrete(drop=FALSE) + ggtitle(paste(objs[k], "Poisson model"))
    maxlist[[j]][[k]] <- maxplot
    options(scipen = -2)
    ht <- 6
    wd <- 3
    ggsave(paste("maxplot", j, k, ".png", sep=""), maxplot, width = wd, height = ht)
  }
}


## poster plots
psosum <- ddply(subset(psoout, iteration == 1000),
                .(obj, ndelta, algorithm, nbhd, nsubs), summarise,
                mean=mean(maxes), median=median(maxes), sd=sd(maxes),
                min=min(maxes), q025=quantile(maxes, 0.1),
                q25=quantile(maxes, 0.25), q75=quantile(maxes, 0.75),
                q975=quantile(maxes, 0.9), max=max(maxes))
psosum$nbhd <- as.character(psosum$nbhd)
psosum$algorithm <- mapvalues(psosum$algorithm,
          from = c("pso", "bbpso-mc", "bbpsoxp-mc", "at-bbpso-mc-3", "at-bbpsoxp-mc-3",
                   "at-bbpso-mc-5", "at-bbpsoxp-mc-5", "at-bbpso-mc-Inf", "at-bbpsoxp-mc-Inf"),
          to = c("PSO", "BBPSO-MC", "BBPSOxp-MC", "AT-BBPSO-MC-3", "AT-BBPSOxp-MC-3",
                 "AT-BBPSO-MC-5", "AT-BBPSOxp-MC-5", "AT-BBPSO-MC-Inf", "AT-BBPSOxp-MC-Inf"))
psosum$obj <- as.character(psosum$obj)
psosum$algorithm <- factor(psosum$algorithm,
                           levels = c("PSO", "BBPSO-MC", "BBPSOxp-MC",
                                      "AT-BBPSO-MC-3", "AT-BBPSOxp-MC-3",
                                      "AT-BBPSO-MC-5", "AT-BBPSOxp-MC-5",
                                      "AT-BBPSO-MC-Inf", "AT-BBPSOxp-MC-Inf"))
psosum$q025  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "PSO" &
              psosum$nbhd == "global"] <- NA
psosum$q25   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "PSO" &
              psosum$nbhd == "global"] <- NA
psosum$q75   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "PSO" &
              psosum$nbhd == "global"] <- NA
psosum$q975  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "PSO" &
              psosum$nbhd == "global"] <- NA
psosum$median[psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "PSO" &
              psosum$nbhd == "global"] <- NA
psosum$q025  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "PSO" &
              psosum$nbhd == "ring-1"] <- NA
psosum$q25   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "PSO" &
              psosum$nbhd == "ring-1"] <- NA
psosum$q75   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "PSO" &
              psosum$nbhd == "ring-1"] <- NA
psosum$q975  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "PSO" &
              psosum$nbhd == "ring-1"] <- NA
psosum$median[psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "PSO" &
              psosum$nbhd == "ring-1"] <- NA
psosum$q025  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "BBPSO-MC" &
              psosum$nbhd == "global" ] <- NA
psosum$q25   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "BBPSO-MC" &
              psosum$nbhd == "global" ] <- NA
psosum$q75   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "BBPSO-MC" &
              psosum$nbhd == "global" ] <- NA
psosum$q975  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "BBPSO-MC" &
              psosum$nbhd == "global" ] <- NA
psosum$median[psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "BBPSO-MC" &
              psosum$nbhd == "global" ] <- NA
psosum$q025  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-3" &
              psosum$nbhd == "global" ] <- NA
psosum$q25   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-3" &
              psosum$nbhd == "global" ] <- NA
psosum$q75   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-3" &
              psosum$nbhd == "global" ] <- NA
psosum$q975  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-3" &
              psosum$nbhd == "global" ] <- NA
psosum$median[psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-3" &
              psosum$nbhd == "global" ] <- NA
psosum$q025  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-5" &
              psosum$nbhd == "global" ] <- NA
psosum$q25   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-5" &
              psosum$nbhd == "global" ] <- NA
psosum$q75   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-5" &
              psosum$nbhd == "global" ] <- NA
psosum$q975  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-5" &
              psosum$nbhd == "global" ] <- NA
psosum$median[psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-5" &
              psosum$nbhd == "global" ] <- NA
psosum$q025  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-Inf" &
              psosum$nbhd == "global" ] <- NA
psosum$q25   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-Inf" &
              psosum$nbhd == "global" ] <- NA
psosum$q75   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-Inf" &
              psosum$nbhd == "global" ] <- NA
psosum$q975  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-Inf" &
              psosum$nbhd == "global" ] <- NA
psosum$median[psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-Inf" &
              psosum$nbhd == "global" ] <- NA
psosum$q025  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "BBPSO-MC" &
              psosum$nbhd == "ring-1" ] <- NA
psosum$q25   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "BBPSO-MC" &
              psosum$nbhd == "ring-1" ] <- NA
psosum$q75   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "BBPSO-MC" &
              psosum$nbhd == "ring-1" ] <- NA
psosum$q975  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "BBPSO-MC" &
              psosum$nbhd == "ring-1" ] <- NA
psosum$median[psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "BBPSO-MC" &
              psosum$nbhd == "ring-1" ] <- NA
psosum$q025  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-3" &
              psosum$nbhd == "ring-1" ] <- NA
psosum$q25   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-3" &
              psosum$nbhd == "ring-1" ] <- NA
psosum$q75   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-3" &
              psosum$nbhd == "ring-1" ] <- NA
psosum$q975  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-3" &
              psosum$nbhd == "ring-1" ] <- NA
psosum$median[psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-3" &
              psosum$nbhd == "ring-1" ] <- NA
psosum$q025  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-5" &
              psosum$nbhd == "ring-1" ] <- NA
psosum$q25   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-5" &
              psosum$nbhd == "ring-1" ] <- NA
psosum$q75   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-5" &
              psosum$nbhd == "ring-1" ] <- NA
psosum$q975  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-5" &
              psosum$nbhd == "ring-1" ] <- NA
psosum$median[psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-5" &
              psosum$nbhd == "ring-1" ] <- NA
psosum$q025  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-Inf" &
              psosum$nbhd == "ring-1" ] <- NA
psosum$q25   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-Inf" &
              psosum$nbhd == "ring-1" ] <- NA
psosum$q75   [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-Inf" &
              psosum$nbhd == "ring-1" ] <- NA
psosum$q975  [psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-Inf" &
              psosum$nbhd == "ring-1" ] <- NA
psosum$median[psosum$ndelta == 30 & psosum$nsubs == 3142 & psosum$algorithm == "AT-BBPSO-MC-Inf" &
              psosum$nbhd == "ring-1" ] <- NA
psosum$obj[psosum$obj == "full"] <- "Full"
psosum$obj[psosum$obj == "iid"] <- "IID"
psosum$nbhd[psosum$nbhd == "ring-1"] <- "Ring-1"
psosum$nbhd[psosum$nbhd == "global"] <- "Global"


temp <- subset(psosum, nsubs == max(nsubs) & ndelta == 30)
maxplot <- ggplot(temp) +
  geom_boxplot(aes(x=algorithm, lower = q25, upper = q75, middle = median, ymin = q025,
                   ymax = q975), stat="identity") +
  facet_wrap(obj~nbhd, scales = "fixed", nrow = 1,
             labeller = label_bquote(paste(.(obj), " model;   ", .(nbhd), " nbhd"))) +
  theme(axis.text = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 11)) +
  scale_x_discrete(drop=FALSE) + xlab("PSO algorithm") + ylab("log posterior") +
  ggtitle("Log posterior at the mode estimated by each PSO algorithm over ten replications; r = 30 random effects") 



options(scipen = -2)
ht <- 3
wd <- 12
ggsave(paste("maxplotpost", ".png", sep=""), maxplot, width = wd, height = ht)


load("mcmcout.RData")
load("mcmcoutimh.RData")
load("mcmcoutstan.RData")
rownames(mcmcout) <- NULL
rownames(mcmcoutstan) <- NULL
rownames(mcmcoutimh) <- NULL
mcmcout <- subset(mcmcout, !(alg %in% c("IMH", "IMHwGibbs")))
mcmcout <- rbind(mcmcout, mcmcoutimh)

mcmcoutsub <- subset(mcmcout, nsub == 3142)[,c(1,2,3,4,5,7,10)]
mcmcoutsub$nefftime <- mcmcoutsub$itertime/mcmcoutsub$neff*10000
mcmcoutstansub <- mcmcoutstan[,c(1,2,3,4,5,7,10)]
mcmcoutstansub$nefftime <- mcmcoutstansub$itertime/mcmcoutstansub$neff*1000
mcmcoutsub$df[mcmcoutsub$df == 0] <- ""
mcmcoutstansub$df[mcmcoutstansub$df == 0] <- ""
mcmcoutsub$alg <- mapvalues(mcmcoutsub$alg, c("IMHwGibbs", "blockRWwG"), c("IMHwG", "B-RWwG"))
mcmcoutsub$alg2 <- mcmcoutsub$alg
mcmcoutstansub$alg2 <- mcmcoutstan$alg
mcmcoutsub2 <- rbind(mcmcoutsub, mcmcoutstansub)
mcmcoutsub2$alg2 <- factor(mcmcoutsub2$alg2, levels(mcmcoutsub2$alg2)[c(6,3,5,4,2,1)])

mcmcoutmelt <- melt(mcmcoutsub2, id.vars = c("model", "ranef", "ndelta", "alg", "df", "alg2"))
poisiidneff <- dcast(mcmcoutmelt, ndelta ~ alg2,
                        subset = .(model == "Poisson" & ranef == "iid" & variable == "neff"))
poisiidnefftime <- dcast(mcmcoutmelt, ndelta ~ alg2,
                        subset = .(model == "Poisson" & ranef == "iid" & variable == "nefftime"))
poisfullneff <- dcast(mcmcoutmelt, ndelta ~ alg2,
                        subset = .(model == "Poisson" & ranef == "full" & variable == "neff"))
poisfullnefftime <- dcast(mcmcoutmelt, ndelta ~ alg2,
                        subset = .(model == "Poisson" & ranef == "full" & variable == "nefftime"))
lnormiidneff <- dcast(mcmcoutmelt, ndelta ~ alg2,
                        subset = .(model == "lognormal" & ranef == "iid" & variable == "neff"))
lnormiidnefftime <- dcast(mcmcoutmelt, ndelta ~ alg2,
                        subset = .(model == "lognormal" & ranef == "iid" & variable == "nefftime"))
lnormfullneff <- dcast(mcmcoutmelt, ndelta ~ alg2,
                        subset = .(model == "lognormal" & ranef == "full" & variable == "neff"))
lnormfullnefftime <- dcast(mcmcoutmelt, ndelta ~ alg2,
                        subset = .(model == "lognormal" & ranef == "full" & variable == "nefftime"))

print(xtable(rbind(poisiidneff, poisfullneff), digits = 0), include.rownames = FALSE)

print(xtable(rbind(lnormiidneff, lnormfullneff), digits = 0), include.rownames = FALSE)

print(xtable(rbind(poisiidnefftime, poisfullnefftime), digits = c(0,0,rep(0,5))),
      include.rownames = FALSE)

print(xtable(rbind(lnormiidnefftime, lnormfullnefftime), digits = c(0,0,0,rep(0,5))),
      include.rownames = FALSE)











iidneff <- matrix(mcmcoutsumiid$neff, ncol = 3)
colnames(iidneff) <- c(5, 15, 25)
rownames(iidneff) <- c("PSO-IMHwGibbs", "PSO-IMH", "blockRWwG", "RWwG")
iidneff <- iidneff[c(2,1,4,3),]

iidnefftime <- matrix(mcmcoutsumiid$nefftime, ncol = 3)
colnames(iidnefftime) <- c(5, 15, 25)
rownames(iidnefftime) <- c("PSO-IMHwGibbs", "PSO-IMH", "blockRWwG", "RWwG")
iidnefftime <- iidnefftime[c(2,1,4,3),]

iidtab <- cbind(iidneff, iidnefftime)
xtable(iidtab)

mcmcoutsumfull <- subset(mcmcoutsub, model == "full")
mcmcoutsumfull <- mcmcoutsumfull[order(mcmcoutsumfull$ndelta),]

fullneff <- matrix(mcmcoutsumfull$neff[-c(1,2)], ncol = 3)
colnames(fullneff) <- c(5, 15, 25)
fullneff <- rbind(rbind(c(0,0,0), c(0,0,0)), fullneff)
fullneff[c(1,2),1] <- mcmcoutsumfull$neff[1:2]
rownames(fullneff) <- c("PSO-IMHwGibbs", "PSO-IMH", "B-RWwG", "RWwG")
fullneff <- fullneff[c(2,1,4,3),]

fullnefftime <- matrix(mcmcoutsumfull$nefftime[-c(1,2)], ncol = 3)
colnames(fullnefftime) <- c(5, 15, 25)
fullnefftime <- rbind(rbind(c(0,0,0), c(0,0,0)), fullnefftime)
fullnefftime[c(1,2),1] <- mcmcoutsumfull$nefftime[1:2]
rownames(fullnefftime) <- c("PSO-IMHwGibbs", "PSO-IMH", "B-RWwG", "RWwG")
fullnefftime <- fullnefftime[c(2,1,4,3),]

fulltab <- cbind(fullneff, fullnefftime)
xtable(fulltab)


xtable(iidneff)
xtable(iidnefftime)
xtable(fullneff)
xtable(fullnefftime)


load("mcmcout3.RData")

mcmcout3sumfull <- subset(mcmcout3, model == "full")[3:1,c(1,6,9)]

mcmcout3sumiid <- subset(mcmcout3, model == "iid")[3:1,c(1,6,9)]
