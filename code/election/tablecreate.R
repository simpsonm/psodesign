library(xtable)
library(plyr)
library(ggplot2)
library(reshape2)


accpsoout <- read.csv("accpollpsoout.csv")


accout <- read.csv("accout.csv")

accout$pso <-mapvalues(accout$pso,
                       c("AT-PSO-0.5-0.1", "AT-BBPSO-MC-5-0.5-0.1", "AT-BBPSOxp-MC-5-0.5-0.1"),
                       c("AT-PSO", "AT-BBPSO-MC", "AT-BBPSOxp-MC"))

small <- subset(accout, model == "small")[,c(6,7,8,5,9)-2]
tempmelt <- melt(small, id.vars = c("pso", "nbhd", "mcmc", "niter"))
smallcast <- dcast(tempmelt, pso + nbhd ~ mcmc + niter)

poll <- subset(accout, model == "poll")[,c(6,7,8,5,9)-2]
tempmelt <- melt(poll, id.vars = c("pso", "nbhd", "mcmc", "niter"))
pollcast <- dcast(tempmelt, pso + nbhd ~ mcmc + niter)

subset(smallcast, pso == "PSO")




psoout <- read.csv("pollpsoout.csv")

sumlast <- subset(psoout, iteration %in% c(1000, 5000, 10000))

psosum <- ddply(sumlast, .(model, algorithm, nbhd, iteration), summarise,
                mean = mean(maxes), median = median(maxes), sd = sd(maxes),
                min=min(maxes), q10=quantile(maxes, 0.1), q25=quantile(maxes, 0.25),
                q75=quantile(maxes, 0.75), q90=quantile(maxes, 0.9), max=max(maxes))

psosum$algorithm <-
  mapvalues(psosum$algorithm,
            c("AT-PSO-0.5-0.1", "AT-BBPSO-MC-5-0.5-0.1", "AT-BBPSOxp-MC-5-0.5-0.1"),
            c("AT-PSO", "AT-BBPSO-MC", "AT-BBPSOxp-MC"))

algids <- c(7,4,5,6,1,2,3)
cutoff <- -10000

small <- subset(psosum, model == "small")[,c(4,5,6,8)-2]
small$mean <- round(small$mean - max(subset(small, algorithm == "PSO")$mean), 2)
small$sd[small$mean < cutoff] <- NA
small$mean[small$mean < cutoff] <- NA
small <- rbind(subset(small, nbhd == "global")[algids,],
               subset(small, nbhd == "ring-3")[algids,],
               subset(small, nbhd == "ring-1")[algids,])
small


