library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)

load("psoout.RData")

maxiter <- max(psoout$iteration)

psoend <- subset(psoout, subset = iteration == maxiter)

psosum <- ddply(psoend, .(model, ranef, ndelta, algorithm, init, nbhd), summarise,
                mean=mean(maxes), sd=sd(maxes), min=min(maxes), q05=quantile(maxes, 0.05),
                q10=quantile(maxes, 0.1),
                q25=quantile(maxes, 0.25), median=median(maxes), q75=quantile(maxes, 0.75),
                q90=quantile(maxes, 0.90),
                q95=quantile(maxes, 0.95), max=max(maxes))
psosum$randelta <- factor(paste(psosum$ranef, ": ", psosum$ndelta, sep = ""),
                          levels = c(paste("iid: ", c(10, 30), sep = ""),
                                     paste("full: ", c(5, 15), sep = "")))
psosum$algorithm <- factor(psosum$algorithm, levels(psosum$algorithm)[c(1:4,6,8,10,5,7,9,11)])
psosum$init <- mapvalues(psosum$init, from = c("zero", "bfgs"), to = c("no BFGS", "BFGS"))
psosum <- subset(psosum, algorithm %in% levels(psosum$algorithm)[c(1:4,6,8,10)])
psosum$algorithm <- mapvalues(psosum$algorithm,
                              from = c("BBPSOxp-MC", "AT-PSO-0.3-0.1", "AT-PSO-0.5-0.1",
                                       "AT-BBPSOxp-MC-1-0.3-0.1", "AT-BBPSOxp-MC-1-0.5-0.1"),
                              to = c("BBPSO", "AT-PSO-0.3", "AT-PSO-0.5",
                                       "AT-BBPSO-0.3", "AT-BBPSO-0.5"))

models <- c("pois", "lnorm")
models2 <- c("Poisson", "Lognormal")
ndeltas <- c(10,30,5,15)

options(scipen = -2)
ht <- 8
wd <- 4
for(i in 1:2){
  for(j in 1:4){
    ranef <- ifelse(j<3, "iid", "full")
    modelname <- paste(models2[i], " ", ranef, ", r = ", ndeltas[j], sep="")
    maxplot <- ggplot(subset(psosum, (model == models[i] & ndelta == ndeltas[j]))) +
      geom_boxplot(aes(x = algorithm, lower = q25, upper = q75, middle = median,
                       ymin = q10, ymax = q90), stat = "identity") +
      facet_grid(init~nbhd) +
      theme(axis.text = element_text(angle = 45, hjust = 1),
            strip.text = element_text(size = 11)) +
      ggtitle(modelname)
    ggsave(paste("maxplot", i, j, ".png", sep=""), maxplot, width = wd, height = ht)
  }
}


k
