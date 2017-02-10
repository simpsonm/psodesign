library(plyr)
library(dplyr)
library(reshape2)
library(xtable)
library(ggplot2)
library(gridExtra)

psoout <- read.csv("psosimsout.csv")
psoout <- subset(psoout, rep <= 40)

convpercent <- function(x, cutoff){
  return(mean(abs(x)<=cutoff))
}

niter <- max(psoout$time)
psoout$k1 <- abs(psoout$logpost)<=0.01
psoout$k2 <- abs(psoout$logpost)<=0.0001
psoout$obj[psoout$obj > 2] <- psoout$obj[psoout$obj > 2] - 1
psoout$obj[psoout$obj > 5] <- psoout$obj[psoout$obj > 5] - 1


psotimesum <- ddply(psoout, .(obj, algid, nbhd, time, type, parset, CF, style), summarise, mean=mean(logpost),
                median=median(logpost), sd=sd(logpost), min=min(logpost),
                q05=quantile(logpost, 0.05), q25=quantile(logpost, 0.25),
                q75=quantile(logpost, 0.75), q95=quantile(logpost, 0.95),
                max=max(logpost), con1=convpercent(logpost, 0.01),
                con2=convpercent(logpost, 0.0001))

save(psotimesum, file = "psotimesum.Rdata")

psorepsum <- ddply(psoout, .(obj, algid, nbhd, rep, type, parset, CF, style), summarise,
                 kk1=niter-sum(k1)+1, kk2=niter-sum(k2)+1)

save(psorepsum, file = "psorepsum.Rdata")

sumlast <- subset(psotimesum, time == niter)

sumrep <- ddply(psorepsum, .(obj, algid, nbhd, type, parset, CF, style), summarise,
                k1 = median(kk1), k2 = median(kk2))

psosum <- merge(sumlast, sumrep, by=c("obj", "algid", "nbhd", "type", "parset", "CF", "style"))
psosum$k1[psosum$k1 > niter] <- Inf
psosum$k2[psosum$k2 > niter] <- Inf

save(psosum, file = "psosum.Rdata")

namefun <- function(x){
  if(x$style == "AT1"){
    style <- "AT3-"
  } else if(x$style == "AT2"){
    style <- "AT5-"
  } else if(x$style == "CI"){
    style <- ""
  } else {
    style <- paste(x$style, "-", sep="")
  }
  if(x$type == "BBPSO"){
    parset <- switch(x$parset, "", "xp")
  } else {
    parset <- x$parset
  }
  type <- paste(x$type, parset, sep = "")
  cf <- switch(x$CF, "-CF")
  out <- paste(style, type, cf, sep = "")
  return(out)
}

nbhdnamefun <- function(x){
  nbhd <- paste("n", x$nbhd, sep="")
  out <- switch(nbhd, n1 = "SS1", n3 = "SS3", n40 = "Global")
  return(out)
}

namesfun <- function(x){
  nbhd <- nbhdnamefun(x)
  alg <- namefun(x)
  return(c(alg, nbhd))
}

psosumsmol <- ddply(psosum, .(obj, algid, nbhd, mean, sd, con1, k1), namesfun)[,c(1,8,9,4,5,6,7)]
colnames(psosumsmol) <- c("Obj", "Algorithm", "Nbhd", "Mean", "SD", "$\\widehat{P}$", "$\\widehat{K}$")


algorder <- c("PSO1", "PSO2", "PSO1-CF", "PSO2-CF",
              "DI-PSO1", "DI-PSO2", "DI-PSO1-CF", "DI-PSO2-CF",
              "AT3-PSO1", "AT3-PSO2", "AT3-PSO1-CF", "AT3-PSO2-CF",
              "AT5-PSO1", "AT5-PSO2", "AT5-PSO1-CF", "AT5-PSO2-CF",
              "AT3-BBPSO", "AT3-BBPSOxp", "AT3-BBPSO-CF", "AT3-BBPSOxp-CF",
              "AT5-BBPSO", "AT5-BBPSOxp", "AT5-BBPSO-CF", "AT5-BBPSOxp-CF")
nbhdorder <- c("Global", "SS3", "SS1")

psosumsmol$Algorithm <- factor(psosumsmol$Algorithm, levels = algorder)

psosumsmol$Nbhd <- factor(psosumsmol$Nbhd, levels = nbhdorder)

psosumorder <- arrange(psosumsmol, Obj, Algorithm, Nbhd)

## this prints the tables, however it's missing the first header.
## paste the following line (without the ##) just after the \begin{tabular} line
## make sure to change OF1 to OF2, OF3, etc., in the those tables
## \multicolumn{1}{l}{\begin{tabular}{lr} OF1; & Nbhd: \end{tabular}} & \multicolumn{4}{c}{Global} & \multicolumn{4}{c}{SS3} & \multicolumn{4}{c}{SS1} \\
for(k in 1:6){
  ofk <- subset(psosumorder, Obj == k)[,-1]
  ## remove large values from table
  bigidx <- which(ofk$Mean > 10000)
  ofk$Mean[bigidx] <- NA
  ofk$SD[bigidx] <- NA
  ## pull out each type of nbhd and combine
  ofk1 <- subset(ofk, Nbhd == "Global")[,-2]
  ofk2 <- subset(ofk, Nbhd == "SS3")[,-2]
  ofk3 <- subset(ofk, Nbhd == "SS1")[,-2]
  ofkout <- cbind(ofk1, ofk2[,-1], ofk3[,-1])
  ## round
  ofkout2 <- ofkout
  ofkout2[,-1] <- round(signif(ofkout2[,-1], 5), 4)
  ## properly label \widehat{K}'s that are large
  ofkout2[ofkout2 == "Inf"] <- "$> 1000$"
  print(xtable(ofkout2, digits = c(0,0,rep(c(2,2,2,2),3)), align = "ll|rrrr|rrrr|rrrr",
               caption = paste("Simulation results for OF", k, ". See text for description", sep = ""),
               label = paste("tab:psosim", k, sep="")),
        include.rownames=FALSE, sanitize.text.function=identity, hline.after=c(-1,seq(0, 24, 4)),
        size = "\\scriptsize")
}




rr <- 2
algidstr <- "PSO-2-CF-AT2"
atplotout <- subset(psoout, (obj %in% c(1, 6) & nbhd == 3 & algid == algidstr &
                             time > 0 & rep == rr))[,c(1,4,5,12)]
atplotout <- rbind(atplotout, data.frame(obj = 1, time = 1:niter, algid = "DI-PSO", inertias = 1/(1 + ((1:niter + 1)/200)^1)))
atplotout$Algorithm <- mapvalues(atplotout$algid, algidstr, "AT2-PSO2-CF")
atplotout$Algorithm <- paste(atplotout$Algorithm, ", Obj = ", atplotout$obj, sep = "")
atplotout$Algorithm <- mapvalues(atplotout$Algorithm, "DI-PSO, Obj = 1", "DI-PSO")
atplotout$Algorithm <- factor(atplotout$Algorithm, unique(atplotout$Algorithm)[c(3,1,2)])
inertiaplot <- qplot(time, inertias, data = atplotout, geom = "line", color = Algorithm) +
  xlab("iteration") + ylab("inertia")

algidstr2 <- "BBPSO-1-CF-AT1"
atplotout2 <- subset(psoout, (obj %in% 1:6 & nbhd == 3 & algid == algidstr2 &
                             time > 0 & rep == rr))[,c(1,4,5,12)]
atplotout2$Obj <- paste("OF", atplotout2$obj, sep="")
scaleplot <- qplot(time, inertias, data = atplotout2, geom = "line", color = Obj) +
  xlab("iteration") + ylab("scale")

ht <- 4
wd <- 8
ggsave("inertiaplot.png", inertiaplot, width = wd, height = ht)

ht <- 4
wd <- 8
ggsave("scaleplot.png", scaleplot, width = wd, height = ht)
