library(plyr)
library(dplyr)
library(reshape2)
library(xtable)

niter <- 500

psoout <- read.csv("psoout.csv")

convpercent <- function(x, cutoff){
  return(mean(abs(x)<=cutoff))
}

psosum <- ddply(psoout, .(obj, type, nbhd, time), summarise, mean=mean(-logpost),
                median=median(-logpost), sd=sd(logpost), min=min(-logpost),
                q05=quantile(-logpost, 0.05), q25=quantile(-logpost, 0.25),
                q75=quantile(-logpost, 0.75), q95=quantile(-logpost, 0.95),
                max=max(-logpost), con1=convpercent(-logpost, 0.01),
                con2=convpercent(-logpost, 0.0001))

psosum$obj[psosum$obj > 2] <- psosum$obj[psosum$obj > 2] - 1
psosum$obj[psosum$obj > 5] <- psosum$obj[psosum$obj > 5] - 1

convtime <- function(lpost, time, cutoff){
  idx <- which(abs(lpost) <= cutoff)
  if(length(idx)>0){
    return(min(time[idx]))
  } else {
    return(Inf)
  }
}

psosum2 <- ddply(psoout, .(obj, type, nbhd, rep), summarise,
                 con1 = convtime(logpost, time, 0.01),
                 con2 = convtime(logpost, time, 0.0001))
psosum2$obj[psosum2$obj > 2] <- psosum2$obj[psosum2$obj > 2] - 1
psosum2$obj[psosum2$obj > 5] <- psosum2$obj[psosum2$obj > 5] - 1


psosum3 <- ddply(psosum2, .(obj, type, nbhd), summarise,
                 median1 = median(con1),
                 median2 = median(con2))

sumlast <- subset(psosum, time == niter)

psosumout <- merge(sumlast, psosum3, by = c("obj", "type", "nbhd"))

save(psosumout, file = "psosumout.RData")



sumlast <- psosumout[,c(1:3, 5, 7, 14,16,15,17)]

for(k in 1:6){
sumlastk <- subset(sumlast, subset=(obj == k))
sumlastk <- sumlastk[,c(2,3,4,5,6,7)]
tabout <- cbind(rbind(subset(sumlastk, nbhd == "global" & type == "PSO"),
                      subset(sumlastk, nbhd == "global" & type %in% c("BBPSO-MC", "BBPSOxp-MC")),
                      subset(sumlastk, nbhd == "global" & type %in%
                                     paste("AT-BBPSO-MC",
                                           c(rep(1, 4), rep(3, 4), rep(5, 4), rep("Inf", 4)),
                                           rep(c(0.1, 0.3, 0.5, 0.7), 4), sep = "-")),
                      subset(sumlastk, nbhd == "global" & type %in%
                                     paste("AT-BBPSOxp-MC",
                                           c(rep(1, 4), rep(3, 4), rep(5, 4), rep("Inf", 4)),
                                           rep(c(0.1, 0.3, 0.5, 0.7), 4), sep = "-")),
                      subset(sumlastk, nbhd == "global"& type %in%
                                       paste("DI-PSO", c(rep(50, 3)),
                                             c(rep(c(1,2,4), 1)), sep = "-")),
                      subset(sumlastk, nbhd == "global"& type %in%
                                       paste("DI-PSO", c(rep(100, 3), rep(200, 3)),
                                             c(rep(c(1,2,4), 2)), sep = "-")),
                      subset(sumlastk, nbhd == "global"& type %in%
                                       paste("AT-PSO", c(0.1, 0.3, 0.5, 0.7), sep = "-")))[,-2],
                rbind(subset(sumlastk, nbhd == "ring-3" & type == "PSO"),
                      subset(sumlastk, nbhd == "ring-3" & type %in% c("BBPSO-MC", "BBPSOxp-MC")),
                      subset(sumlastk, nbhd == "ring-3" & type %in%
                                       paste("AT-BBPSO-MC",
                                           c(rep(1, 4), rep(3, 4), rep(5, 4), rep("Inf", 4)),
                                           rep(c(0.1, 0.3, 0.5, 0.7), 4), sep = "-")),
                      subset(sumlastk, nbhd == "ring-3" & type %in%
                                     paste("AT-BBPSOxp-MC",
                                           c(rep(1, 4), rep(3, 4), rep(5, 4), rep("Inf", 4)),
                                           rep(c(0.1, 0.3, 0.5, 0.7), 4), sep = "-")),
                      subset(sumlastk, nbhd == "ring-3"& type %in%
                                       paste("DI-PSO", c(rep(50, 3)),
                                             c(rep(c(1,2,4), 1)), sep = "-")),
                      subset(sumlastk, nbhd == "ring-3"& type %in%
                                       paste("DI-PSO", c(rep(100, 3), rep(200, 3)),
                                             c(rep(c(1,2,4), 2)), sep = "-")),
                      subset(sumlastk, nbhd == "ring-3"& type %in%
                                       paste("AT-PSO", c(0.1, 0.3, 0.5, 0.7), sep = "-")))[-c(1,2)],
                rbind(subset(sumlastk, nbhd == "ring-1" & type == "PSO"),
                      subset(sumlastk, nbhd == "ring-1" & type %in% c("BBPSO-MC", "BBPSOxp-MC")),
                      subset(sumlastk, nbhd == "ring-1" & type %in%
                                     paste("AT-BBPSO-MC",
                                           c(rep(1, 4), rep(3, 4), rep(5, 4), rep("Inf", 4)),
                                           rep(c(0.1, 0.3, 0.5, 0.7), 4), sep = "-")),
                      subset(sumlastk, nbhd == "ring-1" & type %in%
                                     paste("AT-BBPSOxp-MC",
                                           c(rep(1, 4), rep(3, 4), rep(5, 4), rep("Inf", 4)),
                                           rep(c(0.1, 0.3, 0.5, 0.7), 4), sep = "-")),
                      subset(sumlastk, nbhd == "ring-1"& type %in%
                                       paste("DI-PSO", c(rep(50, 3)),
                                             c(rep(c(1,2,4), 1)), sep = "-")),
                      subset(sumlastk, nbhd == "ring-1"& type %in%
                                       paste("DI-PSO", c(rep(100, 3), rep(200, 3)),
                                             c(rep(c(1,2,4), 2)), sep = "-")),
                      subset(sumlastk, nbhd == "ring-1"& type %in%
                                       paste("AT-PSO", c(0.1, 0.3, 0.5, 0.7), sep = "-")))[-c(1,2)])
tabout$type <- mapvalues(tabout$type, from = paste("AT-BBPSO-MC", c(rep(1, 4), rep(3, 4), rep(5, 4), rep("Inf", 4)), rep(c(0.1, 0.3, 0.5, 0.7), 4), sep = "-"),
          to = paste("$df = ", c(rep(paste(1, ",\\enspace", sep=""), 4), rep(paste(3, ",\\enspace", sep=""), 4), rep(paste(5, ",\\enspace", sep=""), 4), rep("\\infty,", 4)), "$ $R^* =", rep(c(0.1, 0.3, 0.5, 0.7), 4), "$", sep = ""))
tabout$type <- mapvalues(tabout$type, from = paste("AT-BBPSOxp-MC", c(rep(1, 4), rep(3, 4), rep(5, 4), rep("Inf", 4)), rep(c(0.1, 0.3, 0.5, 0.7), 4), sep = "-"),
                              to = paste("$df = ", c(rep(paste(1, ",\\enspace", sep=""), 4), rep(paste(3, ",\\enspace", sep=""), 4), rep(paste(5, ",\\enspace", sep=""), 4), rep("\\infty,", 4)), "$ $R^* =", rep(c(0.1, 0.3, 0.5, 0.7), 4), "$", sep = ""))
tabout$type <- mapvalues(tabout$type,
                              from = paste("DI-PSO", c(rep(50, 3), rep(100, 3), rep(200, 3)), rep(c(1, 2, 4), 3), sep = "-"),
                              to = paste("$\\alpha = ", c(rep(paste(50, ",\\enspace", sep=""), 3), rep("100,", 3), rep("200,", 3)), "$ $\\beta =", rep(c(1, 2, 4), 3), "$", sep = ""))
tabout$type <- mapvalues(tabout$type,
                              from = paste("AT-PSO", c(0.1, 0.3, 0.5, 0.7), sep="-"),
                              to = paste("$R^* = ", c(0.1, 0.3, 0.5, 0.7), "$", sep=""))
tabout$type <- mapvalues(tabout$type,
                              from = c("PSO", "BBPSO-MC", "BBPSOxp-MC"),
                              to = c("\\multicolumn{1}{l|}{PSO}", "\\multicolumn{1}{l|}{BBPSO-MC}",
                                     "\\multicolumn{1}{l|}{BBPSOxp-MC}"))
print(xtable(tabout), include.rownames=FALSE, sanitize.text.function=identity)
}





library(ggplot2)

atplotout <- subset(psoout, (obj %in% c(1, 8) & nbhd == "ring-1" & type == "AT-PSO-0.5" &
                             time > 0 & rep == 40))
atplotout <- rbind(atplotout, data.frame(obj = 1, logpost = 0, argnorm = 0, time = 1:niter, type = "DI-PSO-200-1", nbhd = "ring-1", rep = 1, inertias = 1/(1 + ((1:niter + 1)/200)^1)))
atplotout$obj[atplotout$obj == 8] <- 6
atplotout$Algorithm <- paste(atplotout$type, ", Obj = ", atplotout$obj, sep = "")
atplotout$Algorithm <- mapvalues(atplotout$Algorithm, "DI-PSO-200-1, Obj = 1", "DI-PSO-200-1")
atplotout$Algorithm <- factor(atplotout$Algorithm, unique(atplotout$Algorithm)[c(3,1,2)])
inertiaplot <- qplot(time, inertias, data = atplotout, geom = "line", linetype = Algorithm) +
  xlab("iteration") + ylab("inertia")
ht <- 3
wd <- 6
ggsave("inertiaplot.png", inertiaplot, width = wd, height = ht)

