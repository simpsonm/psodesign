load("psoout.RData")

library(dplyr)
library(reshape2)
library(plyr)

convpercent <- function(x, cutoff){
  return(mean(abs(x)<=cutoff))
}

psosum <- ddply(psoout, .(obj, type, nbhd, time), summarise, mean=mean(logpost),
                median=median(logpost), sd=sd(logpost), min=min(logpost),
                q05=quantile(logpost, 0.05), q25=quantile(logpost, 0.25),
                q75=quantile(logpost, 0.75), q95=quantile(logpost, 0.95),
                max=max(logpost), con1=convpercent(logpost, 0.01),
                con2=convpercent(logpost, 0.0001))
posum$obj[psosum$obj > 2] <- psosum$obj[psosum$obj > 2] - 1
psosum$obj[psosum$obj > 5] <- psosum$obj[psosum$obj > 5] - 1


## psoargsum <- ddply(psoout, .(obj, type, nbhd, time), summarise, mean=mean(argnorm), median=median(argnorm), sd=sd(argnorm), min=min(argnorm), q05=quantile(argnorm, 0.05), q25=quantile(argnorm, 0.25), q75=quantile(argnorm, 0.75), q95=quantile(argnorm, 0.95), max=max(argnorm))

save(psosum, file = "psosum.Rdata")
## save(psoargsum, file = "psoargsum.Rdata")

load("psosum.Rdata")
## load("psoargsum.Rdata")

psosum[,5:13] <- abs(psosum[,5:13])


library(xtable)

sumlast <- subset(psosum, time == niter)
colnames(sumlast)[c(2,3,5,6,7)] <- c("Algorithm", "Nbhd", "Mean", "Median", "SD")

sumlast <- sumlast[,c(1:3, 5, 7, 14:15)]

k <- 2

for(k in 1:6){
tabout <- cbind(rbind(subset(sumlast, subset=(obj == k & Nbhd == "global" & Algorithm %in% c("PSO", "BBPSO-MC", "BBPSOxp-MC"))),
      subset(sumlast, subset=(obj == k & Nbhd == "global" &
                              Algorithm %in% paste("AT-BBPSO-MC", c(rep(1, 4), rep(3, 4), rep(5, 4), rep("Inf", 4)), rep(c(0.1, 0.3, 0.5, 0.7), 4), sep = "-"))),
      subset(sumlast, subset=(obj == k & Nbhd == "global" &
                              Algorithm %in% paste("AT-BBPSOxp-MC", c(rep(1, 4), rep(3, 4), rep(5, 4), rep("Inf", 4)), rep(c(0.1, 0.3, 0.5, 0.7), 4), sep = "-"))),
      subset(sumlast, subset=(obj == k & Nbhd == "global"))[35 + 1:13,])[,c(2,4,5,6,7)],
rbind(subset(sumlast, subset=(obj == k & Nbhd == "ring-3" & Algorithm %in% c("PSO", "BBPSO-MC", "BBPSOxp-MC"))),
      subset(sumlast, subset=(obj == k & Nbhd == "ring-3" &
                              Algorithm %in% paste("AT-BBPSO-MC", c(rep(1, 4), rep(3, 4), rep(5, 4), rep("Inf", 4)), rep(c(0.1, 0.3, 0.5, 0.7), 4), sep = "-"))),
      subset(sumlast, subset=(obj == k & Nbhd == "ring-3" &
                              Algorithm %in% paste("AT-BBPSOxp-MC", c(rep(1, 4), rep(3, 4), rep(5, 4), rep("Inf", 4)), rep(c(0.1, 0.3, 0.5, 0.7), 4), sep = "-"))),
      subset(sumlast, subset=(obj == k & Nbhd == "ring-3"))[35 + 1:13,])[,c(4,5,6,7)],
rbind(subset(sumlast, subset=(obj == k & Nbhd == "ring-1" & Algorithm %in% c("PSO", "BBPSO-MC", "BBPSOxp-MC"))),
      subset(sumlast, subset=(obj == k & Nbhd == "ring-1" &
                              Algorithm %in% paste("AT-BBPSO-MC", c(rep(1, 4), rep(3, 4), rep(5, 4), rep("Inf", 4)), rep(c(0.1, 0.3, 0.5, 0.7), 4), sep = "-"))),
      subset(sumlast, subset=(obj == k & Nbhd == "ring-1" &
                              Algorithm %in% paste("AT-BBPSOxp-MC", c(rep(1, 4), rep(3, 4), rep(5, 4), rep("Inf", 4)), rep(c(0.1, 0.3, 0.5, 0.7), 4), sep = "-"))),
      subset(sumlast, subset=(obj == k & Nbhd == "ring-1"))[35 + 1:13,])[,c(4,5,6,7)])
tabout$Algorithm <- mapvalues(tabout$Algorithm, from = paste("AT-BBPSO-MC", c(rep(1, 4), rep(3, 4), rep(5, 4), rep("Inf", 4)), rep(c(0.1, 0.3, 0.5, 0.7), 4), sep = "-"),
          to = paste("$df = ", c(rep(paste(1, ",\\enspace", sep=""), 4), rep(paste(3, ",\\enspace", sep=""), 4), rep(paste(5, ",\\enspace", sep=""), 4), rep("\\infty,", 4)), "$ $R^* =", rep(c(0.1, 0.3, 0.5, 0.7), 4), "$", sep = ""))
tabout$Algorithm <- mapvalues(tabout$Algorithm, from = paste("AT-BBPSOxp-MC", c(rep(1, 4), rep(3, 4), rep(5, 4), rep("Inf", 4)), rep(c(0.1, 0.3, 0.5, 0.7), 4), sep = "-"),
                              to = paste("$df = ", c(rep(paste(1, ",\\enspace", sep=""), 4), rep(paste(3, ",\\enspace", sep=""), 4), rep(paste(5, ",\\enspace", sep=""), 4), rep("\\infty,", 4)), "$ $R^* =", rep(c(0.1, 0.3, 0.5, 0.7), 4), "$", sep = ""))
tabout$Algorithm <- mapvalues(tabout$Algorithm,
                              from = paste("DI-PSO", c(rep(50, 3), rep(100, 3), rep(200, 3)), rep(c(1, 2, 4), 3), sep = "-"),
                              to = paste("$\\alpha = ", c(rep(paste(50, ",\\enspace", sep=""), 3), rep("100,", 3), rep("200,", 3)), "$ $\\beta =", rep(c(1, 2, 4), 3), "$", sep = ""))
tabout$Algorithm <- mapvalues(tabout$Algorithm,
                              from = paste("AT-PSO", c(0.1, 0.3, 0.5, 0.7), sep="-"),
                              to = paste("$R^* = ", c(0.1, 0.3, 0.5, 0.7), "$", sep=""))
tabout$Algorithm <- mapvalues(tabout$Algorithm,
                              from = c("PSO", "BBPSO-MC", "BBPSOxp-MC"),
                              to = c("\\multicolumn{1}{l|}{PSO}", "\\multicolumn{1}{l|}{BBPSO-MC}",
                                     "\\multicolumn{1}{l|}{BBPSOxp-MC}"))
print(xtable(tabout), include.rownames=FALSE, sanitize.text.function=identity)
}

library(ggplot2)

atplotout <- subset(psoout, (obj %in% c(1, 8) & nbhd == "ring-1" & type == "AT-PSO-0.5" &
                             time > 0 & rep == 45))
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











#################### old under hear ########################




sumlastmod <- sumlast
sumlastmod <- sumlastmod[!(sumlastmod$Algorithm == "PSO" & sumlastmod$Nbhd == "Global" & sumlastmod$obj==1),]
sumlastmod$q05[(sumlastmod$q05 < sumlastmod$q25[(sumlastmod$obj == 1 & sumlastmod$Algorithm == "AT-BBPSOxp-MC-2" & sumlastmod$Nbhd == "Ring-3")] & sumlastmod$obj == 1)] <- sumlastmod$q25[(sumlastmod$obj == 1 & sumlastmod$Algorithm == "AT-BBPSOxp-MC-2" & sumlastmod$Nbhd == "Ring-3")]
sumlastmod <- sumlastmod[!(sumlastmod$Algorithm == "PSO" & sumlastmod$obj == 3) & !(sumlastmod$Algorithm ==  "AT-BBPSO-MC-1" & sumlastmod$obj == 3),]




sumlastmod$q05[(sumlastmod$Algorithm == "BBPSOxp-MC" & sumlastmod$Nbhd == "Global" & sumlastmod$obj==1)] <- sumlastmod$q25[(sumlastmod$Algorithm == "BBPSOxp-MC" & sumlastmod$Nbhd == "Global" & sumlastmod$obj==1)]
sumlastmod$q05[(sumlastmod$Algorithm == "AT-BBPSOxp-MC-2" & sumlastmod$Nbhd == "Global" & sumlastmod$obj==1)] <- sumlastmod$q25[(sumlastmod$Algorithm == "BBPSOxp-MC" & sumlastmod$Nbhd == "Global" & sumlastmod$obj==1)]
sumlastmod$q05[(sumlastmod$Algorithm == "AT-BBPSOxp-MC-2" & sumlastmod$Nbhd == "Global" & sumlastmod$obj==2)] <- sumlastmod$q25[(sumlastmod$Algorithm == "AT-BBPSOxp-MC-2" & sumlastmod$Nbhd == "Global" & sumlastmod$obj==2)]
sumlastmod$q05[(sumlastmod$Algorithm == "BBPSOxp-MC" & sumlastmod$Nbhd == "Global" & sumlastmod$obj==2)] <- sumlastmod$q25[(sumlastmod$Algorithm == "AT-BBPSOxp-MC-2" & sumlastmod$Nbhd == "Global" & sumlastmod$obj==2)]
sumlastmod$q05[(sumlastmod$Algorithm == "PSO" & sumlastmod$Nbhd == "Global" & sumlastmod$obj==2)] <- sumlastmod$q25[(sumlastmod$Algorithm == "AT-BBPSOxp-MC-2" & sumlastmod$Nbhd == "Global" & sumlastmod$obj==2)]
sumlastmod$q05[(sumlastmod$Algorithm == "AT-BBPSO-MC-1" & sumlastmod$Nbhd == "Global" & sumlastmod$obj==2)] <- sumlastmod$q25[(sumlastmod$Algorithm == "AT-BBPSOxp-MC-2" & sumlastmod$Nbhd == "Global" & sumlastmod$obj==2)]
sumlastmod$q05[(sumlastmod$Algorithm == "AT-BBPSO-MC-1" & sumlastmod$Nbhd == "Global" & sumlastmod$obj==3)] <- sumlastmod$q25[(sumlastmod$Algorithm == "AT-BBPSO-MC-1" & sumlastmod$Nbhd == "Global" & sumlastmod$obj==3)]
sumlastmod$q05[(sumlastmod$Algorithm == "PSO" & sumlastmod$Nbhd == "Global" & sumlastmod$obj==4)] <- sumlastmod$q25[(sumlastmod$Algorithm == "PSO" & sumlastmod$Nbhd == "Global" & sumlastmod$obj==4)]
sumlastmod$q95[(sumlastmod$Algorithm == "AT-BBPSO-MC-1" & sumlastmod$Nbhd == "Global" & sumlastmod$obj==6)] <- -18
sumlastmod$q95[(sumlastmod$Algorithm == "AT-BBPSOxp-MC-1" & sumlastmod$Nbhd == "Global" & sumlastmod$obj==6)] <- -18
sumlastmod$q05[(sumlastmod$Algorithm == "PSO" & sumlastmod$Nbhd == "Ring-1" & sumlastmod$obj==1)] <- -100
sumlastmod$q05[(sumlastmod$Algorithm == "PSO" & sumlastmod$Nbhd == "Ring-1" & sumlastmod$obj==3)] <- sumlastmod$q25[(sumlastmod$Algorithm == "PSO" & sumlastmod$Nbhd == "Ring-1" & sumlastmod$obj==3)]
sumlastmod$q05[(sumlastmod$Algorithm == "PSO" & sumlastmod$Nbhd == "Ring-3" & sumlastmod$obj==3)] <- -21000
sumlastmod$q05[(sumlastmod$Algorithm == "PSO" & sumlastmod$Nbhd == "Ring-3" & sumlastmod$obj==1)] <- -400
sumlastmod$q05[(sumlastmod$Algorithm == "BBPSOxp-MC" & sumlastmod$Nbhd == "Ring-3" & sumlastmod$obj==1)] <- -400
sumlastmod$q05[(sumlastmod$Algorithm == "AT-BBPSOxp-MC-2" & sumlastmod$Nbhd == "Ring-3" & sumlastmod$obj==1)] <- -400

library(ggplot2)

maxplot1 <- ggplot(subset(sumlastmod, obj == 1)) +
  geom_boxplot(aes(x=Algorithm, lower = q25, upper = q75, middle = Median, ymin = q05,
                   ymax = q95), stat="identity") +
  facet_wrap(Nbhd~obj, ncol = 1,
             labeller = label_bquote(paste(.(Nbhd), ", ", f[.(obj)]))) +
  theme(axis.text = element_text(angle = 45, hjust = 1),
                strip.text = element_text(size = 11)) +
  scale_x_discrete(drop=FALSE)
maxplot2 <- ggplot(subset(sumlastmod, obj == 2)) +
  geom_boxplot(aes(x=Algorithm, lower = q25, upper = q75, middle = Median, ymin = q05,
                   ymax = q95), stat="identity") +
  facet_wrap(Nbhd~obj, ncol = 1,
             labeller = label_bquote(paste(.(Nbhd), ", ", f[.(obj)]))) +
  theme(axis.text = element_text(angle = 45, hjust = 1),
                strip.text = element_text(size = 11)) +
  scale_x_discrete(drop=FALSE)
maxplot3 <- ggplot(subset(sumlastmod, obj == 3)) +
  geom_boxplot(aes(x=Algorithm, lower = q25, upper = q75, middle = Median, ymin = q05,
                   ymax = q95), stat="identity") +
  facet_wrap(Nbhd~obj, ncol = 1,
             labeller = label_bquote(paste(.(Nbhd), ", ", f[.(obj)]))) +
  theme(axis.text = element_text(angle = 45, hjust = 1),
                strip.text = element_text(size = 11)) +
  scale_x_discrete(drop=FALSE)
maxplot4 <- ggplot(subset(sumlastmod, obj == 4)) +
  geom_boxplot(aes(x=Algorithm, lower = q25, upper = q75, middle = Median, ymin = q05,
                   ymax = q95), stat="identity") +
  facet_wrap(Nbhd~obj, ncol = 1,
             labeller = label_bquote(paste(.(Nbhd), ", ", f[.(obj)]))) +
  theme(axis.text = element_text(angle = 45, hjust = 1),
                strip.text = element_text(size = 11)) +
  scale_x_discrete(drop=FALSE)
maxplot5 <- ggplot(subset(sumlastmod, obj == 5)) +
  geom_boxplot(aes(x=Algorithm, lower = q25, upper = q75, middle = Median, ymin = q05,
                   ymax = q95), stat="identity") +
  facet_wrap(Nbhd~obj, ncol = 1,
             labeller = label_bquote(paste(.(Nbhd), ", ", f[.(obj)]))) +
  theme(axis.text = element_text(angle = 45, hjust = 1),
                strip.text = element_text(size = 11)) +
  scale_x_discrete(drop=FALSE)
maxplot6 <- ggplot(subset(sumlastmod, obj == 6)) +
  geom_boxplot(aes(x=Algorithm, lower = q25, upper = q75, middle = Median, ymin = q05,
                   ymax = q95), stat="identity") +
  facet_wrap(Nbhd~obj, ncol = 1,
             labeller = label_bquote(paste(.(Nbhd), ", ", f[.(obj)]))) +
  theme(axis.text = element_text(angle = 45, hjust = 1),
                strip.text = element_text(size = 11)) +
  scale_x_discrete(drop=FALSE)


options(scipen = -2)
ht <- 6
wd <- 3
ggsave("maxplot11.png", maxplot1, width = wd, height = ht)
ggsave("maxplot12.png", maxplot2, width = wd, height = ht)
ggsave("maxplot13.png", maxplot3, width = wd, height = ht)
ggsave("maxplot14.png", maxplot4, width = wd, height = ht)
ggsave("maxplot15.png", maxplot5, width = wd, height = ht)
ggsave("maxplot16.png", maxplot6, width = wd, height = ht)



options(scipen = -1)
ggsave("argmaxplot21.png", argmaxplot1, width = 12, height = 7)
ggsave("argmaxplot22.png", argmaxplot2, width = 12, height = 7)
