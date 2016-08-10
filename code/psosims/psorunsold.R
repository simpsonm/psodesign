library(ggplot2)
library(gridExtra)
library(mnormt)
source('psofunold.R')
library(plyr)
source('testfuns.R')

nswarm <- 20
np <- 30
niter <- 500
nrep <- 1000
inertias <- 0.7298
cognitives <- 1.496
socials <- 1.496
sig0 <- 1
rates <- c(.2, .4)
dfs <- c(1, Inf)
sigs <- 1
eta <- 2
maxfail <- 5

nbhd <- list()
nbhd[[1]] <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)})
nbhd[[2]] <- matrix(1:nswarm, ncol=nswarm, nrow=nswarm)
time <- 0:niter

normvec <- function(x){
  sqrt(mean(x^2))
}

set.seed(342)

initfun <- function(n, nswarm, np){
  inits <- switch(n,
         matrix(runif(nswarm*np, 50, 100), ncol = nswarm),
         matrix(runif(nswarm*np, 50, 100), ncol = nswarm),
         matrix(runif(nswarm*np, 50, 100), ncol = nswarm),
         matrix(runif(nswarm*np, 15, 30), ncol = nswarm),
         matrix(runif(nswarm*np, 2.56, 5.12), ncol = nswarm),
         matrix(runif(nswarm*np, -500, -250), ncol = nswarm),
         matrix(runif(nswarm*np, 300, 600), ncol = nswarm),
         matrix(runif(nswarm*np, 16, 32), ncol = nswarm))
  return(inits)
}



## PSO
{
psoout <- NULL
bbpsoout <- NULL
bbpsojout <- NULL
atbbpsoout <- NULL
atbbpsojout <- NULL
for(n in c(1:8)[-c(3,6)]){
  for(m in 1:2){
    for(l in 1:nrep){
      print(c(n,m,l))
      inits <- initfun(n, nswarm, np)
      temp <- pso(niter, nswarm, inertias, cognitives, socials, inits, nbhd[[m]], fwrap, opt=n)
      psoout <- rbind(psoout, data.frame(obj = n, logpost = temp[["maxes"]], argnorm = apply(temp[["argmaxes"]], 2, normvec), time = time, type = "pso", nbhd = m, rep = l))
      temp <- bbpso(niter, nswarm, sigs, eta, 0, Inf, inits, nbhd[[m]], fwrap, Inf, FALSE, opt = n)
      bbpsoout <- rbind(bbpsoout, data.frame(obj = n, logpost= temp[["maxes"]], argnorm = apply(temp[["argmaxes"]], 2, normvec), time = time, type = "bbpso", nbhd = m, rep = l))
      temp <- bbpso(niter, nswarm, sigs, eta, 0, maxfail, inits, nbhd[[m]], fwrap, Inf, FALSE, opt = n)
      bbpsojout <- rbind(bbpsojout, data.frame(obj = n, logpost= temp[["maxes"]], argnorm = apply(temp[["argmaxes"]], 2, normvec), time = time, type = "bbpsoj", nbhd = m, rep = l))
      for(i in 1:length(dfs)){
        temp <- bbpso(niter, nswarm, sig0, 1, rates[i], Inf, inits, nbhd[[m]], fwrap, dfs[i], TRUE, opt = n)
        atbbpsoout <- rbind(atbbpsoout, data.frame(obj = n, logpost= temp[["maxes"]], argnorm = apply(temp[["argmaxes"]], 2, normvec), time = time, type = "atbbpso", nbhd = m, df = dfs[i], R = rates[i], rep = l))
        temp <- bbpso(niter, nswarm, sig0, 1, rates[i], Inf, inits, nbhd[[m]], fwrap, dfs[i], TRUE, opt = n)
        atbbpsojout <- rbind(atbbpsojout, data.frame(obj = n, logpost= temp[["maxes"]], argnorm = apply(temp[["argmaxes"]], 2, normvec), time = time, type = "atbbpsoj", nbhd = m, df = dfs[i], R = rates[i], rep = l))
      }
    }
  }
}
}


save(psoout, file = "psoout.Rdata")
save(bbpsoout, file = "bbpsoout.Rdata")
save(bbpsojout, file = "bbpsojout.Rdata")
save(atbbpsoout, file = "atbbpsoout.Rdata")
save(atbbpsojout, file = "atbbpsojout.Rdata")

load("psoout.Rdata")
load("bbpsoout.Rdata")
load("bbpsojout.Rdata")
load("atbbpsoout.Rdata")
load("atbbpsojout.Rdata")

psoout$nbhd[psoout$nbhd == 1] <- "Ring-1"
psoout$nbhd[psoout$nbhd == 2] <- "Gbest"
psoout$type <- "PSO"
psosum <- ddply(psoout, .(obj, type, nbhd, time), summarise, mean=mean(logpost), median=median(logpost), sd=sd(logpost), min=min(logpost), q025=quantile(logpost, 0.025), q25=quantile(logpost, 0.25), q75=quantile(logpost, 0.75), q975=quantile(logpost, 0.975), max=max(logpost))
psoargsum <- ddply(psoout, .(obj, type, nbhd, time), summarise, mean=mean(argnorm), median=median(argnorm), sd=sd(argnorm), min=min(argnorm), q025=quantile(argnorm, 0.025), q25=quantile(argnorm, 0.25), q75=quantile(argnorm, 0.75), q975=quantile(argnorm, 0.975), max=max(argnorm))
rm(psoout)
bbpsoout$nbhd[bbpsoout$nbhd == 1] <- "Ring-1"
bbpsoout$nbhd[bbpsoout$nbhd == 2] <- "Gbest"
bbpsoout$type <- "BBPSO"
bbpsosum <- ddply(bbpsoout, .(obj, type, nbhd, time), summarise, mean=mean(logpost), median=median(logpost), sd=sd(logpost), min=min(logpost), q025=quantile(logpost, 0.025), q25=quantile(logpost, 0.25), q75=quantile(logpost, 0.75), q975=quantile(logpost, 0.975), max=max(logpost))
bbpsoargsum <- ddply(bbpsoout, .(obj, type, nbhd, time), summarise, mean=mean(argnorm), median=median(argnorm), sd=sd(argnorm), min=min(argnorm), q025=quantile(argnorm, 0.025), q25=quantile(argnorm, 0.25), q75=quantile(argnorm, 0.75), q975=quantile(argnorm, 0.975), max=max(argnorm))
rm(bbpsoout)
bbpsojout$nbhd[bbpsojout$nbhd == 1] <- "Ring-1"
bbpsojout$nbhd[bbpsojout$nbhd == 2] <- "Gbest"
bbpsojout$type <- "BBPSO-CJ"
bbpsojsum <- ddply(bbpsojout, .(obj, type, nbhd, time), summarise, mean=mean(logpost), median=median(logpost), sd=sd(logpost), min=min(logpost), q025=quantile(logpost, 0.025), q25=quantile(logpost, 0.25), q75=quantile(logpost, 0.75), q975=quantile(logpost, 0.975), max=max(logpost))
bbpsojargsum <- ddply(bbpsojout, .(obj, type, nbhd, time), summarise, mean=mean(argnorm), median=median(argnorm), sd=sd(argnorm), min=min(argnorm), q025=quantile(argnorm, 0.025), q25=quantile(argnorm, 0.25), q75=quantile(argnorm, 0.75), q975=quantile(argnorm, 0.975), max=max(argnorm))
rm(bbpsojout)
atbbpsoout$nbhd[atbbpsoout$nbhd == 1] <- "Ring-1"
atbbpsoout$nbhd[atbbpsoout$nbhd == 2] <- "Gbest"
atbbpsoout$type <- as.character(atbbpsoout$type)
atbbpsoout$type[atbbpsoout$df == Inf] <- "AT-BBPSO-1"
atbbpsoout$type[atbbpsoout$df == 1] <- "AT-BBPSO-2"
atbbpsosum <- ddply(atbbpsoout, .(obj, type, nbhd, time), summarise, mean=mean(logpost), median=median(logpost), sd=sd(logpost), min=min(logpost), q025=quantile(logpost, 0.025), q25=quantile(logpost, 0.25), q75=quantile(logpost, 0.75), q975=quantile(logpost, 0.975), max=max(logpost))
atbbpsoargsum <- ddply(atbbpsoout, .(obj, type, nbhd, time), summarise, mean=mean(argnorm), median=median(argnorm), sd=sd(argnorm), min=min(argnorm), q025=quantile(argnorm, 0.025), q25=quantile(argnorm, 0.25), q75=quantile(argnorm, 0.75), q975=quantile(argnorm, 0.975), max=max(argnorm))
rm(atbbpsoout)
atbbpsojout$nbhd[atbbpsojout$nbhd == 1] <- "Ring-1"
atbbpsojout$nbhd[atbbpsojout$nbhd == 2] <- "Gbest"
atbbpsojout$type <- as.character(atbbpsojout$type)
atbbpsojout$type[atbbpsojout$df == Inf] <- "AT-BBPSO-CJ-1"
atbbpsojout$type[atbbpsojout$df == 1] <- "AT-BBPSO-CJ-2"
atbbpsojsum <- ddply(atbbpsojout, .(obj, type, nbhd, time), summarise, mean=mean(logpost), median=median(logpost), sd=sd(logpost), min=min(logpost), q025=quantile(logpost, 0.025), q25=quantile(logpost, 0.25), q75=quantile(logpost, 0.75), q975=quantile(logpost, 0.975), max=max(logpost))
atbbpsojargsum <- ddply(atbbpsojout, .(obj, type, nbhd, time), summarise, mean=mean(argnorm), median=median(argnorm), sd=sd(argnorm), min=min(argnorm), q025=quantile(argnorm, 0.025), q25=quantile(argnorm, 0.25), q75=quantile(argnorm, 0.75), q975=quantile(argnorm, 0.975), max=max(argnorm))
rm(atbbpsojout)
psosums <- rbind(psosum, bbpsosum, bbpsojsum, atbbpsosum, atbbpsojsum)
psoargsums <- rbind(psoargsum, bbpsoargsum, bbpsojargsum, atbbpsoargsum, atbbpsojargsum)
save(psosums, file = "psosums.Rdata")
save(psoargsums, file = "psoargsums.Rdata")

load("psosums.Rdata")
load("psoargsums.Rdata")

ggplot(psosums, aes(time, mean, color = type, linetype = nbhd)) + geom_line(size = I(1), subset=.(time > 400)) + facet_wrap(~obj, scale="free")

sumlast <- subset(psosums, time == niter)
sumlast$type <- factor(sumlast$type, levels = c("PSO", "BBPSO", "BBPSO-CJ", "AT-BBPSO-1", "AT-BBPSO-2", "AT-BBPSO-CJ-1", "AT-BBPSO-CJ-2"))
sumlast <- sumlast[order(sumlast$type),]
colnames(sumlast)[c(2,3,5,6,7)] <- c("Algorithm", "Nbhd", "Mean", "Median", "SD")

argsumlast <- subset(psoargsums, time == niter)
argsumlast$type <- factor(argsumlast$type, levels = c("PSO", "BBPSO", "BBPSO-CJ", "AT-BBPSO-1", "AT-BBPSO-2", "AT-BBPSO-CJ-1", "AT-BBPSO-CJ-2"))
argsumlast <- argsumlast[order(argsumlast$type),]
colnames(argsumlast)[c(2,3,5,6,7)] <- c("Algorithm", "Nbhd", "Mean", "Median", "SD")


head(sumlast[order(sumlast$Median, decreasing = TRUE),c(2,3,5,6,7,8,9)], 20)

library(xtable)


for(onum in c(1:8)[-c(3,6)]){
  print(xtable(subset(sumlast, obj==onum)[,c(2,3,5,7)]), include.rownames=FALSE)
}


out1 <- cbind(subset(sumlast, obj==1)[,c(2,3,5,7)], subset(sumlast, obj==2)[,c(5,7)], cbind(subset(sumlast, obj==4)[,c(5,7)]))
out2 <- cbind(subset(sumlast, obj==5)[,c(2,3,5,7)], subset(sumlast, obj==7)[,c(5,7)], cbind(subset(sumlast, obj==8)[,c(5,7)]))
colnames(out1)[5:6] <- c("Mean", "SD")

digits1 <- matrix(c(rep(rep(0, 14),3), rep(c(0,1,5,2,5,1,7,5,4,2,7,5,4,2),2), rep(c(-2,-2,-2,-3,-2,-3,-2,-2,-2,-3,-2,-2,-2,-3),2), rep(c(-5,-2,-1,-3,-1,-3,-1,-3,-1,-1,-1,-3,-1,-1),2)), ncol = 9)
digits1[digits1>4] <- -2
digits1[digits1< 0] <- -2
digits2 <- matrix(c(rep(rep(0, 14),3), rep(c(1,rep(2,13)),2), rep(c(1,2,4,3,4,3,4,3,3,3,4,3,3,3),2), rep(c(3,3,2,2,2,3,3,3,3,4,3,3,3,4),2)), ncol = 9)
digits2[digits2>4] <- -2
digits2[digits2< 0] <- -2

print(xtable(out1, digits=digits1), include.rownames=FALSE)
print(xtable(out2, digits=digits2), include.rownames=FALSE)



out3 <- cbind(subset(argsumlast, obj==1)[,c(2,3,5,7)], subset(argsumlast, obj==2)[,c(5,7)], cbind(subset(argsumlast, obj==4)[,c(5,7)]))
out4 <- cbind(subset(argsumlast, obj==5)[,c(2,3,5,7)], subset(argsumlast, obj==7)[,c(5,7)], cbind(subset(argsumlast, obj==8)[,c(5,7)]))
colnames(out3)[5:6] <- c("Mean", "SD")

digits3 <- matrix(c(rep(rep(0, 14),3), rep(c(2,3,4,3,4,3,5,5,4,3,6,5,4,3),2), rep(rep(1,14),2), rep(c(2,2,rep(1,12)),2)), ncol = 9)
digits3[digits3>4] <- -2
digits3[digits3< 0] <- -2
digits4 <- matrix(c(rep(rep(0, 14),3), rep(c(3,3,4,4,3,4,3,3,3,4,3,3,3,4),2), rep(c(1,2,3,2,3,2,2,2,3,2,2,2,3,2),2), rep(c(2,2,1,-2,-1,-3,1,-2,-6,-5,1,-2,-5,-8),2)), ncol = 9)
digits4[digits4>4] <- -2
digits4[digits4< 0] <- -2

print(xtable(out3, digits=digits3), include.rownames=FALSE)
print(xtable(out4, digits=digits4), include.rownames=FALSE)


sumlastmod <- sumlast
sumlastmod <- sumlastmod[!(sumlastmod$Algorithm == "PSO" & sumlastmod$Nbhd == "Gbest" & sumlastmod$obj==1),]
sumlastmod <- sumlastmod[!(sumlastmod$Algorithm == "PSO" & sumlastmod$Nbhd == "Gbest" & sumlastmod$obj==4),]
sumlastmod <- sumlastmod[!(sumlastmod$Algorithm == "PSO" & sumlastmod$Nbhd == "Gbest" & sumlastmod$obj==5),]
sumlastmod <- sumlastmod[!(sumlastmod$Algorithm == "PSO" & sumlastmod$Nbhd == "Gbest" & sumlastmod$obj==7),]
sumlastmod <- sumlastmod[!(sumlastmod$Algorithm == "BBPSO" & sumlastmod$Nbhd == "Ring-1" & sumlastmod$obj==4),]
sumlastmod <- sumlastmod[!(sumlastmod$Algorithm == "BBPSO-CJ" & sumlastmod$Nbhd == "Ring-1" & sumlastmod$obj==4),]

sumlastmod$q025[(sumlastmod$Algorithm == "AT-BBPSO-2" & sumlastmod$Nbhd == "Gbest" & sumlastmod$obj==1)] <- sumlastmod$q25[(sumlastmod$Algorithm == "AT-BBPSO-2" & sumlastmod$Nbhd == "Gbest" & sumlastmod$obj==1)] 
sumlastmod$q025[(sumlastmod$Algorithm == "AT-BBPSO-CJ-2" & sumlastmod$Nbhd == "Gbest" & sumlastmod$obj==1)] <- sumlastmod$q25[(sumlastmod$Algorithm == "AT-BBPSO-2" & sumlastmod$Nbhd == "Gbest" & sumlastmod$obj==1)]
sumlastmod$q025[sumlastmod$Nbhd == "Gbest" & sumlastmod$obj == 4] <- -5000
sumlastmod$q975[sumlastmod$Nbhd == "Gbest" & sumlastmod$obj == 8 & sumlastmod$Algorithm %in% c("BBPSO", "BBPSO-CJ")] <- -17
sumlastmod$q025[(sumlastmod$Algorithm == "BBPSO-CJ" & sumlastmod$Nbhd == "Ring-1" & sumlastmod$obj==1)] <- sumlastmod$q025[(sumlastmod$Algorithm == "BBPSO" & sumlastmod$Nbhd == "Ring-1" & sumlastmod$obj==1)]
sumlastmod$q025[(sumlastmod$Algorithm == "PSO" & sumlastmod$Nbhd == "Ring-1" & sumlastmod$obj==4)] <- sumlastmod$q025[(sumlastmod$Algorithm == "AT-BBPSO-CJ-1" & sumlastmod$Nbhd == "Ring-1" & sumlastmod$obj==1)]

sumlastmod$Nbhd[sumlastmod$Nbhd == "Gbest"] <- "Global"
sumlastmod$obj[sumlastmod$obj > 2] <- sumlastmod$obj[sumlastmod$obj > 2] - 1
sumlastmod$obj[sumlastmod$obj > 5] <- sumlastmod$obj[sumlastmod$obj > 5] - 1

maxplot <- ggplot(sumlastmod) +
  geom_boxplot(aes(x=Algorithm, lower = q25, upper = q75, middle = Median, ymin = q025,
                   ymax = q975), stat="identity") +
  facet_wrap(Nbhd~obj, scales="free_y", ncol = 6,
             labeller = label_bquote(paste(.(Nbhd), ", ", f[.(obj)]))) +
  theme(axis.text = element_text(angle = 45, hjust = 1),
                strip.text = element_text(size = 11)) +
  ggtitle("Maximum value found")

argsumlastmod <- argsumlast
argsumlastmod <- argsumlastmod[!(argsumlastmod$Algorithm == "PSO" & argsumlastmod$Nbhd == "Gbest" & argsumlastmod$obj==1),]
argsumlastmod <- argsumlastmod[!(argsumlastmod$Algorithm == "PSO" & argsumlastmod$Nbhd == "Gbest" & argsumlastmod$obj==7),]
argsumlastmod <- argsumlastmod[!(argsumlastmod$Algorithm == "AT-BBPSO-2" & argsumlastmod$Nbhd == "Gbest" & argsumlastmod$obj==8),]
argsumlastmod <- argsumlastmod[!(argsumlastmod$Algorithm == "AT-BBPSO-CJ-2" & argsumlastmod$Nbhd == "Gbest" & argsumlastmod$obj==8),]
argsumlastmod <- argsumlastmod[!(argsumlastmod$Algorithm == "AT-BBPSO-2" & argsumlastmod$Nbhd == "Ring-1" & argsumlastmod$obj==8),]
argsumlastmod <- argsumlastmod[!(argsumlastmod$Algorithm == "AT-BBPSO-CJ-2" & argsumlastmod$Nbhd == "Ring-1" & argsumlastmod$obj==8),]

argsumlastmod$q975[(argsumlastmod$Algorithm == "BBPSO-CJ" & argsumlastmod$Nbhd == "Gbest" & argsumlastmod$obj==8)] <- argsumlastmod$q975[(argsumlastmod$Algorithm == "BBPSO" & argsumlastmod$Nbhd == "Gbest" & argsumlastmod$obj==8)]
argsumlastmod$q975[(argsumlastmod$Algorithm == "BBPSO-CJ" & argsumlastmod$Nbhd == "Ring-1" & argsumlastmod$obj==8)] <- argsumlastmod$q75[(argsumlastmod$Algorithm == "BBPSO-CJ" & argsumlastmod$Nbhd == "Ring-1" & argsumlastmod$obj==8)]
argsumlastmod$q975[(argsumlastmod$Algorithm == "BBPSO" & argsumlastmod$Nbhd == "Ring-1" & argsumlastmod$obj==8)] <- argsumlastmod$q75[(argsumlastmod$Algorithm == "BBPSO-CJ" & argsumlastmod$Nbhd == "Ring-1" & argsumlastmod$obj==8)]
argsumlastmod$q975[(argsumlastmod$Algorithm == "AT-BBPSO-CJ-1" & argsumlastmod$Nbhd == "Ring-1" & argsumlastmod$obj==8)] <- argsumlastmod$q75[(argsumlastmod$Algorithm == "BBPSO-CJ" & argsumlastmod$Nbhd == "Ring-1" & argsumlastmod$obj==8)]
argsumlastmod$q975[(argsumlastmod$Algorithm == "AT-BBPSO-1" & argsumlastmod$Nbhd == "Ring-1" & argsumlastmod$obj==8)] <- argsumlastmod$q75[(argsumlastmod$Algorithm == "BBPSO-CJ" & argsumlastmod$Nbhd == "Ring-1" & argsumlastmod$obj==8)]

argsumlastmod$Nbhd[argsumlastmod$Nbhd == "Gbest"] <- "Global"
argsumlastmod$obj[argsumlastmod$obj > 2] <- argsumlastmod$obj[argsumlastmod$obj > 2] - 1
argsumlastmod$obj[argsumlastmod$obj > 5] <- argsumlastmod$obj[argsumlastmod$obj > 5] - 1
##argsumlastmod$obj <- paste("f[", argsumlastmod$obj, "]", sep="")


argmaxplot <- ggplot(argsumlastmod) +
  geom_boxplot(aes(x=Algorithm, lower = q25, upper = q75, middle = Median, ymin = q025,
                   ymax = q975), stat="identity") +
  facet_wrap(Nbhd~obj, scales="free_y", ncol = 6,
             labeller = label_bquote(paste(.(Nbhd), ", ", f[.(obj)]))) +
  theme(axis.text = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 11)) +
  ggtitle("Euclidean distance from the true maximum")

options(scipen = -2)
ggsave("maxplot.png", maxplot, width = 12, height = 7)
options(scipen = -1)
ggsave("argmaxplot.png", argmaxplot, width = 12, height = 7)
