library(plyr)
library(dplyr)
library(ggplot2)
library(xtable)
library(gridExtra)
library(ggmap)
library(reshape2)
library(geoR)
library(shapefiles)
library(maptools)
library(rgdal)

load("chrispsoouts.RData")
load("noelpsoouts.RData")
load("noelpsooutsextra.RData")
load("homepsoouts.RData")

valueout <- homepsoouts[[1]][[1]]
parout <- homepsoouts[[1]][[2]]
for(i in 2:length(homepsoouts)){
  valueout <- rbind(valueout, homepsoouts[[i]][[1]])
  parout <- rbind(parout, homepsoouts[[i]][[2]])
}
for(i in 1:length(chrispsoouts)){
  valueout <- rbind(valueout, chrispsoouts[[i]][[1]])
  parout <- rbind(parout, chrispsoouts[[i]][[2]])
}
for(i in 1:length(noelpsoouts)){
  valueout <- rbind(valueout, noelpsoouts[[i]][[1]])
  parout <- rbind(parout, noelpsoouts[[i]][[2]])
}
for(i in 1:length(noelpsooutsextra)){
  valueout <- rbind(valueout, noelpsooutsextra[[i]][[1]])
  parout <- rbind(parout, noelpsooutsextra[[i]][[2]])
}


save(valueout, file = "valueout.RData")
save(parout, file = "parout.RData")

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



valueclean <- ddply(valueout, .(obj, algid, nbhd, time, logpost, inertias), namesfun)[, c(1, 7, 8, 4, 5, 6)]

names(valueclean) <- c("Obj", "Algorithm", "Nbhd", "Time", "logpost", "inertia")

valueclean$Obj <- mapvalues(valueclean$Obj, c("sig2fuk.max", "sig2fuk.mean"),
                            c("sig2puk.max", "sig2puk.mean"))


algorder <- c("PSO1", "PSO2", "PSO1-CF", "PSO2-CF",
              "DI-PSO1", "DI-PSO2", "DI-PSO1-CF", "DI-PSO2-CF",
              "AT3-PSO1", "AT3-PSO2", "AT3-PSO1-CF", "AT3-PSO2-CF",
              "AT5-PSO1", "AT5-PSO2", "AT5-PSO1-CF", "AT5-PSO2-CF",
              "AT3-BBPSO", "AT3-BBPSOxp", "AT3-BBPSO-CF", "AT3-BBPSOxp-CF",
              "AT5-BBPSO", "AT5-BBPSOxp", "AT5-BBPSO-CF", "AT5-BBPSOxp-CF")
nbhdorder <- c("Global", "SS3", "SS1")

valueclean$Algorithm <- factor(valueclean$Algorithm, levels = algorder)
valueclean$Nbhd <- factor(valueclean$Nbhd, levels = nbhdorder)
valueorder <- arrange(valueclean, Obj, Algorithm, Nbhd, Time)


nlook <- floor(max(valueorder$Time))
valuelast <- filter(valueorder, Time == nlook)
## this prints the tables, however it's missing the first header.
## paste the following line (without the ##) just after the \begin{tabular} line
## \multicolumn{1}{l}{\begin{tabular}{lr} ; & Nbhd: \end{tabular}} & \multicolumn{2}{1}{Global} & \multicolumn{1}{c}{SS3}\\

for(k in unique(valuelast$Obj)){
  ofk <- filter(valuelast, Obj == k & !(Algorithm %in% c("DI-PSO1", "DI-PSO2", "DI-PSO1-CF", "DI-PSO2-CF")))[,-1]
  ## remove large values from table
  bigidx <- which(ofk$logpost > 10000)
  ofk$logpost[bigidx] <- NA
  ## pull out each type of nbhd and combine
  ofk1 <- filter(ofk, Nbhd == "Global")[,-c(2,3,5)]
  ofk2 <- filter(ofk, Nbhd == "SS3")[,-c(2,3,5)]
  ofkout <- cbind(ofk1, ofk2[,-1])
  colnames(ofkout)[3] <- "logpost"
  ##iidd <- c(3:4, 7:8, 11:12, 15:16, 19:20, 23:24)
  ##ofprint <- ofkout[-iidd,-c(3,5,7)]
  print(xtable(ofkout, digits = c(0,0,rep(c(2),2)), align = "ll|r|r",
               caption = paste("Simulation results for ", k, ". See text for description", sep = ""),
               label = paste("tab:psosim", k, sep="")),
        include.rownames=FALSE, sanitize.text.function=identity, hline.after=c(-1,seq(0, 20, 4)),
        size = "\\scriptsize")
}






## map initialization
houston <- read.csv("houstonout.csv")
coshape <- readOGR(dsn="shape", layer="05000")
names(coshape)[1] <- "ID"
cogeom <- fortify(coshape)
ids <- sapply(coshape@polygons, function(x){x@ID})
housgeom <- filter(cogeom, id == ids[117])

d <- data.frame( x = housgeom$long, y = housgeom$lat)
coordinates(d) <- c("x", "y")
proj4string(d) <- proj4string(coshape)
CRS.new <- CRS("+init=epsg:4326") # WGS 84, i.e. world
d2 <- spTransform(d, CRS.new)
housgeom$longitude <- d2@coords[,1]
housgeom$latitude <- d2@coords[,2]

basemap <- get_map(location = "houston", zoom = 9, maptype = 'terrain')
p <- ggmap(basemap) + geom_point(aes(Longitude, Latitude), data = houston, size = I(3), alpha=0.6)

## prediction locations map
load("datlist.RData")

p + geom_polygon(aes(longitude, latitude, group = id), data = housgeom, fill = NA, colour = "black") +
  geom_point(aes(u, v), data = data.frame(u = datlist$tt[,1]/(pi/180*6371),
                                          v = datlist$tt[,2]/(pi/180*6371)),
             colour = "blue", alpha = 0.6, shape="+", size = I(3))


### now use fits too
meandat <- filter(parout, obj == "sig2fuk.mean" & algid == "PSO-2-notCF-AT2" & nbhd > 3)
meanout <- data.frame(matrix(meandat$par, ncol = 2)/(pi/180*6371))
colnames(meanout) <- c("x", "y")

p.mean <- p + geom_polygon(aes(longitude, latitude, group = id), data = housgeom, fill = NA, colour = "black") +
  geom_point(aes(x, y), data = meanout, size = I(3), alpha = 0.6, colour = "red") +
  xlab("longitude") + ylab("latitude")

maxdat <- filter(parout, obj == "sig2fuk.max" & algid == "PSO-1-notCF-AT1" & nbhd > 3)
maxout <- data.frame(matrix(maxdat$par, ncol = 2)/(pi/180*6371))
colnames(maxout) <- c("x", "y")

p.max <- p + geom_polygon(aes(longitude, latitude, group = id), data = housgeom, fill = NA, colour = "black") +
  geom_point(aes(x, y), data = maxout, size = I(3), alpha = 0.6, colour = "red") +
  xlab("longitude") + ylab("latitude")

w <- 7
h <- 7
ggsave(filename = "sig2pukmean.png", plot = p.mean, width = w, height = h)
ggsave(filename = "sig2pukmax.png", plot = p.max, width = w, height = h)







objname <- "sig2puk.mean"
q1 <- qplot(Time, logpost, color = Algorithm, geom = "line", facets = Nbhd~., size = I(1),
            data = filter(valueorder, Obj == objname & Nbhd %in% c("Global", "SS3") &
                                       Algorithm %in% c("PSO1", "PSO2", "AT3-PSO1", "AT3-PSO2")))
q2 <- qplot(Time, logpost, color = Algorithm, geom = "line", facets = Nbhd~., size = I(1),
            data = filter(valueorder, Obj == objname & Nbhd %in% c("Global", "SS3") &
                                       Algorithm %in% c("AT3-BBPSO", "AT3-BBPSOxp", "AT5-BBPSO", "AT5-BBPSOxp")))
grid.arrange(q1, q2, ncol = 2)

objname <- "sig2puk.max"
qplot(Time, logpost, color = Algorithm, geom = "line", facets = Nbhd~., size = I(1),
      data = filter(valueorder, Obj ==objname & Nbhd %in% c("Global", "SS3") &
      Algorithm %in% c("PSO1-CF", "PSO2-CF", "AT3-PSO1-CF", "AT3-PSO2-CF",
                       "AT5-PSO1-CF", "AT5-PSO2-CF",
                       "AT3-BBPSO-CF", "AT3-BBPSOxp-CF",
                       "AT5-BBPSO-CF", "AT5-BBPSOxp-CF")))

objname <- "sig2puk.mean"
qplot(Time, logpost, color = Algorithm, geom = "line", facets = Nbhd~., size = I(1),
      data = filter(valueorder, Obj ==objname & Nbhd %in% c("Global", "SS3") &
      Algorithm %in% c("PSO1", "PSO2", "AT3-PSO1", "AT3-PSO2",
                       "AT5-PSO1", "AT5-PSO2",
                       "AT3-BBPSO", "AT3-BBPSOxp",
                       "AT5-BBPSO", "AT5-BBPSOxp")))

objname <- "sig2puk.new.mean"
p1 <- qplot(Time, logpost, color = Algorithm, geom = "line", size = I(1), facets = Nbhd ~.,
            data = filter(valueorder, Obj == objname))
objname <- "sig2puk.new.max"
p2 <- qplot(Time, logpost, color = Algorithm, geom = "line", size = I(1), facets = Nbhd ~.,
            data = filter(valueorder, Obj == objname))
grid.arrange(p1, p2, ncol = 2)


objname <- "sig2puk.mean"
p1 <- qplot(Time, logpost, color = Algorithm, geom = "line", facets = Nbhd~., size = I(1),
      data = subset(valueorder, Obj == objname & Nbhd %in% c("Global", "SS3") &
      Algorithm %in% c("PSO1", "AT3-PSO1", "AT5-PSO1")))


p2 <- qplot(Time, logpost, color = Algorithm, geom = "line", facets = Nbhd~., size = I(1),
      data = subset(valueorder, Obj == objname & Nbhd %in% c("Global", "SS3") &
      Algorithm %in% c("PSO2", "AT3-PSO2", "AT5-PSO2")))


grid.arrange(p1, p2, ncol = 2)

init <- c(spsample(datlist$sppoly, ndesign, "random")@coords)

optimtest <- optim(init, sig2fuk.mean2, datlist = datlist, method = "BFGS",
                   control = list(reltol = .Machine$double.eps, maxit = 1000))




sig2fuk.mean2 <- function(dd, datlist){
  dd2 <- matrix(dd, ncol = 2)
  poly <- datlist$sppoly@polygons$b@Polygons$a@coords
  check <- all(point.in.polygon(dd2[,1], dd2[,2], poly[,1], poly[,2]))
  if(check)
    out <- sig2fuk.mean(dd, datlist)
  else
    out <- Inf
  return(out)
}

sig2fuk.max2 <- function(dd, datlist){
  dd2 <- matrix(dd, ncol = 2)
  poly <- datlist$sppoly@polygons$b@Polygons$a@coords
  check <- all(point.in.polygon(dd2[,1], dd2[,2], poly[,1], poly[,2]))
  if(check)
    out <- sig2fuk.max(dd, datlist)
  else
    out <- Inf
  return(out)
}
##########################################################################

psoout <- read.csv("psosimsout.csv")
psoout <- psoout[psoout$rep < 5,]
psoout$obj <- mapvalues(psoout$obj, c("sig2fuk.max", "sig2fuk.mean"), c("sig2puk.max", "sig2puk.mean"))
niter <- max(psotimesum$time)

gaout <- read.csv("gasimsout.csv")
exout <- read.csv("exsimsout.csv")


psotimesum <- ddply(psoout, .(obj, algid, nbhd, time, type, parset, CF, style), summarise, mean=mean(logpost),
                median=median(logpost), sd=sd(logpost), min=min(logpost),
                q05=quantile(logpost, 0.05), q25=quantile(logpost, 0.25),
                q75=quantile(logpost, 0.75), q95=quantile(logpost, 0.95),
                max=max(logpost))




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

psosumsmol <- ddply(psotimesum, .(obj, algid, nbhd, time, mean, sd), namesfun)[,c(1,7,8,4,5,6)]

colnames(psosumsmol) <- c("Obj", "Algorithm", "Nbhd", "Time", "Mean", "SD")


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

nlook <- niter
psosumlast <- psosumorder[psosumorder$Time == nlook,]
## this prints the tables, however it's missing the first header.
## paste the following line (without the ##) just after the \begin{tabular} line
## \multicolumn{1}{l}{\begin{tabular}{lr} ; & Nbhd: \end{tabular}} & \multicolumn{2}{c}{Global} & \multicolumn{2}{c}{SS3} & \multicolumn{2}{c}{SS1} \\
for(k in unique(psosumlast$Obj)){
  ofk <- subset(psosumlast, Obj == k & !(Algorithm %in% c("DI-PSO1", "DI-PSO2", "DI-PSO1-CF", "DI-PSO2-CF")))[,-1]
  ## remove large values from table
  bigidx <- which(ofk$Mean > 10000)
  ofk$Mean[bigidx] <- NA
  ofk$SD[bigidx] <- NA
  ## pull out each type of nbhd and combine
  ofk1 <- subset(ofk, Nbhd == "Global")[,-c(2,3)]
  ofk2 <- subset(ofk, Nbhd == "SS3")[,-c(2,3)]
  ofk3 <- subset(ofk, Nbhd == "SS1")[,-c(2,3)]
  ofkout <- cbind(ofk1, ofk2[,-1]) #, ofk3[,-1])
  ##iidd <- c(3:4, 7:8, 11:12, 15:16, 19:20, 23:24)
  ##ofprint <- ofkout[-iidd,-c(3,5,7)]
  print(xtable(ofkout, digits = c(0,0,rep(c(2,2),2)), align = "ll|rr|rr",
               caption = paste("Simulation results for ", k, ". See text for description", sep = ""),
               label = paste("tab:psosim", k, sep="")),
        include.rownames=FALSE, sanitize.text.function=identity, hline.after=c(-1,seq(0, 24, 4)),
        size = "\\scriptsize")
}



objname <- "sig2puk.max"
q1 <- qplot(Time, Mean, color = Algorithm, geom = "line", facets = Nbhd~., size = I(1),
            data = subset(psosumorder, Time > 100 & Obj == objname & Nbhd %in% c("Global", "SS3") &
                                       Algorithm %in% c("PSO1", "PSO2", "AT3-PSO1", "AT3-PSO2")))
q2 <- qplot(Time, Mean, color = Algorithm, geom = "line", facets = Nbhd~., size = I(1),
            data = subset(psosumorder, Time > 100 & Obj == objname & Nbhd %in% c("Global", "SS3") &
                                       Algorithm %in% c("AT3-BBPSO", "AT3-BBPSOxp", "AT5-BBPSO", "AT5-BBPSOxp")))
grid.arrange(q1, q2, ncol = 2)

qplot(Time, Mean, color = Algorithm, geom = "line", facets = Nbhd~., size = I(1),
      data = subset(psosumorder, Time > 100 & Obj =="sig2puk.mean" & Nbhd %in% c("Global", "SS3") &
      Algorithm %in% c("PSO1-CF", "PSO2-CF", "AT3-PSO1-CF", "AT3-PSO2-CF",
                       "AT5-PSO1-CF", "AT5-PSO2-CF",
                       "AT3-BBPSO-CF", "AT3-BBPSOxp-CF",
                       "AT5-BBPSO-CF", "AT5-BBPSOxp-CF")))


p1 <- qplot(Time, Mean, color = Algorithm, geom = "line", facets = Nbhd~., size = I(1),
      data = subset(psosumorder, Time > 10 & Obj =="sig2puk.mean" & Nbhd %in% c("Global", "SS3") &
      Algorithm %in% c("PSO1", "AT3-PSO1", "AT5-PSO1")))


p2 <- qplot(Time, Mean, color = Algorithm, geom = "line", facets = Nbhd~., size = I(1),
      data = subset(psosumorder, Time > 10 & Obj =="sig2puk.mean" & Nbhd %in% c("Global", "SS3") &
      Algorithm %in% c("PSO2", "AT3-PSO2", "AT5-PSO2")))

library(gridExtra)

grid.arrange(p1, p2, ncol = 2)


psosumorder

