library(geoR)
library(dplyr)
library(ggplot2)
library(rgdal)
library(ggmap)
source("krigingfun.R")
library(shapefiles)
library(maptools)
source("poly_coords_function.R")

houston <- read.csv("houstonout.csv")
houston$se <- 0

## county level shapefile
coshape <- readOGR(dsn="shape", layer="05000")
names(coshape)[1] <- "ID"
cogeom <- poly_coords(coshape)
housgeom <- filter(cogeom, NAME =="Harris")

## shapefile is in a different projection; convert to lon + lat
d <- data.frame( x = housgeom$PolyCoordsY, y = housgeom$PolyCoordsX)
coordinates(d) <- c("x", "y")
proj4string(d) <- proj4string(coshape)
CRS.new <- CRS("+init=epsg:4326") # WGS 84, i.e. world
d2 <- spTransform(d, CRS.new)
housgeom$longitude <- d2@coords[,1]
housgeom$latitude <- d2@coords[,2]

## plot of houston area, Harris county, and current monitoring sites
basemap <- get_map(location = "houston", zoom = 9, maptype = 'terrain')
p <- ggmap(basemap) + geom_point(aes(Longitude, Latitude), data = houston, size = I(3), alpha=0.6)
p + geom_polygon(aes(longitude,latitude, group = Poly_Name), data = housgeom, fill = NA,
                 colour = "black")

## convert longitude & lattitude to approximate kilometers
houston$u <- houston$Longitude*pi/180*6371
houston$v <- houston$Latitude*pi/180*6371

houstongeo <- list(coords = as.matrix(select(houston, u, v)),
                   data = houston$avg)
colnames(houstongeo[[1]]) <- NULL


## initial fits for good starting values, but these ignore S and the measurement error variance
ml.lin <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = TRUE, nugget = 0, trend = "1st")
ml.cte <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = TRUE, nugget = 0, trend = "cte")

parslin <- c(ml.lin$beta, 0.0001, log(ml.lin$sigmasq), ml.lin$phi)
parscte <- c(ml.cte$beta, 0.0001, log(ml.cte$sigmasq), ml.cte$phi)
parsquad <- c(ml.lin$beta, c(0,0,0), 0.0001, log(ml.lin$sigmasq), ml.lin$phi)
                       
confit <- mylikefit(parscte, "cte", houston)
linfit <- mylikefit(parslin, "lin", houston)
quadfit <- mylikefit(parsquad, "quad", houston)

confit$ics
linfit$ics
quadfit$ics

### AIC/BIC are basically agnotstic between linear vs constant and spatial vs no.
### will assume spatial + linear


##### now estimate the variogram nonparametrically using two methods
trendtype <- "1st"
cloud1 <- variog(houstongeo, option = "cloud", trend=trendtype)
cloud2 <- variog(houstongeo, option = "cloud", estimator.type = "modulus",
                 trend=trendtype)
bin1 <- variog(houstongeo, trend=trendtype)
bin2 <- variog(houstongeo, estimator.type = "modulus", trend=trendtype)

par(mfrow = c(2, 2))
plot(cloud1, main = "classical estimator")
plot(cloud2, main = "modulus estimator")
plot(bin1, main = "classical estimator")
plot(bin2, main = "modulus estimator")

# Now, plotting fitted models against empirical variogram
par(mfrow = c(2,1))
plot(bin1, main = "classical variogram")
seqs <- seq(0, 150, length.out=1000)
lines(seqs, linfit$par[4] + linfit$par[5]*(1 - exp(-seqs/linfit$par[6])), col = "red")
plot(bin2, main = "modulus variogram")
lines(seqs, linfit$par[4] + linfit$par[5]*(1 - exp(-seqs/linfit$par[6])), col = "red")

betahat <- linfit$par[1:3]
tau2hat <- linfit$par[4]
sig2hat <- linfit$par[5]
phihat <- linfit$par[6]
thetahat <- c(sig2hat, phihat)

basemap <- get_map(location = "houston", zoom = 9, maptype = 'terrain')
p <- ggmap(basemap) + geom_point(aes(Longitude, Latitude), data = houston, size = I(3), alpha=0.6)


## for use in objective function - for checking whether candidate points are in the poly
harrispoly <- cbind(housgeom$longitude, housgeom$latitude)*pi/180*6371
currloc <- cbind(houston$u, houston$v)



N.sim <- 1000
datlist <- list()
datlist$poly <- Polygon(harrispoly)
datlist$theta <- thetahat
datlist$sig2z <- tau2hat
datlist$ss <- currloc
datlist$tt <- spsample(datlist$poly, 1000, "random")@coords
N.t <- nrow(datlist$tt)
D.t <- matrix(0, N.t, N.t)
tt <- datlist$tt
for(i in 2:N.t){
  for(j in 1:(i-1)){
    D.t[i,j] <- sqrt(sum((tt[i,] - tt[j,])^2))
    D.t[j,i] <- D.t[i,j]
  }
}
datlist$D.t <- D.t
N.s <- nrow(datlist$ss)
D.s <- matrix(0, N.s, N.s)
ss <- datlist$ss
for(i in 2:N.s){
  for(j in 1:(i-1)){
    D.s[i,j] <- sqrt(sum((ss[i,] - ss[j,])^2))
    D.s[j,i] <- D.s[i,j]
  }
}
datlist$D.s <- D.s
D.s.t <- matrix(0, N.s, N.t)
for(i in 1:N.s){
  for(j in 1:N.t){
    D.s.t[i,j] <- sqrt(sum((datlist$ss[i,] - tt[j,])^2))
  }
}
datlist$D.s.t <- D.s.t
Cz.s <- diag(datlist$sig2z, N.s) + datlist$theta[1]*exp(-D.s/datlist$theta[2])
datlist$invCz.s <- chol2inv(chol(Cz.s))
datlist$Cyy.s.t <- datlist$theta[1]*exp(-D.s.t/datlist$theta[2])
datlist$Cy.t <- datlist$theta[1]

save(datlist, file = "datlist.Rdata")


p + geom_polygon(aes(longitude,latitude, group = Poly_Name), data = housgeom, fill = NA,
                 colour = "black") +
  geom_point(aes(u, v), data = data.frame(u = datlist$tt[,1]/(pi/180*6371),
                                          v = datlist$tt[,2]/(pi/180*6371)),
             colour = "blue", alpha = 0.6, shape="+", size = I(3)) 



source("../psofun.R")

niter <- 100
nswarm <- 20
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496
nbhdnames <- c("ring-1", "ring-3", "global")
rates <- c(0.5)
dfs <- c(5)
ccc <- c(0.1)
alpha <- .2*niter
beta <- 1
nbhd <- list()
nbhd[[1]] <- sapply(1:nswarm, function(x){return( (x + -1:1 - 1)%%nswarm + 1)}) ## ring-1
nbhd[[2]] <- sapply(1:nswarm, function(x){return( (x + -3:3 - 1)%%nswarm + 1)}) ## ring-3
nbhd[[3]] <- matrix(1:nswarm, ncol=nswarm, nrow=nswarm) ## global

vals <- apply(t(datlist$tt), 2, negsig2sk.mean, datlist = datlist)

ndesign <- 5
idxs2 <- order(vals, decreasing = TRUE)[1:ndesign]
npar <- 2*ndesign
inits <- list()
inits[[1]] <- replicate(nswarm, c(spsample(datlist$poly, ndesign, "random")@coords))
idxs <- matrix(replicate(nswarm, sample(1:nrow(datlist$tt), ndesign)), nrow = ndesign)
inits2 <- matrix(0, npar, nswarm)
for(i in 1:nswarm){
  inits2[,i] <- c(datlist$tt[idxs[,i],])
}
inits[[2]] <- inits2
inits[[2]][,1] <- c(datlist$tt[idxs2,])


system.time(meanpso <- pso(niter, nswarm, inertia, cognitive, social, inits[[1]], nbhd[[1]],
                           negsig2uk.mean, datlist = datlist))

system.time(minpso <- pso(niter, nswarm, inertia, cognitive, social, inits[[1]], nbhd[[1]],
                          negsig2uk.min, datlist = datlist))


system.time(meanpso2 <- pso(niter, nswarm, inertia, cognitive, social, inits[[2]], nbhd[[1]],
                           negsig2uk.mean, datlist = datlist))

system.time(minpso2 <- pso(niter, nswarm, inertia, cognitive, social, inits[[2]], nbhd[[1]],
                          negsig2uk.min, datlist = datlist))


par(mfrow=c(2,1))
plot(ts(meanpso$maxes), ylim = c(min(meanpso$maxes, meanpso2$maxes),
                                 max(meanpso$maxes, meanpso2$maxes)))
lines(ts(meanpso2$maxes), col = "red")
plot(ts(minpso$maxes), ylim = c(min(minpso$maxes, minpso2$maxes),
                                max(minpso$maxes, minpso2$maxes)))
lines(ts(minpso2$maxes), col = "red")


p + geom_polygon(aes(longitude,latitude, group = Poly_Name), data = housgeom, fill = NA,
                       colour = "black") +
  geom_point(aes(u, v), data = data.frame(u = datlist$tt[,1]/(pi/180*6371),
                                          v = datlist$tt[,2]/(pi/180*6371)),
             colour = "blue", alpha = 0.6, shape="+", size = I(3)) +
  geom_point(aes(u, v), data = data.frame(u = meanpso$argmax[1:ndesign]/(pi/180*6371),
                                          v = meanpso$argmax[ndesign + 1:ndesign]/(pi/180*6371)),
             colour = "red", shape="X", size = I(5)) +
  geom_point(aes(u, v), data = data.frame(u = minpso$argmax[1:ndesign]/(pi/180*6371),
                                          v = minpso$argmax[ndesign + 1:ndesign]/(pi/180*6371)),
             colour = "blue", shape="X", size = I(5)) + 
  geom_point(aes(u, v), data = data.frame(u = meanpso2$argmax[1:ndesign]/(pi/180*6371),
                                          v = meanpso2$argmax[ndesign + 1:ndesign]/(pi/180*6371)),
             colour = "red", shape="O", size = I(5)) +
  geom_point(aes(u, v), data = data.frame(u = minpso2$argmax[1:ndesign]/(pi/180*6371),
                                          v = minpso2$argmax[ndesign + 1:ndesign]/(pi/180*6371)),
             colour = "blue", shape="O", size = I(5))





ncand <- 10000
nnbor <- 5
npoints <- 5
poly <- datlist$poly@coords

system.time(exchtest1 <- exch(ncand, sig2uk.mean, poly, nnbor, npoints, datlist = datlist))

nnbor <- 10
system.time(exchtest2 <- exch(ncand, sig2uk.mean, poly, nnbor, npoints, datlist = datlist))
