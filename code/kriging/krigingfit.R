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

### AIC/BIC are basically agnotstic between linear vs constant and spatial vs no.
### will assume spatial + linear


##### now estimate the variogram nonparametrically using two methods (ignores S)
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
datlist$sig2z <- 0
datlist$covfun <- expcov2
datlist$ss <- currloc
datlist$tt <- spsample(datlist$poly, 1000, "random")@coords
Cz.s <- Czfun(currloc, nrow(currloc), thetahat, 0, expcov2)
##+ diag(mean(S), length(S)) ## should we use that...?
datlist$invCz.s <- chol2inv(chol(Cz.s))
datlist$Cyy.s.t <- Cyyfun(currloc, datlist$tt, nrow(currloc), nrow(datlist$tt), thetahat, expcov2)
datlist$Cy.t <- drop(Cyyfun(matrix(c(0,0), nrow=1), matrix(c(0,0), nrow=1), 1, 1, thetahat, expcov2))

save(datlist, file = "datlist.Rdata")

p + geom_polygon(aes(longitude,latitude, group = Poly_Name), data = housgeom, fill = NA,
                 colour = "black") +
  geom_point(aes(u, v), data = data.frame(u = datlist$tt[,1]/(pi/180*6371),
                                          v = datlist$tt[,2]/(pi/180*6371)),
             colour = "blue", alpha = 0.6, shape="+", size = I(3)) 


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




source("../psofun.R")

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


