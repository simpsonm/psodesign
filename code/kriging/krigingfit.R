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

trendtype <- "1st"

## Fitting models with nugget fixed to zero
ml1 <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = TRUE, trend = trendtype)
ml2 <- likfit(houstongeo, ini = c(.1,00.5), fix.nugget = TRUE, trend = trendtype)
ml3 <- likfit(houstongeo, ini = c(10,5), fix.nugget = TRUE, trend = trendtype)
ml4 <- likfit(houstongeo, ini = c(10,0.5), fix.nugget = TRUE, trend = trendtype)
	 
## Fitting models estimated nugget
ml.n1 <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = FALSE, nugget = 1, trend = trendtype)
ml.n2 <- likfit(houstongeo, ini = c(.1,0.05), fix.nugget = FALSE, nugget = 1, trend = trendtype)
ml.n3 <- likfit(houstongeo, ini = c(10,5), fix.nugget = FALSE, nugget = 1, trend = trendtype)
ml.n4 <- likfit(houstongeo, ini = c(10,0.05), fix.nugget = FALSE, nugget = 1, trend = trendtype)
ml.n5 <- likfit(houstongeo, ini = c(.1,5), fix.nugget = FALSE, nugget = 1, trend = trendtype)

ml.n1 <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = FALSE, nugget = .1, trend = trendtype)
ml.n2 <- likfit(houstongeo, ini = c(.1,0.05), fix.nugget = FALSE, nugget = .1, trend = trendtype)
ml.n3 <- likfit(houstongeo, ini = c(10,5), fix.nugget = FALSE, nugget = .1, trend = trendtype)
ml.n4 <- likfit(houstongeo, ini = c(10,0.05), fix.nugget = FALSE, nugget = .1, trend = trendtype)
ml.n5 <- likfit(houstongeo, ini = c(.1,5), fix.nugget = FALSE, nugget = .1, trend = trendtype)

ml.n1 <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = FALSE, nugget = 10, trend = trendtype)
ml.n2 <- likfit(houstongeo, ini = c(.1,0.05), fix.nugget = FALSE, nugget = 10, trend = trendtype)
ml.n3 <- likfit(houstongeo, ini = c(10,5), fix.nugget = FALSE, nugget = 10, trend = trendtype)
ml.n4 <- likfit(houstongeo, ini = c(10,0.05), fix.nugget = FALSE, nugget = 10, trend = trendtype)
ml.n5 <- likfit(houstongeo, ini = c(.1,5), fix.nugget = FALSE, nugget = 10, trend = trendtype)

## basically, nugget is zero 

ml <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = TRUE, trend = trendtype)
ml.n <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = FALSE, nugget = 1, trend = trendtype)

cbind(ml$beta +  qnorm(0.025) * sqrt(diag(as.matrix(ml$beta.var))),
      ml$beta +  qnorm(0.975) * sqrt(diag(as.matrix(ml$beta.var))))

cbind(ml$beta +  qnorm(0.005) * sqrt(diag(as.matrix(ml$beta.var))),
      ml$beta +  qnorm(0.995) * sqrt(diag(as.matrix(ml$beta.var))))

cbind(ml$beta +  qnorm(0.1) * sqrt(diag(as.matrix(ml$beta.var))),
      ml$beta +  qnorm(0.9) * sqrt(diag(as.matrix(ml$beta.var))))

## confidence intervals suggest linear in y axis is reasonable even at 99% level
## but x axis is not a significant predictor even at 80% level

## not shown: 2nd order polynomial clearly seems to be an overfit




##### now estimate the variogram nonparametrically using two methods


## for use in objective function - for checking whether candidate points are in the poly
harrispoly <- cbind(housgeom$longitude, housgeom$latitude)*pi/180*6371
currloc <- cbind(houston$u, houston$v)

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

## obtain ML estimates
ml <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = TRUE, trend = trendtype)
ml.n <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = FALSE, nugget = 0, trend = trendtype)

# Now, plotting fitted models against empirical variogram
par(mfrow = c(2,1))
plot(bin1, main = "classical variogram")
lines(ml)
lines(ml.n, lty = 2)
legend(40, 15, legend=c("ML","ML.N"),lty=c(1,2),lwd=c(1,1), cex=0.7)
plot(bin2, main = "modulus variogram")
lines(ml)
lines(ml.n, lty = 2)
legend(40, 15, legend=c("ML","ML.N"),lty=c(1,2),lwd=c(1,1), cex=0.7)

betahat <- ml$beta
tau2hat <- ml$tausq
sig2hat <- ml$sigmasq
phihat <- ml$phi
varbetahat <- ml$beta.var
thetahat <- c(phihat, 0, sig2hat)

basemap <- get_map(location = "houston", zoom = 9, maptype = 'terrain')
p <- ggmap(basemap) + geom_point(aes(Longitude, Latitude), data = houston, size = I(3), alpha=0.6)

N.sim <- 1000
datlist <- list()
datlist$poly <- Polygon(harrispoly)
datlist$theta <- thetahat
datlist$sig2z <- 0
datlist$covfun <- expcov
datlist$ss <- currloc
datlist$tt <- spsample(datlist$poly, 1000, "random")@coords
datlist$invCz.s <- chol2inv(chol(Czfun(currloc, nrow(currloc), thetahat, 0, expcov)))
datlist$Cyy.s.t <- Cyyfun(currloc, datlist$tt, nrow(currloc), nrow(datlist$tt), thetahat, expcov)
datlist$Cy.t <- drop(Cyyfun(matrix(c(0,0), nrow=1), matrix(c(0,0), nrow=1), 1, 1, thetahat, expcov))

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

ndesign <- 1
npar <- 2*ndesign
inits <- list()
inits[[1]] <- replicate(nswarm, c(spsample(datlist$poly, ndesign, "random")@coords))
idxs <- matrix(replicate(nswarm, sample(1:nrow(datlist$tt), ndesign)), nrow = ndesign)
inits2 <- matrix(0, npar, nswarm)
for(i in 1:nswarm){
  inits2[,i] <- c(datlist$tt[idxs[,i],])
}
inits[[2]] <- inits2
vals <- apply(t(datlist$tt), 2, negsig2sk.mean, datlist = datlist)
idxs2 <- order(vals, decreasing = TRUE)[1:ndesign]
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


par(mfrow=c(1,1))
plot(ts(meanpso$maxes), ylim = c(min(meanpso$maxes, meanpso2$maxes),
                                 max(meanpso$maxes, meanpso2$maxes)))
lines(ts(meanpso2$maxes), col = "red")


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


