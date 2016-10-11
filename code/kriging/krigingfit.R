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



N.grid <- 5000
ndim <- ceiling(sqrt(N.grid))
N.grid <- ndim^2
mins <- apply(harrispoly, 2, min)
maxes <- apply(harrispoly, 2, max)
by.x <- (maxes[1] - mins[1])/ndim
by.y <- (maxes[2] - mins[2])/ndim
cand <- expand.grid(seq(from = mins[1], by = by.x, length.out = ndim),
                    seq(from = mins[2], by = by.y, length.out = ndim))
checks <- apply(cand, 1, function(x, poly) {
  point.in.polygon(x[1], x[2], poly[,1], poly[,2])}, poly = harrispoly)
grid <- as.matrix(cand[checks==1,])




datlist <- list()
datlist$poly <- Polygon(harrispoly)
datlist$theta <- thetahat
datlist$sig2z <- tau2hat
datlist$ss <- currloc
datlist$tt <- grid
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
nswarm <- 40
inertia <- 0.7298
cognitive <- 1.496
social <- 1.496
nnbor <- 3
ndesign <- 5
lower <- rep(apply(datlist$poly@coords, 2, min), each = ndesign)
upper <- rep(apply(datlist$poly@coords, 2, max), each = ndesign)


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

spsoCI <- spso(niter, nswarm, nnbor, inertia, cognitive, social, sig2fuk.mean, lower, upper,
               style = "CI", CF = FALSE, datlist = datlist)

spsoCI2 <- spso(niter, nswarm, nnbor, inertia, cognitive, social, sig2fuk.mean,
                lower, upper, style = "AT", CF = FALSE, datlist = datlist)

spsoCI.max <- spso(niter, nswarm, nnbor, inertia, cognitive, social, sig2fuk.max, lower, upper,
                   style = "CI", CF = FALSE, datlist = datlist)

spsoCI.max2 <- spso(niter, nswarm, nnbor, inertia, cognitive, social, sig2fuk.max,
                   lower, upper, style = "AT", CF = FALSE, datlist = datlist)

nbatch <- 1
nchrome <- 2
nrun <- ndesign
mutvar <- 2
mutrate <- .1
mutvar2 <- 1
mutrate2 <- .01
## these seem like good parameter values!

ga1 <- ga(niter, nbatch, floor(nswarm/2), nchrome, nrun, mutvar, mutrate, lower, upper,
          sig2fuk.mean, datlist=datlist)

ga2 <- ga(niter, nbatch, floor(nswarm/2), nchrome, nrun, mutvar2, mutrate2, lower, upper,
          sig2fuk.mean, datlist=datlist)

ga3 <- ga(niter, nbatch, floor(nswarm/2), nchrome, nrun, mutvar2, mutrate, lower, upper,
          sig2fuk.mean, datlist=datlist)

ga4 <- ga(niter, nbatch, floor(nswarm/2), nchrome, nrun, mutvar, mutrate2, lower, upper,
          sig2fuk.mean, datlist=datlist)

c(ga1$value, ga2$value, ga3$value, ga4$value)

par(mfrow=c(2,1))
plot(ts(spsoCI$values), ylim = c(min(spsoCI$values, spsoCI2$values),
                                 max(spsoCI2$values)))
lines(ts(spsoCI2$values), col = "red")
plot(ts(spsoCI.max$values), ylim = c(min(spsoCI.max$values, spsoCI.max2$values),
                                     max(spsoCI.max$values)))
lines(ts(spsoCI.max2$values), col = "red")

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
  geom_point(aes(u, v), data = data.frame(u = spsoCI$par[1:ndesign]/(pi/180*6371),
                                          v = spsoCI$par[ndesign + 1:ndesign]/(pi/180*6371)),
             colour = "red", shape="X", size = I(5)) +
  geom_point(aes(u, v), data = data.frame(u = spsoCI2$par[1:ndesign]/(pi/180*6371),
                                          v = spsoCI2$par[ndesign + 1:ndesign]/(pi/180*6371)),
             colour = "red", shape="O", size = I(5)) +
  geom_point(aes(u, v), data = data.frame(u = spsoCI.max$par[1:ndesign]/(pi/180*6371),
                                          v = spsoCI.max$par[ndesign + 1:ndesign]/(pi/180*6371)),
             colour = "blue", shape="X", size = I(5)) +
  geom_point(aes(u, v), data = data.frame(u = spsoCI.max2$par[1:ndesign]/(pi/180*6371),
                                          v = spsoCI.max2$par[ndesign + 1:ndesign]/(pi/180*6371)),
             colour = "blue", shape="O", size = I(5)) 







ncand <- 10000
nnbor <- 5
npoints <- 5
poly <- datlist$poly@coords

system.time(exchtest1 <- exch(ncand, sig2uk.mean, poly, nnbor, npoints, datlist = datlist))

nnbor <- 10
system.time(exchtest2 <- exch(ncand, sig2uk.mean, poly, nnbor, npoints, datlist = datlist))
