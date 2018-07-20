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
p2 <- p +
    geom_polygon(aes(longitude,latitude, group = Poly_Name),
                 data = housgeom, fill = NA,
                 colour = "black")

ggsave("houston.png", p2, scale = 1)

## convert longitude & lattitude to approximate kilometers
houston$u <- houston$Longitude*pi/180*6371
houston$v <- houston$Latitude*pi/180*6371

houstongeo <- list(coords = as.matrix(select(houston, u, v)),
                   data = houston$avg)
colnames(houstongeo[[1]]) <- NULL


## MLEs, trying different starting values
ml.lin <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = TRUE, nugget = 0, trend = "1st")
ml.cte <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = TRUE, nugget = 0, trend = "cte")

ml.lin2 <- likfit(houstongeo, ini = c(1,0.05), fix.nugget = TRUE, nugget = 0, trend = "1st")
ml.cte2 <- likfit(houstongeo, ini = c(1,0.05), fix.nugget = TRUE, nugget = 0, trend = "cte")

ml.lin3 <- likfit(houstongeo, ini = c(10,0.05), fix.nugget = TRUE, nugget = 0, trend = "1st")
ml.cte3 <- likfit(houstongeo, ini = c(10,0.05), fix.nugget = TRUE, nugget = 0, trend = "cte")

ml.lin4 <- likfit(houstongeo, ini = c(.0010,0.05), fix.nugget = TRUE, nugget = 0, trend = "1st")
ml.cte4 <- likfit(houstongeo, ini = c(.0010,0.05), fix.nugget = TRUE, nugget = 0, trend = "cte")


ml.lin.nug <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = FALSE, trend = "1st")
ml.cte.nug <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = FALSE, trend = "cte")

ml.lin.nug2 <- likfit(houstongeo, ini = c(1,0.005), fix.nugget = FALSE, trend = "1st")
ml.cte.nug2 <- likfit(houstongeo, ini = c(1,0.005), fix.nugget = FALSE, trend = "cte")

ml.lin.nug3 <- likfit(houstongeo, ini = c(100,0.5), fix.nugget = FALSE, trend = "1st")
ml.cte.nug3 <- likfit(houstongeo, ini = c(100,0.5), fix.nugget = FALSE, trend = "cte")

ml.lin.nug4 <- likfit(houstongeo, ini = c(.001,0.5), fix.nugget = FALSE, trend = "1st")
ml.cte.nug4 <- likfit(houstongeo, ini = c(.001,0.5), fix.nugget = FALSE, trend = "cte")

ml.lin.nug5 <- likfit(houstongeo, ini = c(.001,0.5), fix.nugget = FALSE, nugget = 100, trend = "1st")
ml.cte.nug5 <- likfit(houstongeo, ini = c(.001,0.5), fix.nugget = FALSE, nugget = 100, trend = "cte")

## tried a bunch of different starting values

ml <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = TRUE, trend = "1st")
ml.n <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = FALSE, nugget = 1, trend = "1st")

## AIC/BIC favors no nugget, which is estimated to be zero when its included in the model.

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


## Now, plotting fitted models against empirical variogram
par(mfrow = c(2,1))
plot(bin1, main = "classical variogram")
lines(ml)
lines(ml.n, lty = 2)
legend(40, 15, legend=c("ML","ML.N"),lty=c(1,2),lwd=c(1,1), cex=0.7)
plot(bin2, main = "modulus variogram")
lines(ml)
lines(ml.n, lty = 2)
legend(40, 15, legend=c("ML","ML.N"),lty=c(1,2),lwd=c(1,1), cex=0.7)

## store parameter estimates
betahat <- ml$beta
tau2hat <- ml$tausq
sig2hat <- ml$sigmasq
phihat <- ml$phi
varbetahat <- ml$beta.var
thetahat <- c(sig2hat, phihat)


## for use in objective function - for checking whether candidate points are in the poly
harrispoly <- cbind(housgeom$longitude, housgeom$latitude)*pi/180*6371
currloc <- cbind(houston$u, houston$v)


## create prediction grid
N.grid <- 2000
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

## create list of inputs to objective functions
datlist <- list()
datlist$poly <- Polygon(harrispoly) ## Harris County polygon
datlist$theta <- thetahat           ## sig2 and phi estimates
datlist$sig2z <- tau2hat            ## sig2z estimate (set to 0)
datlist$ss <- currloc               ## current observation locations
datlist$tt <- grid                  ## prediction locations
N.t <- nrow(datlist$tt)             ## number of prediction locations
D.t <- matrix(0, N.t, N.t)          ## prediction location distance matrix
tt <- datlist$tt
for(i in 2:N.t){
  for(j in 1:(i-1)){
    D.t[i,j] <- sqrt(sum((tt[i,] - tt[j,])^2))
    D.t[j,i] <- D.t[i,j]
  }
}
datlist$D.t <- D.t
N.s <- nrow(datlist$ss)             ## number of current observation locations
D.s <- matrix(0, N.s, N.s)          ## observation location distance matrix
ss <- datlist$ss
for(i in 2:N.s){
  for(j in 1:(i-1)){
    D.s[i,j] <- sqrt(sum((ss[i,] - ss[j,])^2))
    D.s[j,i] <- D.s[i,j]
  }
}
datlist$D.s <- D.s
D.s.t <- matrix(0, N.s, N.t)        ## observation vs. prediction distance matrix
for(i in 1:N.s){
  for(j in 1:N.t){
    D.s.t[i,j] <- sqrt(sum((datlist$ss[i,] - tt[j,])^2))
  }
}
datlist$D.s.t <- D.s.t
Cz.s <- diag(datlist$sig2z, N.s) + datlist$theta[1]*exp(-D.s/datlist$theta[2])
datlist$invCz.s <- chol2inv(chol(Cz.s))  ## inverse of cov(Z) on observation locations
datlist$Cyy.s.t <- datlist$theta[1]*exp(-D.s.t/datlist$theta[2])  ## cov(Y(s), Y(t)) matrix
datlist$Cy.t <- datlist$theta[1]         ## variance at a point

save(datlist, file = "datlist.Rdata")


## map initialization
basemap <- get_map(location = "houston", zoom = 9, maptype = 'terrain')
p <- ggmap(basemap) + geom_point(aes(Longitude, Latitude), data = houston, size = I(3), alpha=0.6)

## prediction locations map
p + geom_polygon(aes(longitude,latitude, group = Poly_Name), data = housgeom, fill = NA,
                 colour = "black") +
  geom_point(aes(u, v), data = data.frame(u = datlist$tt[,1]/(pi/180*6371),
                                          v = datlist$tt[,2]/(pi/180*6371)),
             colour = "blue", alpha = 0.6, shape="+", size = I(3))



## test some optimzation routines
source("../psofun.R")

niter <- 500
nswarm <- 40
inertias <- c(0.7298, 1/(log(2)*2))
cognitives <- c(1.496, log(2) + 1/2)
socials <- c(1.496, log(2) + 1/2)
nnbor <- 3
ndesign <- 5
lower <- rep(apply(datlist$poly@coords, 2, min), each = ndesign)
upper <- rep(apply(datlist$poly@coords, 2, max), each = ndesign)


psotime <- system.time(spsoCI1 <- spso(niter, nswarm, nnbor, inertias[1], cognitives[1], socials[1], sig2fuk.mean, lower, upper, style = "CI", CF = FALSE, datlist = datlist))

spsoCI2 <- spso(niter, nswarm, nnbor, inertias[2], cognitives[2], socials[2], sig2fuk.mean,
                lower, upper, style = "CI", CF = FALSE, datlist = datlist)

spsoCI1.max <- spso(niter, nswarm, nnbor, inertias[1], cognitives[1], socials[1], sig2fuk.max,
                lower, upper, style = "CI", CF = FALSE, datlist = datlist)

spsoCI2.max <- spso(niter, nswarm, nnbor, inertias[2], cognitives[2], socials[2], sig2fuk.max,
                lower, upper, style = "CI", CF = FALSE, datlist = datlist)

par(mfrow=c(2,1))
plot(ts(spsoCI1$values))
lines(ts(spsoCI2$values), col="red")
plot(ts(spsoCI1.max$values))
lines(ts(spsoCI2.max$values), col="red")



spsoCI2 <- spso(niter, nswarm, nnbor, inertia, cognitive, social, sig2fuk.mean,
                lower, upper, style = "AT", CF = FALSE, datlist = datlist)

spsoCI.max <- spso(niter, nswarm, nnbor, inertia, cognitive, social, sig2fuk.max, lower, upper,
                   style = "CI", CF = FALSE, datlist = datlist)

spsoCI.max2 <- spso(niter, nswarm, nnbor, inertia, cognitive, social, sig2fuk.max,
                   lower, upper, style = "AT", CF = FALSE, datlist = datlist)

nbatch <- 1
nbatch2 <- 2
nchrome <- 2
nrun <- ndesign
mutvar <- 2
mutrate <- .1
mutvar2 <- 1
mutrate2 <- .01
## these seem like good parameter values!

gatime1 <- system.time(ga1 <- ga(niter, nbatch, floor(nswarm/2), nchrome, nrun,
                                 mutvar, mutrate, lower, upper,
                                 sig2fuk.mean, datlist=datlist))

gatime2 <- system.time(ga2 <- ga(niter/nbatch2, nbatch2, floor(nswarm/2), nchrome, nrun,
                                 mutvar, mutrate, lower, upper,
                                 sig2fuk.mean, datlist=datlist))

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







ncand <- 2000
nnbor <- 5
npoints <- 5
poly <- datlist$poly@coords

system.time(exchtest1 <- exch(ncand, sig2fuk.mean, poly, nnbor, npoints, datlist = datlist))

nnbor <- 7
system.time(exchtest2 <- exch(ncand, sig2fuk.mean, poly, nnbor, npoints, datlist = datlist))


nnbor <- 9
system.time(exchtest3 <- exch(ncand, sig2fuk.mean, poly, nnbor, npoints, datlist = datlist))
