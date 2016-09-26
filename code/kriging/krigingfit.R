library(geoR)
library(dplyr)
library(ggplot2)
library(rgdal)
library(ggmap)
library(gridExtra)
source("krigingfun.R")

houston <- read.csv("houstonout.csv")

## convert coordinates to a projection suited for this area of texas
x <- houston$Longitude
y <- houston$Latitude
d <- data.frame(lon=x, lat=y)
coordinates(d) <- c("lon", "lat")
proj4string(d) <- CRS("+init=epsg:4326") # WGS 84, i.e. world
CRS.new <- CRS("+init=epsg:2256") # Texas 
d2 <- spTransform(d, CRS.new)

## rescale so betas & phi are easier to fit
houston$u <- d2@coords[,1]/1000000  
houston$v <- d2@coords[,2]/1000000  

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

## basically, nugget is zero and estiamtes consistent across starting values

ml <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = TRUE, trend = trendtype)
ml.n <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = FALSE, nugget = 1, trend = trendtype)

cbind(ml.n$beta +  qnorm(0.025) * sqrt(diag(as.matrix(ml.n$beta.var))),
      ml.n$beta +  qnorm(0.975) * sqrt(diag(as.matrix(ml.n$beta.var))))

cbind(ml$beta +  qnorm(0.025) * sqrt(diag(as.matrix(ml$beta.var))),
      ml$beta +  qnorm(0.975) * sqrt(diag(as.matrix(ml$beta.var))))

## confidence intervals are basically the same - linear in X, though not Y

############################################################
## now try 2nd order polynomial
trendtype <- "2nd"

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


## now higly sensitive to starting values, but nugget=0 is still MLE

ml <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = TRUE, trend = trendtype)
ml.n <- likfit(houstongeo, ini = c(10,0.05), fix.nugget = FALSE, nugget = 0, trend = trendtype)

cbind(ml.n$beta +  qnorm(0.025) * sqrt(diag(as.matrix(ml.n$beta.var))),
      ml.n$beta +  qnorm(0.975) * sqrt(diag(as.matrix(ml.n$beta.var))))

cbind(ml$beta +  qnorm(0.025) * sqrt(diag(as.matrix(ml$beta.var))),
      ml$beta +  qnorm(0.975) * sqrt(diag(as.matrix(ml$beta.var))))

## confidence intervals now think nothing is significant - so stay linear

## now estimate the variogram nonparametrically using two methods
trendtype <- "1st"
cloud1 <- variog(houstongeo, option = "cloud", max.dist = 1, trend=trendtype)
cloud2 <- variog(houstongeo, option = "cloud", estimator.type = "modulus",
                 max.dist = 1, trend=trendtype)
bin1 <- variog(houstongeo, uvec = seq(0, 1, l = 11), trend=trendtype)
bin2 <- variog(houstongeo, uvec = seq(0, 1, l = 11),
               estimator.type = "modulus", trend=trendtype)
par(mfrow = c(2, 2))
plot(cloud1, main = "classical estimator")
plot(cloud2, main = "modulus estimator")
plot(bin1, main = "classical estimator")
plot(bin2, main = "modulus estimator")

## obtain ML estimates
ml <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = TRUE, trend = trendtype)
ml.n <- likfit(houstongeo, ini = c(10,0.05), fix.nugget = FALSE, nugget = 0, trend = trendtype)

# Now, plotting fitted models against empirical variogram
par(mfrow = c(2,1))
plot(bin1, main = "classical variogram", ylim=c(0,30))
lines(ml, max.dist = 1)
lines(ml.n, lty = 2, max.dist = 1)
legend(0.2, 15, legend=c("ML","ML.N"),lty=c(1,2),lwd=c(1,1), cex=0.7)
plot(bin2, main = "modulus variogram", ylim=c(0,30))
lines(ml, max.dist = 1)
lines(ml.n, lty = 2, max.dist = 1)
legend(0.2, 15, legend=c("ML","ML.N"),lty=c(1,2),lwd=c(1,1), cex=0.7)

## looks reasonable on the classical variogram, though a little off on the modulus one



basemap1 <- get_map(location = c(lon = -95.27566, lat = 29.75722), zoom = 9, maptype = 'terrain')
p1 <- ggmap(basemap1) + geom_point(aes(Longitude, Latitude), data = houston,
                                  size = I(3), alpha=0.4)


basemap2 <- get_map(location = c(lon = -95.5, lat = 29.7), zoom = 9, maptype = 'terrain')
p2 <- ggmap(basemap2) + geom_point(aes(Longitude, Latitude), data = houston,
                                  size = I(3), alpha=0.4)

grid.arrange(p1, p2, ncol = 2)

xs <- c(-95.75, -95.75, -95, -95)
ys <- c(29.5, 30, 30, 29.5)
dat <- data.frame(xs = xs, ys = ys)
p2 + geom_polygon(aes(xs, ys), data = dat, fill = NA, colour = "black")


## todo: convet bounds to transformed coordinates
## program up objective functions for a new site
## let PSO rip!

