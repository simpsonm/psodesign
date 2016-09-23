library(geoR)
library(dplyr)
source("krigingfun.R")

houston <- read.csv("houstonout.csv")

houstongeo <- list(coords = as.matrix(select(houston, Latitude, Longitude)),
                   data = houston$avg)
colnames(houstongeo[[1]]) <- NULL

cloud1 <- variog(houstongeo, option = "cloud", max.dist = 1)
cloud2 <- variog(houstongeo, option = "cloud", estimator.type = "modulus", max.dist = 1)
bin1 <- variog(houstongeo, uvec = seq(0, 1, l = 11))
bin2 <- variog(houstongeo, uvec = seq(0, 1, l = 11),
               estimator.type = "modulus")
par(mfrow = c(2, 2))
plot(cloud1, main = "classical estimator")
plot(cloud2, main = "modulus estimator")
plot(bin1, main = "classical estimator")
plot(bin2, main = "modulus estimator")


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
## Fitting models with nugget fixed to zero
ml <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = TRUE, trend = trendtype)
reml <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = TRUE, lik.method = "RML", trend = trendtype)
ols <- variofit(bin1, ini = c(1,0.5), fix.nugget = TRUE, weights="equal")
wls <- variofit(bin1, ini = c(1,0.5), fix.nugget = TRUE)
	 
## Fitting models estimated nugget
ml.n <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = FALSE, nugget = 1, trend = trendtype)
reml.n <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = FALSE, nugget = 1,
                 lik.method = "RML", trend = trendtype)
ols.n <- variofit(bin1, ini = c(1,0.5), fix.nugget = FALSE, weights="equal")
wls.n <- variofit(bin1, ini = c(1,0.5), fix.nugget = FALSE)


cbind(ml.n$beta +  qnorm(0.025) * sqrt(diag(as.matrix(ml.n$beta.var))),
      ml.n$beta +  qnorm(0.975) * sqrt(diag(as.matrix(ml.n$beta.var))))

cbind(ml$beta +  qnorm(0.025) * sqrt(diag(as.matrix(ml$beta.var))),
      ml$beta +  qnorm(0.975) * sqrt(diag(as.matrix(ml$beta.var))))


# Now, plotting fitted models against empirical variogram
par(mfrow = c(2,1))
plot(bin1, main = expression(paste("fixed ", tau^2 == 0)))
lines(ml, max.dist = 1)
lines(reml, lwd = 2, max.dist = 1)
lines(ols, lty = 2, max.dist = 1)
lines(wls, lty = 2, lwd = 2, max.dist = 1)
legend(0.2, .04, legend=c("ML","REML","OLS","WLS"),lty=c(1,1,2,2),lwd=c(1,2,1,2), cex=0.7)
plot(bin1, main = expression(paste("estimated  ", tau^2)))
lines(ml.n, max.dist = 1)
lines(reml.n, lwd = 2, max.dist = 1)
lines(ols.n, lty = 2, max.dist = 1)
lines(wls.n, lty =2, lwd = 2, max.dist = 1)
legend(0.2, .04, legend=c("ML","REML","OLS","WLS"), lty=c(1,1,2,2), lwd=c(1,2,1,2), cex=0.7)





trendtype <- "cte"
## Fitting models with nugget fixed to zero
ml <- likfit(houstongeo, ini = c(5,.5), fix.nugget = TRUE, trend = trendtype)
reml <- likfit(houstongeo, ini = c(5,.5), fix.nugget = TRUE, method = "RML", trend = trendtype)
ols <- variofit(bin1, ini = c(5,.5), fix.nugget = TRUE, weights="equal")
wls <- variofit(bin1, ini = c(5,.5), fix.nugget = TRUE)
	 
## Fitting models estimated nugget
ml.n <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = FALSE, nugget = 1, trend = trendtype)
reml.n <- likfit(houstongeo, ini = c(1,0.5), fix.nugget = FALSE, nugget = 1,
                 method = "RML", trend = trendtype)
ols.n <- variofit(bin1, ini = c(1,0.5), fix.nugget = FALSE, weights="equal")
wls.n <- variofit(bin1, ini = c(1,0.5), fix.nugget = FALSE)


# Now, plotting fitted models against empirical variogram
par(mfrow = c(2,1))
plot(bin1, main = expression(paste("fixed ", tau^2 == 0)))
lines(ml, max.dist = 1)
lines(reml, lwd = 2, max.dist = 1)
lines(ols, lty = 2, max.dist = 1)
lines(wls, lty = 2, lwd = 2, max.dist = 1)
legend(0.2, 15, legend=c("ML","REML","OLS","WLS"),lty=c(1,1,2,2),lwd=c(1,2,1,2), cex=0.7)
plot(bin1, main = expression(paste("estimated  ", tau^2)))
lines(ml.n, max.dist = 1)
lines(reml.n, lwd = 2, max.dist = 1)
lines(ols.n, lty = 2, max.dist = 1)
lines(wls.n, lty =2, lwd = 2, max.dist = 1)
legend(0.2, 15, legend=c("ML","REML","OLS","WLS"), lty=c(1,1,2,2), lwd=c(1,2,1,2), cex=0.7)








par(mfrow = c(2, 2))
points(houstongeo, xlab = "Coord X", ylab = "Coord Y")
points(houstongeo, xlab = "Coord X", ylab = "Coord Y",
       pt.divide = "rank.prop")
points(houstongeo, xlab = "Coord X", ylab = "Coord Y",
       cex.max = 1.7, col = gray(seq(1, 0.1, l = 100)),
       pt.divide = "equal")
points(houstongeo, pt.divide = "quintile", xlab = "Coord X",
       ylab = "Coord Y")

cloud1 <- variog(s100, option = "cloud", max.dist = 1)
cloud2 <- variog(s100, option = "cloud", estimator.type = "modulus",
                 max.dist = 1)
bin1 <- variog(s100, uvec = seq(0, 1, l = 11))
bin2 <- variog(s100, uvec = seq(0, 1, l = 11),
               estimator.type = "modulus")

par(mfrow = c(2, 2))
plot(cloud1, main = "classical estimator")
plot(cloud2, main = "modulus estimator")
plot(bin1, main = "classical estimator")
plot(bin2, main = "modulus estimator")


trendtype <- "1st"
cloud1 <- variog(s100, option = "cloud", max.dist = 1, trend=~x1 + x2)
cloud2 <- variog(s100, option = "cloud", estimator.type = "modulus",
                 max.dist = 1, trend=trendtype)
bin1 <- variog(s100, uvec = seq(0, 1, l = 11), trend=trendtype)
bin2 <- variog(s100, uvec = seq(0, 1, l = 11),
               estimator.type = "modulus", trend=trendtype)
par(mfrow = c(2, 2))
plot(cloud1, main = "classical estimator")
plot(cloud2, main = "modulus estimator")
plot(bin1, main = "classical estimator")
plot(bin2, main = "modulus estimator")

