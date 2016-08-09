library(shapefiles)
library(spdep)
library(maptools)

popdat <- read.csv("ACS_14_5YR_S0101_with_ann.csv", header=TRUE, colClasses = "character")
popdat <- popdat[,1:4]
colnames(popdat)[3:4] <- c("County", "Population")
popdat$Population <- as.numeric(popdat$Population)

coshape <- readShapePoly("county/050_00.shp")
names(coshape)[1] <- "ID"

coshape@data <- merge(coshape@data, popdat, by.x = "ID", by.y = "Id")

A <- poly2nb(coshape)
adjmat <- nb2mat(A, zero.policy=TRUE, style="B")
adjmat <- adjmat[-2917, -2917]
## shape file includes bedford city, VA as separate from bedford county
## as of 2013 this is no longer the case.
## since bedford city is an island in bedford county, this will not affect nbhd structure

xmat <- matrix(1, nrow = nrow(popdat), ncol = 1)
pxmat <- xmat%*%tcrossprod(solve(crossprod(xmat)), xmat)
Imat <- diag(nrow(popdat))
G <- (Imat - pxmat)%*%adjmat%*%(Imat - pxmat)
eigG <- eigen(G)


R <- 100
smat <- eigG$vectors[,1:R]
dat <- cbind(coshape$Population, xmat, smat)
colnames(dat) <- c("pop", "x", paste("s", 1:R, sep = "."))

popdat <- dat

save(popdat, file = "popdat.RData")
