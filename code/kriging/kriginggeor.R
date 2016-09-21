library(geoR)

data(s100)

summary(s100)

plot(s100)

par(mfrow = c(2, 2))
points(s100, xlab = "Coord X", ylab = "Coord Y")
points(s100, xlab = "Coord X", ylab = "Coord Y",
       pt.divide = "rank.prop")
points(s100, xlab = "Coord X", ylab = "Coord Y",
       cex.max = 1.7, col = gray(seq(1, 0.1, l = 100)),
       pt.divide = "equal")
points(s100, pt.divide = "quintile", xlab = "Coord X",
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
