library(geoR)
library(plyr)
library(dplyr)
library(lubridate)

houston8 <- read.csv("houston8-16eighthouravgs.csv")
for(i in 1:31){
  houston8[,1 + i] <- as.character(houston8[,1 + i])
}
houston8[houston8 == "NV"] <- NA
for(i in 1:31){
  houston8[,1 + i] <- as.numeric(houston8[,1 + i])
}
houston8$nvalid <- apply(!is.na(as.matrix(houston8[,2:32])), 1, sum)
houston8$avg <- apply((as.matrix(houston8[,2:32])), 1, mean, na.rm=TRUE)
houston8$sd <- apply((as.matrix(houston8[,2:32])), 1, sd, na.rm=TRUE)
houston8$se <- houston8$sd/sqrt(houston8$nvalid)

houstongis <- read.table("houstongis.csv", TRUE, "\t")
houstongis <- distinct(houstongis, Name, .keep_all=TRUE)

houston8$Name <- unlist(strsplit(as.character(houston8$Site), " C[0-9]"))[2*(0:43) + 1]

houstonsmol <- select(houston8, Name, avg, sd, nvalid, se)


matches <- inner_join(houstonsmol, houstongis)
misses <- anti_join(houstonsmol, houstongis)
candidates <- anti_join(houstongis, houstonsmol)
misses$Name[grep("St.", misses$Name)] <- paste(unlist(strsplit(misses$Name[grep("St.", misses$Name)], "St.")), "Street", sep = "")
misses$Name[grep("Co.", misses$Name)] <- paste(unlist(strsplit(misses$Name[grep("Co.", misses$Name)], "Co.")), "County", sep = "")
misses$Name[grep("HRM-3", misses$Name)] <- "HRM 3 Haden Rd"
misses$Name[grep("U of H", misses$Name)] <- "UH Sugarland"
misses$Name[grep("Hou.Deer", misses$Name)] <- "Houston Deer Park 2"

matches2 <- inner_join(misses, candidates)
misses2 <- anti_join(misses, candidates)

houstonout <- select(full_join(matches, matches2), Name, AQS, Latitude, Longitude, avg, sd, nvalid, se)

houston <- read.csv("houstonfull.csv", skip = 10)
houston$Date <- ymd(houston$Date)
houston$ID <- paste(houston$County.Cd, houston$Site.ID, sep = "")
houston$dateID <- paste(houston$Date, houston$ID, sep = "--")

houstonlen <- summarise(group_by(houston, dateID), len = length(Value))
target <- filter(houstonlen, len >= 24)$dateID
offtarget <- filter(houstonlen, len < 24)$dateID
houstonsub <- filter(houston, dateID %in% target)

houstonavgs <- ddply(houstonsub, .(ID, Date, dateID),
                     function(x)stats::filter(x$Value, rep(1/8, 8),
                                              sides = 1, method = "conv")[-c(1:7)])

houstonmaxes <- houstonavgs %>% group_by(dateID) %>%
  transmute(ID = ID, date = Date, month = month(Date),
            n = sum(!is.na(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11,
                             V12, V13, V14, V15, V16, V17))),
            max = max(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11,
                        V12, V13, V14, V15, V16, V17),
                      na.rm=TRUE))


houstonaug <- read.csv("houstonaugfull.csv", skip = 10)
houstonaug$Date <- ymd(houstonaug$Date)
houstonaug$ID <- paste(houstonaug$County.Cd, houstonaug$Site.ID, sep = "")
houstonaug$dateID <- paste(houstonaug$Date, houstonaug$ID, sep = "--")

houstonauglen <- summarise(group_by(houstonaug, dateID), len = length(Value))
target <- filter(houstonauglen, len >= 24)$dateID
offtarget <- filter(houstonauglen, len < 24)$dateID
houstonaugsub <- filter(houstonaug, dateID %in% target)

houstonaugavgs <- ddply(houstonaugsub, .(ID, Date, dateID),
                     function(x)stats::filter(x$Value, rep(1/8, 8),
                                              sides = 1, method = "conv")[-c(1:7)])

houstonaugmaxes <- houstonaugavgs %>% group_by(dateID) %>%
  transmute(ID = ID, date = Date, month = month(Date),
            n = sum(!is.na(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11,
                             V12, V13, V14, V15, V16, V17))),
            max = max(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11,
                        V12, V13, V14, V15, V16, V17),
                      na.rm=TRUE))

houstonaugmonths <- houstonaugmaxes %>% group_by(ID) %>%
  summarise(mean(max[max>-Inf]), sd(max[max>-Inf]), n = sum(max>-Inf))


houstonmaxes2 <- houstonsub %>% group_by(dateID) %>%
  transmute(ID = ID, date = Date, month = month(Date),
            n = sum(!is.na(Value)),
            max = max(Value, na.rm=TRUE))
           
houstonmonths <- houstonmaxes %>% group_by(ID, month) %>%
  summarise(mean(max), sd(max))

houstonmonths2 <- houstonmaxes2 %>% filter(n >=20) %>% group_by(ID, month) %>%
  summarise(mean(max), sd(max))

ttt <- filter(houstonmonths2, month == 6)

ttt <- filter(houstonmonths, month == 8)


houstonavgs$max <- apply(houstonavgs[,3 + 1:17], 1, max, na.rm=TRUE)
houstonavgs$n <- apply(houstonavgs[,3 + 1:17], 1, function(x)sum(!is.na(x)))



filtest <- filter(values, rep(1/8, 8), sides = 1)
realtest <- rep(0, length(filtest) - 7)
for(i in 1:(length(filtest) - 7)){
  realtest[i] <- mean(values[1:8+i-1])
}

movavg <- function(values, times){
  diffs <- (as.numeric(times))
  ts(values
}

dailyhigh <- ddply(houston, .(Meth.Cd, County.Cd, Site.ID, Date), summarize, high = max(Value), nval = length(Value))

subset(houston, Meth.Cd == 87 & County.Cd == 167 & Site.ID == 1034 & Date == "2016-08-02")

dailyhigh$month <- month(dailyhigh$Date)

monthlyavghigh <- ddply(dailyhigh, .(Meth.Cd, County.Cd, Site.ID, month), summarize,
                        avghigh = mean(high))

subset(monthlyavghigh, month == 6)

subset(dailyhigh, month == 8 & Site.ID == 1034 & County.Cd == 201)

unique(subset(monthlyavghigh, month == 8)$Site.ID)
unique(subset(monthlyavghigh, month == 7)$Site.ID)
unique(subset(monthlyavghigh, month == 6)$Site.ID)

subset(dailyhigh, month == 8)

test <- ddply(houston, .(Meth.Cd, County.Cd, Site.ID), summarize, numpoc = lenunique(POC), avgpoc = mean(POC), nobs = length(Value))

lenunique <- function(x) length(unique(x))

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
