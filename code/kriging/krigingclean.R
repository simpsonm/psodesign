library(dplyr)

houston8 <- read.csv("houston8-16eighthouravgs.csv")
for(i in 1:31){
  houston8[,1 + i] <- as.character(houston8[,1 + i])
}
houston8[houston8 == "NV"] <- NA
for(i in 1:31){
  houston8[,1 + i] <- as.numeric(houston8[,1 + i])
}
houston8$nvalid <- apply(!is.na(as.matrix(houston8[,2:32])), 1, sum)
houston8$avg <- apply((as.matrix((houston8[,2:32]))), 1, mean, na.rm=TRUE)
houston8$sd <- apply((as.matrix((houston8[,2:32]))), 1, sd, na.rm=TRUE)
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

houstonout <- filter(houstonout, avg > -Inf)

write.csv(houstonout, file = "houstonout.csv", row.names=FALSE)
