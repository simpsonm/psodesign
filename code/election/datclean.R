# Load in data for region indicators
# Use "state", an R data file (type ?state from the R command window for info)
#
# Regions:  1=northeast, 2=south, 3=north central, 4=west, 5=d.c.
# We have to insert d.c. (it is the 9th "state" in alphabetical order)

data(state)                  # "state" is an R data file
state.abbr <- c (state.abb[1:8], "DC", state.abb[9:50])
dc <- 9
not.dc <- c(1:8,10:51)
region <- c(3,4,4,3,4,4,1,1,5,3,3,4,4,2,2,2,2,3,3,1,1,1,2,2,3,2,4,2,4,1,1,4,1,3,2,2,3,4,1,1,3,2,3,3,4,1,3,4,1,2,4)

# Load in data from the CBS polls in 1988

library (foreign)
polls <- read.dta("fromwebsite/polls.dta")

# Select just the data from the last survey (#9158)

table(polls$survey)                # look at the survey id's


polls.subset <- polls    # select the subset of interest

polldat <- polls[!is.na(polls.subset$bush),3:10]
ok <- polldat$survey==9158            # define the condition
polldatsmall <- polldat[ok,]

presvote <- read.dta("fromwebsite/presvote.dta")
v.prev <- presvote$g76_84pr
not.dc <- c(1:8,10:51)
candidate.effects <- read.table ("fromwebsite/candidate_effects.dat", header=T)
v.prev[not.dc] <- v.prev[not.dc] +
 (candidate.effects$X76 + candidate.effects$X80 + candidate.effects$X84)/3

statedat <- data.frame(prev = v.prev, region = region, abbr = state.abbr)

nage <- length(unique(polldat$age))
nedu <- length(unique(polldat$edu))
nageedu <- nage*nedu
nstate <- max(unique(polldat$state))
npoll <- length(unique(polldat$survey))

xmat <- model.matrix(~female + black + female*black, data = polldat)
statemat <- matrix(0, nrow(polldat), nstate)
agemat <- matrix(0, nrow(polldat), nage)
edumat <- matrix(0, nrow(polldat), nedu)
ageedumat <- matrix(0, nrow(polldat), nageedu)
pollmat <- matrix(0, nrow(polldat), npoll)
colnames(statemat) <- rep("", ncol(statemat))
colnames(ageedumat) <- rep("", ncol(ageedumat))
colnames(pollmat) <- rep("", ncol(pollmat))
colnames(agemat) <- rep("", ncol(agemat))
colnames(edumat) <- rep("", ncol(edumat))
for(i in 1:nstate){
  statemat[,i] <- as.numeric(polldat$state == i)
  colnames(statemat)[i] <- paste("state", i, sep = ".")
}
for(i in 1:nage){
  agemat[,i] <- as.numeric(polldat$age == i)
  colnames(agemat)[i] <- paste("age", i, sep = ".")
}
for(i in 1:nedu){
  edumat[,i] <- as.numeric(polldat$edu == i)
  colnames(edumat)[i] <- paste("edu", i, sep = ".")
}
ageedumat <- matrix(0, ncol = 16, nrow = nrow(polldat))
colnames(ageedumat) <- rep("", 16)
for(i in 1:nage){
  for(j in 1:nedu){
    ageedumat[,(j - 1)*nage + i] <- as.numeric(polldat$age == i & polldat$edu == j)
    colnames(ageedumat)[(j - 1)*nage + i] <- paste("age", "edu", i, j, sep = ".")
  }
}
for(i in 1:npoll){
  pollmat[,i] <- as.numeric(polldat$survey == unique(polldat$survey)[i])
  colnames(pollmat)[i] <- paste("poll", unique(polldat$survey)[i], sep = ".")
}

nregion <- length(unique(statedat$region))
regionmat <- matrix(0, ncol = nregion, nrow = nstate)
for(i in 1:nregion){
  regionmat[,i] <- as.numeric(statedat$region == i)
}
colnames(regionmat) <- c(paste("region", 1:nregion, sep = "."))

wmat <- cbind(v.prev, regionmat)
xzmat <- cbind(xmat, statemat, ageedumat)

datlistplus <- list(y = polldat$bush, statedat = statedat, xmat = xmat,
                    statemat = statemat, ageedumat = ageedumat,
                    regionmat = regionmat, prev = v.prev, agemat = agemat,
                    edumat = edumat, pollmat = pollmat, wmat = wmat, xzmat = xzmat,
                    betamn = 0, betavar = 1000, sig2a = 1, sig2b = 1,
                    nobs = length(polldat$bush), nstate = nstate,
                    nage = nage, nedu = nedu, nageedu = nageedu,
                    npoll = npoll, nregion = nregion, nbeta = 4)

xmat <- model.matrix(~female + black + female*black, data = polldatsmall)
statemat <- matrix(0, nrow(polldatsmall), nstate)
agemat <- matrix(0, nrow(polldatsmall), nage)
edumat <- matrix(0, nrow(polldatsmall), nedu)
ageedumat <- matrix(0, nrow(polldatsmall), nageedu)
pollmat <- matrix(0, nrow(polldatsmall), npoll)
colnames(statemat) <- rep("", ncol(statemat))
colnames(ageedumat) <- rep("", ncol(ageedumat))
colnames(pollmat) <- rep("", ncol(pollmat))
colnames(agemat) <- rep("", ncol(agemat))
colnames(edumat) <- rep("", ncol(edumat))
for(i in 1:nstate){
  statemat[,i] <- as.numeric(polldatsmall$state == i)
  colnames(statemat)[i] <- paste("state", i, sep = ".")
}
for(i in 1:nage){
  agemat[,i] <- as.numeric(polldatsmall$age == i)
  colnames(agemat)[i] <- paste("age", i, sep = ".")
}
for(i in 1:nedu){
  edumat[,i] <- as.numeric(polldatsmall$edu == i)
  colnames(edumat)[i] <- paste("edu", i, sep = ".")
}
ageedumat <- matrix(0, ncol = 16, nrow = nrow(polldatsmall))
colnames(ageedumat) <- rep("", 16)
for(i in 1:nage){
  for(j in 1:nedu){
    ageedumat[,(j - 1)*nage + i] <- as.numeric(polldatsmall$age == i & polldatsmall$edu == j)
    colnames(ageedumat)[(j - 1)*nage + i] <- paste("age", "edu", i, j, sep = ".")
  }
}
for(i in 1:npoll){
  pollmat[,i] <- as.numeric(polldatsmall$survey == unique(polldatsmall$survey)[i])
  colnames(pollmat)[i] <- paste("poll", unique(polldatsmall$survey)[i], sep = ".")
}

datlistsmall <- list(y = polldatsmall$bush, statedat = statedat, xmat = xmat,
                    statemat = statemat, ageedumat = ageedumat,
                    regionmat = regionmat, prev = v.prev, agemat = agemat,
                    edumat = edumat, wmat = wmat, xzmat = xzmat,
                    betamn = 0, betavar = 1000, sig2a = 1, sig2b = 1,
                    nobs = length(polldatsmall$bush), nstate = nstate,
                    nage = nage, nedu = nedu, nageedu = nageedu,
                    nregion = nregion, nbeta = 4)

save(datlistplus,  file = "datlistplus.RData")
save(datlistsmall, file = "datlistsmall.RData")
