source("../psofun.R")
source("testfuns.R")

nswarm <- 40
np <- 20
niter <- 1000
nrep <- 50
inertias <- c(0.7298, 1/(log(2)*2))
cognitives <- c(1.496, log(2) + 1/2)
socials <- c(1.496, log(2) + 1/2)
lower <- rep(-100, np)
upper <- rep(100, np)
nnbors <- c(1, 3, 40)
alpha <- 0.2*niter
beta <- 2
rates <- c(0.3, 0.5)
ccc <- 0.1
df <- 1
pcuts <- c(0, 0.5)
sig0 <- 1
inertia0 <- 1.2

time <- 0:niter

normvec <- function(x){
  sqrt(mean(x^2))
}

set.seed(342)

psoout <- NULL
cat("\n")
for(repl in 1:nrep){
  for(nnbor in nnbors){
    for(CF in c(TRUE, FALSE)){
      for(parset in 1:2){
        for(fnum in c(1:8)[-c(3,6)]){
          cat("rep = "); cat(repl); cat("; fnum = "); cat(fnum); cat("; nnbor = "); cat(nnbor)
          cat("; CF = "); cat(CF); cat("; parset = "); cat(parset); cat("\n")
          for(style in c("AT1", "AT2")){
            rate <- ifelse(style=="AT1", rates[1], rates[2])
            pcut <- pcuts[parset]
            temp <- sbbpso(niter, nswarm, nnbor, sig0,
                           pcut = pcut, CF = CF, AT = TRUE, rate = rate, df = df, ccc = 0.1,
                           obj = negfwrap, opt = fnum)
            algid <- paste("BBPSO", parset, ifelse(CF, "CF", "notCF"), style, sep = "-")
            tempdat <- data.frame(obj = fnum, logpost = temp[["values"]],
                                  argnorm = apply(temp[["pars"]], 2, normvec),
                                  time = time, algid = algid,
                                  type = "BBPSO", parset = parset, CF = CF,
                                  style = style, nbhd = nnbor, rep = repl,
                                  inertias = temp$sigs)
            psoout <- rbind(psoout, tempdat)
          }
          for(style in c("CI", "DI", "AT1", "AT2")){
            rate <- ifelse(style=="AT1", rates[1], rates[2])
            c.in <- ifelse(style=="CI", inertias[parset], inertia0)
            c.co <- cognitives[parset]
            c.so <- socials[parset]
            temp <- spso(niter, nswarm, nnbor, c.in, c.co, c.so, negfwrap, lower, upper,
                         style = substr(style, 1, 2), CF = CF, alpha = alpha, beta = beta,
                         rate = rate, ccc = ccc, opt = fnum)
            algid <- paste("PSO", parset, ifelse(CF, "CF", "notCF"), style, sep = "-")
            tempdat <- data.frame(obj = fnum, logpost = temp[["values"]],
                                  argnorm = apply(temp[["pars"]], 2, normvec),
                                  time = time, algid = algid,
                                  type = "PSO", parset = parset, CF = CF,
                                  style = style, nbhd = nnbor, rep = repl,
                                  inertias = temp$inertias)
            psoout <- rbind(psoout, tempdat)
          }
          write.csv(psoout, file = "psosimsout.csv", row.names=FALSE)
        }
      }
    }
  }
}
