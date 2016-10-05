negfwrap <- function(x, opt){
  return(-fwrap(x, opt))
}

fwrap <- function(x, opt){
  if(opt == 1){
    out <- f1(x)
  } else if(opt ==2){
    out <- f2(x)
  } else if(opt ==3){
    out <- f3(x)
  } else if(opt ==4){
    out <- f4(x)
  } else if(opt ==5){
    out <- f5(x)
  } else if(opt ==6){
    out <- f6(x)
  } else if(opt ==7){
    out <- f7(x)
  } else if(opt ==8){
    out <- f8(x)
  }
  return(out)  
}

## De Jong
f1 <- function(x){
  if(is.null(dim(x))){
    out <- sum(x^2)
  } else {
    out <- apply(x^2, 2, sum)
  }
  return(-out)
}

## Schwefel 1.2
f2 <- function(x){
    if(is.null(dim(x))){
      out <- sum(cumsum(x)^2)
    } else {
      out <- apply(apply(x, 2, cumsum)^2, 2, sum)
    }
    return(-out)
}

## Schaffer's f6
f3 <- function(x){
  if(is.null(dim(x))){
    out <- 0.5 + ( sin( sqrt( sum(x^2) ) ) - 0.5 )/( 1 + 0.0001*( sum(x^2) ) )^2
  } else {
    x2 <- apply(x^2, 2, sum)
    out <- 0.5 + (sin(sqrt(x2)) - 0.5)/(1 + 0.0001*(x2))^2
  }
  return(-out)
}

## Rosenbrock
f4 <- function(x){
  x <- x - 1
  if(is.null(dim(x))){
    out <- sum( 100*( x[-1] - x[-length(x)]^2 )^2 + (x[-length(x)] - 1)^2)
  } else {
    xtemp <- 100*(x[-1,] - x[-nrow(x),]^2)^2 + (x[-nrow(x),] - 1)^2
    out <- apply(xtemp, 2, sum)
  }
  return(-out)
}

## Rastrigrin
f5 <- function(x){
  if(is.null(dim(x))){
    out <- sum(x^2 - cos(2*pi*x) + 10) - 9*length(x)
  } else {
    xtemp <- x^2 - cos(2*pi*x) + 10
    out <- apply(xtemp, 2, sum) - 9*nrow(x)
  }
  return(-out)  
}

## Schwefel 2.6
f6 <- function(x){
 if(is.null(dim(x))){
    out <- -sum(x*sin(sqrt(abs(x))))
  } else {
    xtemp <- -x*sin(sqrt(abs(x)))
    out <- apply(xtemp, 2, sum)
  }
  return(-out)  
}

## Griewank
f7 <- function(x){
 if(is.null(dim(x))){
    out <- sum(x^2)/4000 - prod(cos(x/sqrt(1:length(x)))) + 1
  } else {
    xtemp <- cos(x/sqrt(1:nrow(x)))
    out <- apply(x^2, 2, sum)/4000 - apply(xtemp, 2, prod) + 1
  }
  return(-out)  
}

## Ackley
f8 <- function(x){
  if(is.null(dim(x))){
    out <- -20*exp(-0.2*sqrt(mean(x^2))) - exp(mean(cos(2*pi*x))) + 20 + exp(1)
  } else {
    xtemp <- cos(2*pi*x)
    out <- -20*exp(-0.2*sqrt(apply(x^2, 2, mean))) - exp(apply(xtemp, 2, mean)) + 20 + exp(1)
  }
  return(-out)  
}

## Normal for stochastic PSO
fnorm <- function(x, partidx, MU, SIG){
  ndim <- ncol(partidx)
  muidx <- matrix(MU[partidx], ncol = ndim)
  xidx <- matrix(x[partidx], ncol = ndim)
  sigidx <- matrix(SIG[partidx], ncol = ndim)
  outmat <- dnorm(xidx, muidx, sigidx, TRUE)
  out <- apply(outmat, 1, sum)
  return(out)
}
