

#---------------------------------------------
# Hazard function of the PE distribution

hpexp <- function(x, tgrid, rates)
{
  id <- as.numeric(cut(x, tgrid, include.lowest = TRUE))
  return(rates[id])
}


#---------------------------------------------
# Cumulative hazard function of the PE distribution

Hpexp <- function(x, tgrid, rates)
{
  h <- diff(tgrid)
  n <- length(x)
  id <- as.numeric(cut(x, tgrid, include.lowest = TRUE))
  area <- cumsum(h*rates)

  if(n==1)
  {
    if(id==1){
      Ht <- x*rates[1]
    }else{
      Ht <- (x-tgrid[id])*rates[id] + area[id-1]
    }
  }else{
    Ht  <- rep(0, n)
    Ht[which(id==1)] <- x[which(id==1)]*rates[1]
    aux <- which(id>1)
    Ht[aux] <- (x[aux]-tgrid[id[aux]])*rates[id[aux]]+ area[id[aux]-1]
  }
  return(Ht)
}


#---------------------------------------------

ppexp <- function(x, tgrid, rates, lower.tail=TRUE, p.log=FALSE)
{
  Ht <- Hpexp(x, tgrid, rates)
  St <- exp(-Ht)

  if(lower.tail==TRUE)
  {
    if(p.log==FALSE)
      return(1-St)
    else
      return(log(1-St))
  }else{
    if(p.log==FALSE)
      return(St)
    else
      return(log(St))
  }
}


##########################################################################################################
# Density function of the PE distribution
dpexp <- function(x, tgrid, rates, log=FALSE)
{
  Ht <- Hpexp(x, tgrid, rates)
  id <- as.numeric(cut(x, tgrid, include.lowest = TRUE))
  log_ft <- log(rates[id]) - Ht
  if(log==FALSE)
    return(exp(log_ft))
  else
    return(log_ft)
}


##########################################################################################################
# returns the quantile of the piecewise exponential distribution:
qpexp <- function(p, tgrid, rates, lower.tail=TRUE, p.log=FALSE)
{
  if(lower.tail==TRUE)
  {
    u <- p
  }else{
    u <- 1-p
  }

  h <- diff(tgrid)
  n <- length(p)
  area <- c(0,cumsum(h*rates))
  Ft <- 1-exp(-area)
  id <- as.numeric(cut(u, Ft, include.lowest = TRUE))

  if(n==1)
  {
    if(id==1){
      q <- -log(1-u)/rates[1]
    }else{
      q <- tgrid[id] - (area[id] + log(1-u))/rates[id]
    }
  }else{
    q  <- rep(0, n)
    q[which(id==1)] <- -log(1-u[which(id==1)])/rates[1]
    aux <- which(id>1)
    q[aux] <- tgrid[id[aux]] - (area[id[aux]] + log(1-u[aux]))/rates[id[aux]]
  }
  return(q)
}

##########################################################################################################
# returns a sample of size n from the piecewise exponential distribution:

rpexp <- function(n, tgrid, rates){
  u <- stats::runif(n)
  x <- qpexp(u, tgrid, rates)
  return(x)
}

#----------------------------------

#' Time grid
#' @export
#' @param time Vector of failure times
#' @param status Vector of failure indicators
#' @param n.int Optional. Number of intervals. If \code{NULL}, the number of intervals is
#' set to be equal to the number of distinct observed failure times.
#' @return Time grid.


timeGrid <- function(time, status, n.int=NULL)
{
  o <- order(time)
  time <- time[o]
  status <- status[o]
  time.aux <- unique(time[status==1])
  if(is.null(n.int))
  {
    n.int <- length(time.aux)
  }

  m <- length(time.aux)
  if(n.int > m)
  {
    tgrid <- c(0,unique(time[status==1]))
    tgrid[length(tgrid)] <- Inf
  }
  else
  {
    b <- min(m,n.int)
    k1 <- trunc(m/b)
    r <- m-b*k1
    k2 <- k1+1
    idf1 <- seq(k1,(b-r)*k1, k1)
    idf2 <- sort(seq(m,max(idf1),-k2))
    idf <- unique(c(idf1,idf2))
    tgrid <- c(0,time.aux[idf])
    tgrid[length(tgrid)] <- Inf
  }
  return(tgrid)
}


findInt <- function(time, tgrid)
{
  id <- as.numeric(cut(time, tgrid, include.lowest = TRUE))
  return(id)
}
