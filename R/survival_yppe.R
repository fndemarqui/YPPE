

yppeSurv <- function(time, z, par, rho, tau, n_int){
  q <- length(z)
  psi <- par[1:q]
  phi <- par[(q+1):(2*q)]
  lambda <- par[(2*q+1):(2*q+n_int)]/tau
  rho <- rho*tau
  theta_S <- exp( as.numeric(z%*%psi) )
  theta_L <- exp( as.numeric(z%*%phi) )
  Ht0 <- Hpexp(time, rho, lambda)
  St0 <- exp(-Ht0)
  Ft0 <- 1-St0
  St <- exp( -theta_L*(log(theta_L*St0 + theta_S*Ft0)-(as.numeric(z%*%phi) - Ht0)  ) )
  class(St) <- "survfit.yppe"
  return(St)
}


yppeSurv2 <- function(time, z, x, par, rho, tau, n_int){
  q <- length(z)
  p <- length(x)
  psi <- par[1:q]
  phi <- par[(q+1):(2*q)]
  beta <- par[(2*q+1):(2*q+p)]
  lambda <- par[(2*q+p+1):(2*q+p+n_int)]/tau
  rho <- rho*tau
  theta_S <- exp( as.numeric(z%*%psi) )
  theta_L <- exp( as.numeric(z%*%phi) )
  theta_C <- exp( as.numeric(x%*%beta) )
  Ht0 <- Hpexp(time, rho, lambda)
  St0 <- exp(-Ht0)
  Ft0 <- 1-St0
  St <- exp( -theta_L*theta_C*(log(theta_L*St0 + theta_S*Ft0)-(as.numeric(z%*%phi) - Ht0)  ) )
  class(St) <- "survfit.yppe"
  return(St)
}

#---------------------------------------------
#' Generic S3 method survfit
#' @aliases survfit
#' @export
#' @param object a fitted model object
#' @param ... further arguments passed to or from other methods.
#' @return the crossing survival time
#'
survfit <- function(object, ...) UseMethod("survfit")

#---------------------------------------------
#' Survival function for the YPPE model
#'
#' @aliases survfit.yppe
#' @rdname survfit-methods
#' @method survfit yppe
#' @export
#' @export survfit
#' @param object an object of the class yppe
#' @param newdata a data frame containing the set of explanatory variables.
#' @param ... further arguments passed to or from other methods.
#' @return  a list containing the estimated survival probabilities.
#' @examples
#' \donttest{
#' # ML approach:
#' library(YPPE)
#' mle <- yppe(Surv(time, status)~arm, data=ipass, approach="mle")
#' summary(mle)
#' ekm <- survfit(Surv(time, status)~arm, data=ipass)
#' newdata <- data.frame(arm=0:1)
#' St <- survfit(mle, newdata)
#' time <- sort(ipass$time)
#' plot(ekm, col=1:2)
#' lines(time, St[[1]])
#' lines(time, St[[2]], col=2)
#'
#' # Bayesian approach:
#' bayes <- yppe(Surv(time, status)~arm, data=ipass, approach="bayes")
#' summary(bayes)
#' ekm <- survfit(Surv(time, status)~arm, data=ipass)
#' newdata <- data.frame(arm=0:1)
#' St <- survfit(bayes, newdata)
#' time <- sort(ipass$time)
#' plot(ekm, col=1:2)
#' lines(time, St[[1]])
#' lines(time, St[[2]], col=2)
#' }
#'
survfit.yppe <- function(object, newdata, ...){
  mf <- object$mf
  labels <- names(mf)[-1]
  time <- sort( stats::model.response(mf)[,1])
  status <- sort( stats::model.response(mf)[,2])
  data <- data.frame(cbind(time, status, mf[,-1]))
  names(data) <- c("time", "status", names(mf)[-1])
  rho <- object$rho
  n_int <- object$n_int
  tau <- object$tau
  labels <- match.arg(names(newdata), labels, several.ok = TRUE)
  formula <- object$formula
  Z <- as.matrix(stats::model.matrix(formula, data = newdata, rhs = 1)[,-1])
  X <- suppressWarnings(try( as.matrix(stats::model.matrix(formula, data = newdata, rhs = 2)[,-1]), TRUE))
  St <- list()


  if(object$approach=="mle"){
    par <- object$fit$par[-grep("log_", names(object$fit$par))]
    if(object$p==0){
      for(i in 1:nrow(newdata)){
        St[[i]] <- yppeSurv(time, Z[i,], par, rho, tau, n_int)
      }
    }else{
      for(i in 1:nrow(newdata)){
        St[[i]] <- yppeSurv2(time, Z[i,], X[i,], par, rho, tau, n_int)
      }
    }
  }else{ # Bayesian approach
    samp <- rstan::extract(object$fit)
    if(object$p==0){
      par <- cbind(samp$psi, samp$phi, samp$lambda)
      for(i in 1:nrow(newdata)){
        aux <- apply(par, 1, yppeSurv, time=time, z=Z[i,], rho=rho, tau=tau, n_int=n_int)
        St[[i]] <- apply(aux, 1, mean)
      }
    }else{
      par <- cbind(samp$psi, samp$phi, samp$beta, samp$lambda)
      for(i in 1:nrow(newdata)){
        aux <- apply(par, 1, yppeSurv2, time=time, z=Z[i,], x=X[i,], rho=rho, tau=tau, n_int=n_int)
        St[[i]] <- apply(aux, 1, mean)
      }
    }
  }


  class(St) <- "survfit.yppe"
  return(St)
}




diffSurv <- function(time, z1, z2, par, rho, tau, n_int){
    St1 <- yppeSurv(time=time, z=z1, par=par, rho=rho, tau=tau, n_int=n_int)
    St2 <- yppeSurv(time=time, z=z2, par=par, rho=rho, tau=tau, n_int=n_int)
  return(St1-St2)
}

diffSurv2 <- function(time, z1, z2, x, par, rho, tau, n_int){
    St1 <- yppeSurv2(time=time, z=z1, x=x, par=par, rho=rho, tau=tau, n_int=n_int)
    St2 <- yppeSurv2(time=time, z=z2, x=x, par=par, rho=rho, tau=tau, n_int=n_int)
  return(St1-St2)
}



yppeCrossSurv <- function(z1, z2, par, rho, tau0, tau, n_int){
  I <- c(tau0, 1.5*tau)
  t <- try(stats::uniroot(diffSurv, interval=I, z1=z1, z2=z2, par=par,
                          rho=rho, tau=tau, n_int=n_int)$root, TRUE)
  if(class(t)=="try-error"){
    return(NA)
  }else{
    return(t)
  }
}

yppeCrossSurv2 <- function(z1, z2, x, par, rho, tau0, tau, n_int){
  I <- c(tau0, 1.5*tau)
  t <- try(stats::uniroot(diffSurv2, interval=I, z1=z1, z2=z2, x=x, par=par,
                          rho=rho, tau=tau, n_int=n_int)$root, TRUE)
  if(class(t)=="try-error"){
    return(NA)
  }else{
    return(t)
  }
}


#---------------------------------------------
#' Generic S3 method crossTime
#' @aliases crossTime
#' @export
#' @param object a fitted model object
#' @param ... further arguments passed to or from other methods.
#' @return the crossing survival time
#'
crossTime <- function(object, ...) UseMethod("crossTime")


#' Computes the crossing survival times
#'
#' @aliases crossTime.yppe
#' @rdname crossTime-methods
#' @method crossTime yppe
#' @export
#' @export crossTime
#' @param object an object of class yppe
#' @param newdata1 a data frame containing the first set of explanatory variables
#' @param newdata2 a data frame containing the second set of explanatory variables
#' @param conf.level level of the confidence/credible intervals
#' @param nboot number of bootstrap samples (default nboot=4000); ignored if approach="bayes".
#' @param ... further arguments passed to or from other methods.
#' @return  the crossing survival time
#' @examples
#' \donttest{
#' # ML approach:
#' library(YPPE)
#' mle <- yppe(Surv(time, status)~arm, data=ipass, approach="mle")
#' summary(mle)
#' newdata1 <- data.frame(arm=0)
#' newdata2 <- data.frame(arm=1)
#' tcross <- crossTime(mle, newdata1, newdata2)
#' tcross
#' ekm <- survfit(Surv(time, status)~arm, data=ipass)
#' newdata <- data.frame(arm=0:1)
#' St <- survfit(mle, newdata)
#' time <- sort(ipass$time)
#' plot(ekm, col=1:2)
#' lines(time, St[[1]])
#' lines(time, St[[2]], col=2)
#' abline(v=tcross, col="blue")
#'
#' # Bayesian approach:
#' bayes <- yppe(Surv(time, status)~arm, data=ipass, approach="bayes")
#' summary(bayes)
#' newdata1 <- data.frame(arm=0)
#' newdata2 <- data.frame(arm=1)
#' tcross <- crossTime(bayes, newdata1, newdata2)
#' tcross
#' ekm <- survfit(Surv(time, status)~arm, data=ipass)
#' newdata <- data.frame(arm=0:1)
#' St <- survfit(bayes, newdata)
#' time <- sort(ipass$time)
#' plot(ekm, col=1:2)
#' lines(time, St[[1]])
#' lines(time, St[[2]], col=2)
#' abline(v=tcross, col="blue")
#' }
#'
crossTime.yppe <- function(object, newdata1, newdata2,
                           conf.level=0.95, nboot=4000, ...){
  q <-object$q
  p <-object$p
  mf <- object$mf
  labels <- names(mf)[-1]
  time <- stats::model.response(mf)[,1]
  status <- stats::model.response(mf)[,2]
  o <- order(time)
  data <- data.frame(cbind(time, status, mf[,-1]))[o,]
  names(data) <- c("time", "status", labels)
  tau0 <- min(time)
  rho <- object$rho
  n_int <- object$n_int
  tau <- object$tau
  labels <- match.arg(names(newdata1), names(newdata2), several.ok=TRUE)
  labels <- match.arg(names(mf)[-1], names(newdata1), several.ok=TRUE)
  z1 <- matrix(stats::model.matrix(object$formula, data = newdata1, rhs = 1)[,-1], ncol=q)
  z2 <- matrix(stats::model.matrix(object$formula, data = newdata2, rhs = 1)[,-1], ncol=q)
  if(p>0){
    x <- matrix(stats::model.matrix(object$formula, data = newdata2, rhs = 2)[,-1], ncol=p)
  }

  I <- c(tau0, 1.5*tau)
  alpha <- 1 - conf.level
  prob <- c(alpha/2, 1-alpha/2)

  if(object$approach=="mle"){
    t <- c()
    par <- object$fit$par[-grep("log_", names(object$fit$par))]
    for(i in 1:nrow(newdata1)){
      if(p==0){
        t[i] <- yppeCrossSurv(z1=z1[i,], z2=z2[i,], par=par, rho=rho, tau0=tau0, tau=tau, n_int=n_int)
      }else{
        t[i,] <- yppeCrossSurv2(z1=z1[i,], z2=z2[i,], x=x[i,], par=par, rho=rho, tau0=tau0, tau=tau, n_int=n_int)
      }
    }
    par <- with(object, yppeBoot(formula=formula, data=data, n_int=n_int,
                                rho=rho, tau=tau, nboot=nboot, prob=prob))

    ci <- matrix(nrow=nrow(newdata1), ncol=2)
    if(object$p==0){
      for(i in 1:nrow(newdata1)){
        aux <- apply(par, 1, yppeCrossSurv, z1=z1[i,], z2=z2[i,], rho=rho, tau0=tau0, tau=tau, n_int=n_int)
        ci[i,] <- stats::quantile(aux, probs=prob, na.rm=TRUE)
      }
    }else{
      for(i in 1:nrow(newdata1)){
        aux <- apply(par, 1, yppeCrossSurv2, z1=z1[i,], z2=z2[i,], x=x[i,], rho=rho, tau0=tau0, tau=tau, n_int=n_int)
        ci[i,] <- stats::quantile(aux, probs=prob, na.rm=TRUE)
      }
    }

    t <- data.frame(cbind(t, ci))
    names(t) <- c("Est.", paste(100*prob, "%", sep=""))
  }else{ # Bayesian approach
    t <- matrix(nrow=nrow(newdata1), ncol=3)
    samp <- rstan::extract(object$fit)
    if(object$p==0){
      par <- cbind(samp$psi, samp$phi, samp$lambda)
      for(i in 1:nrow(newdata1)){
        aux <- apply(par, 1, yppeCrossSurv, z1=z1[i,], z2=z2[i,], rho=rho, tau0=tau0, tau=tau, n_int=n_int)
        ci <- stats::quantile(aux, probs=prob, na.rm=TRUE)
        t[i,] <- c(mean(aux, na.rm=TRUE), ci)
      }
    }else{
      par <- cbind(samp$psi, samp$phi, samp$beta, samp$lambda)
      for(i in 1:nrow(newdata1)){
        aux <- apply(par, 1, yppeCrossSurv2, z1=z1[i,], z2=z2[i,], x=x[i,], rho=rho, tau0=tau0, tau=tau, n_int=n_int)
        ci <- stats::quantile(aux, probs=prob, na.rm=TRUE)
        t[i,] <- c(mean(aux, na.rm=TRUE), ci)
      }
    }
    t <- as.data.frame(t)
    names(t) <- c("Est.", names(ci) )
  }
  #class(t) <- "crossTime.yppe"
  return(t)
}

