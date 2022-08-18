

yppeSurv <- function(time, z, par, rho, tau, n_int){
  q <- length(z)
  psi <- par[1:q]
  phi <- par[(q+1):(2*q)]
  gamma <- par[(2*q+1):(2*q+n_int)]/tau
  theta <- exp( as.numeric(z%*%phi) )
  ratio <- exp( as.numeric(z%*%(psi-phi)) )
  Ht0 <- Hpexp(time, rho, gamma)
  Rt0 = expm1(Ht0)
  aux <- ratio*Rt0
  St <- exp(-theta*log1p(aux))
  return(St)
}


yppeSurv2 <- function(time, z, x, par, rho, tau, n_int){
  q <- length(z)
  p <- length(x)
  psi <- par[1:q]
  phi <- par[(q+1):(2*q)]
  beta <- par[(2*q+1):(2*q+p)]
  gamma <- par[(2*q+p+1):(2*q+p+n_int)]/tau
  theta <- exp( as.numeric(z%*%phi) + as.numeric(x%*%beta) )
  ratio <- exp( as.numeric(z%*%(psi-phi)) )
  Ht0 <- Hpexp(time, rho, gamma)
  Rt0 = expm1(Ht0)
  aux <- ratio*Rt0
  St <- exp(-theta*log1p(aux))
  class(St) <- "survfit.yppe"
  return(St)
}

#---------------------------------------------
#' survfit method for yppe models
#'
#' @aliases survfit.yppe
#' @description Computes the predicted survivor function for a yppe model.
#' @importFrom survival survfit
#' @export
#' @param formula an object of the class yppe
#' @param newdata a data frame containing the set of explanatory variables.
#' @param ... further arguments passed to or from other methods.
#' @return  a list containing the estimated survival probabilities.
#' @examples
#' \donttest{
#' # ML approach:
#' library(YPPE)
#' mle <- yppe(Surv(time, status)~arm, data=ipass, n_int=10, approach="mle")
#' summary(mle)
#' ekm <- survival::survfit(Surv(time, status)~arm, data=ipass)
#' newdata <- data.frame(arm=0:1)
#' St <- survfit(mle, newdata)
#' plot(ekm, col=1:2)
#' with(St, lines(time, surv[[1]]))
#' with(St, lines(time, surv[[2]], col=2))
#'
#' # Bayesian approach:
#' bayes <- yppe(Surv(time, status)~arm, data=ipass, n_int=10, approach="bayes")
#' summary(bayes)
#' ekm <- survival::survfit(Surv(time, status)~arm, data=ipass)
#' newdata <- data.frame(arm=0:1)
#' St <- survfit(bayes, newdata)
#' plot(ekm, col=1:2)
#' with(St, lines(time, surv[[1]]))
#' with(St, lines(time, surv[[2]], col=2))
#' }
#'
survfit.yppe <- function(formula, newdata, ...){
  object <- formula
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
  Z <- stats::model.matrix(formula, data = newdata, rhs = 1)[,-1, drop = FALSE]
  X <- suppressWarnings(try( stats::model.matrix(formula, data = newdata, rhs = 2)[,-1, drop = FALSE], TRUE))
  St <- list()


  if(object$approach=="mle"){
    #par <- object$fit$par[-grep("log_", names(object$fit$par))]
    par <- object$fit$par
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
      par <- cbind(samp$psi, samp$phi, samp$gamma)
      for(i in 1:nrow(newdata)){
        aux <- apply(par, 1, yppeSurv, time=time, z=Z[i,], rho=rho, tau=tau, n_int=n_int)
        St[[i]] <- apply(aux, 1, mean)
      }
    }else{
      par <- cbind(samp$psi, samp$phi, samp$beta, samp$gamma)
      for(i in 1:nrow(newdata)){
        aux <- apply(par, 1, yppeSurv2, time=time, z=Z[i,], x=X[i,], rho=rho, tau=tau, n_int=n_int)
        St[[i]] <- apply(aux, 1, mean)
      }
    }
  }


  out <- list(time = time, surv = St)
  class(out) <- "survfit.yppe"
  return(out)
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
  if(is(t, "try-error")){
    return(NA)
  }else{
    return(t)
  }
}

yppeCrossSurv2 <- function(z1, z2, x, par, rho, tau0, tau, n_int){
  I <- c(tau0, 1.5*tau)
  t <- try(stats::uniroot(diffSurv2, interval=I, z1=z1, z2=z2, x=x, par=par,
                          rho=rho, tau=tau, n_int=n_int)$root, TRUE)
  if(is(t, "try-error")){
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
#' @param nboot number of bootstrap samples (default nboot=1000); ignored if approach="bayes".
#' @param ... further arguments passed to or from other methods.
#' @return  the crossing survival time
#' @examples
#' \donttest{
#' # ML approach:
#' library(YPPE)
#' mle <- yppe(Surv(time, status)~arm, data=ipass, n_int=10, approach="mle", init = 0)
#' summary(mle)
#' newdata1 <- data.frame(arm=0)
#' newdata2 <- data.frame(arm=1)
#' tcross <- crossTime(mle, newdata1, newdata2, nboot = 10)
#' tcross
#' ekm <- survfit(Surv(time, status)~arm, data=ipass)
#' newdata <- data.frame(arm=0:1)
#' St <- survfit(mle, newdata)
#' plot(ekm, col=1:2)
#' with(St, lines(time, surv[[1]]))
#' with(St, lines(time, surv[[2]], col=2))
#' abline(v=tcross, col="blue")
#'
#' # Bayesian approach:
#' bayes <- yppe(Surv(time, status)~arm, data=ipass, n_int=10, approach="bayes", chains=1, iter=10)
#' summary(bayes)
#' newdata1 <- data.frame(arm=0)
#' newdata2 <- data.frame(arm=1)
#' tcross <- crossTime(bayes, newdata1, newdata2)
#' tcross
#' ekm <- survfit(Surv(time, status)~arm, data=ipass)
#' newdata <- data.frame(arm=0:1)
#' St <- survfit(bayes, newdata)
#' plot(ekm, col=1:2)
#' with(St, lines(time, surv[[1]]))
#' with(St, lines(time, surv[[2]], col=2))
#' abline(v=tcross, col="blue")
#' }
#'
crossTime.yppe <- function(object, newdata1, newdata2,
                           conf.level=0.95, nboot=1000, ...){
  q <-object$q
  p <-object$p
  mf <- object$mf
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
  z1 <- stats::model.matrix(object$formula, data = newdata1, rhs = 1)[,-1, drop = FALSE]
  z2 <- stats::model.matrix(object$formula, data = newdata2, rhs = 1)[,-1, drop = FALSE]
  if(p>0){
    x <- stats::model.matrix(object$formula, data = newdata2, rhs = 2)[,-1, drop = FALSE]
  }

  I <- c(tau0, 1.5*tau)
  alpha <- 1 - conf.level
  prob <- c(alpha/2, 1-alpha/2)

  if(object$approach=="mle"){
    t <- c()
    #par <- object$fit$par[-grep("log_", names(object$fit$par))]
    par <- object$fit$par
    for(i in 1:nrow(newdata1)){
      if(p==0){
        t[i] <- yppeCrossSurv(z1=z1[i,], z2=z2[i,], par=par, rho=rho, tau0=tau0, tau=tau, n_int=n_int)
      }else{
        t[i] <- yppeCrossSurv2(z1=z1[i,], z2=z2[i,], x=x[i,], par=par, rho=rho, tau0=tau0, tau=tau, n_int=n_int)
      }
    }
    par <- yppe_boot(object, nboot=nboot)

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
      par <- cbind(samp$psi, samp$phi, samp$gamma)
      for(i in 1:nrow(newdata1)){
        aux <- apply(par, 1, yppeCrossSurv, z1=z1[i,], z2=z2[i,], rho=rho, tau0=tau0, tau=tau, n_int=n_int)
        ci <- stats::quantile(aux, probs=prob, na.rm=TRUE)
        t[i,] <- c(mean(aux, na.rm=TRUE), ci)
      }
    }else{
      par <- cbind(samp$psi, samp$phi, samp$beta, samp$gamma)
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

