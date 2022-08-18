

popeSurv <- function(time, x, par, rho, tau, n_int){
  p <- length(x)
  beta <- par[1:p]
  gamma <- par[(p+1):(p+n_int)]/tau
  Ht0 <- Hpexp(time, rho, gamma)
  lambda <- exp( as.numeric(x%*%beta) )
  Rt0 = expm1(Ht0)
  aux <- lambda*Rt0
  St <- exp(-log1p(aux))
  return(St)
}


#---------------------------------------------
#' survfit method for pope models
#'
#' @aliases survfit.pope
#' @description Computes the predicted survivor function for a pope model.
#' @importFrom survival survfit
#' @export
#' @param formula an object of the class pope
#' @param newdata a data frame containing the set of explanatory variables.
#' @param ... further arguments passed to or from other methods.
#' @return  a list containing the estimated survival probabilities.
#' @examples
#' \donttest{
#' # ML approach:
#' library(YPPE)
#' mle <- pope(Surv(time, status)~arm, data=ipass, n_int=10, approach="mle", init = 0)
#' summary(mle)
#' ekm <- survival::survfit(Surv(time, status)~arm, data=ipass)
#' newdata <- data.frame(arm=0:1)
#' St <- survfit(mle, newdata)
#' plot(ekm, col=1:2)
#' with(St, lines(time, surv[[1]]))
#' with(St, lines(time, surv[[2]], col=2))
#'
#' # Bayesian approach:
#' bayes <- pope(Surv(time, status)~arm, data=ipass, n_int=10, approach="bayes")
#' summary(bayes)
#' ekm <- survival::survfit(Surv(time, status)~arm, data=ipass)
#' newdata <- data.frame(arm=0:1)
#' St <- survfit(bayes, newdata)
#' plot(ekm, col=1:2)
#' with(St, lines(time, surv[[1]]))
#' with(St, lines(time, surv[[2]], col=2))
#' }
#'
survfit.pope <- function(formula, newdata, ...){
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
  X <- stats::model.matrix(formula, data = newdata)[,-1, drop = FALSE]
  St <- list()


  if(object$approach=="mle"){
    par <- object$fit$par
    for(i in 1:nrow(newdata)){
      St[[i]] <- popeSurv(time, X[i,], par, rho, tau, n_int)
    }
  }else{ # Bayesian approach
    samp <- rstan::extract(object$fit)
    par <- cbind(samp$beta, samp$gamma)
    for(i in 1:nrow(newdata)){
      aux <- apply(par, 1, popeSurv, time=time, x=X[i,], rho=rho, tau=tau, n_int=n_int)
      St[[i]] <- apply(aux, 1, mean)
    }
  }


  out <- list(time = time, surv = St)
  class(out) <- "survfit.pope"
  return(out)
}


