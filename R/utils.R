
TTT <- function(time, rho){
  n <- length(time)
  n_int <- length(rho) - 1
  ttt <- matrix(nrow = n, ncol = n_int)
  for(i in 1:n){
    for(j in 1:n_int){
      ttt[i,j] <- (min(time[i], rho[j+1]) - rho[j])*(time[i] - rho[j]>0)
    }
  }
  return(ttt)
}



#---------------------------------------------
#' Variance-covariance matrix for a yppe model
#'
#' @aliases vcov.ypppe
#' @description This function extracts and returns the variance-covariance matrix associated with the regression coefficients when the maximum likelihood estimation approach is used in the model fitting.
#' @export
#' @param object an object of the class yppe.
#' @param ... further arguments passed to or from other methods.
#' @return  the variance-covariance matrix associated with the regression coefficients.
#'
#'
vcov.yppe <- function(object, ...){
  p <- object$p
  q <- object$q
  V <- MASS::ginv(-object$fit$hessian)[1:(2*q+p), 1:(2*q+p)]
  colnames(V) <- names(object$fit$par)[1:(2*q+p)]
  rownames(V) <- names(object$fit$par)[1:(2*q+p)]
  #class(V) <- "vcov.yppe"
  return(V)
}


#---------------------------------------------
#' Estimated regression coefficients
#'
#' @aliases coef.yppe
#' @description This function returns the estimated regression coefficients when the maximum likelihood estimation approach is used in the model fitting.
#' @export
#' @param object an object of the class yppe.
#' @param ... further arguments passed to or from other methods.
#' @return  the estimated regression coefficients.
#' @examples
#' \dontrun{
#' fit <- yppe(Surv(time, status)~arm, data=ipass, init = 0)
#' coef(fit)
#' }
#'
coef.yppe <- function(object, ...){
  p <- object$p
  q <- object$q
  coeffs <- object$fit$par[1:(2*q+p)]
  return(coeffs)
}



#---------------------------------------------
#' Confidence intervals for the regression coefficients
#'
#' @aliases confint.yppe
#' @export
#' @param object an object of the class yppe.
#' @param level the confidence level required.
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param ... further arguments passed to or from other methods.
#' @return  100(1-alpha) confidence intervals for the regression coefficients.
#'
#'
confint.yppe <- function(object, parm = NULL, level=0.95, ...){
  p <- object$p
  q <- object$q
  V <- vcov(object)
  par.hat <- object$fit$par[1:(2*q+p)]
  alpha <- 1-level
  d <- stats::qnorm(1 - alpha/2)*sqrt(diag(V))
  lower <- par.hat - d
  upper <- par.hat + d
  CI <- cbind(lower, upper)
  labels <- round(100*(c(alpha/2, 1-alpha/2)),1)
  colnames(CI) <- paste0(labels, "%")
  return(CI)
}


#---------------------------------------------
#' Model.matrix method for yppe models
#'
#' @aliases model.matrix.yppe
#' @description Reconstruct the model matrix (or matrices if the alternative formulation of the YP model is used) for a yppe model.
#' @export
#' @param object an object of the class yppe.
#' @param ... further arguments passed to or from other methods.
#' @return  The model matrix (or matrices) for the fit.
#' @examples
#' \dontrun{
#' fit <- yppe(Surv(time, status)~arm, data=ipass)
#' model.matrix(fit)
#'}
#'
model.matrix.yppe <- function(object, ...){
  formula <- Formula::Formula(object$formula)
  mf <- object$mf
  Z <- stats::model.matrix(formula, data = mf, rhs = 1)
  X <- suppressWarnings(try( stats::model.matrix(formula, data = mf, rhs = 2), TRUE))
  attrZ <- attributes(Z)
  attrZ$dim[2] <- ncol(Z) - 1
  attrZ$assign <- attrZ$assign[-1]
  attrZ$dimnames[[2]] <- attrZ$dimnames[[2]][-1]
  Z <- matrix(Z[,-1], ncol=ncol(Z) - 1)
  attributes(Z) <- attrZ
  if(ncol(X)>0){
    attrX <- attributes(X)
    attrX$dim[2] <- ncol(X) - 1
    attrX$assign <- attrX$assign[-1]
    attrX$dimnames[[2]] <- attrX$dimnames[[2]][-1]
    X <- matrix(X[,-1], ncol=ncol(X) - 1)
    attributes(X) <- attrX
  }
  if(ncol(X)>0){
    out <- list(Z=Z, X=X)
  }else{
    out <- Z
  }

  return(out)
}


# internal function used to compute AIC (several models)
get_arg_names <- function(...) {
  argnames <- sys.call()
  paste0(lapply(argnames[-1], as.character))
}

#' Akaike information criterion
#' @aliases AIC.yppe
#' @param object an object of the class yppe.
#' @param ... further arguments passed to or from other methods.
#' @param k numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.
#' @return the Akaike information criterion
#' @export
#'

AIC.yppe <- function(object, ..., k = 2){
  objects <- c(as.list(environment()), list(...))
  argnames <- sys.call()
  argnames <- paste0(lapply(argnames[-1], as.character))
  k <- objects[[2]]
  objects <- objects[-2]
  J <- nargs()
  aic <- c()
  for(j in 1:J){
    loglik <- objects[[j]]$fit$value
    npar <- objects[[j]]$npar
    aic[j] <- -2*loglik + k*npar
  }
  if(length(argnames)>1){
    names(aic) <- argnames
    aic <- sort(aic)
  }
  return(aic)
}


#' Extract Log-Likelihood
#' @aliases logLik.yppe
#' @param object an object of the class yppe.
#' @param ... further arguments passed to or from other methods.
#' @return the log-likelihood associated with the fitted model.
#' @export
#'
logLik.yppe <- function(object, ...){
  objects <- c(as.list(environment()), list(...))
  argnames <- sys.call()
  argnames <- paste0(lapply(argnames[-1], as.character))
  J <- nargs()
  loglik <- c()
  for(j in 1:J){
    loglik[j] <- objects[[j]]$fit$value
  }
  if(length(argnames)>1){
    names(loglik) <- argnames
    loglik <- sort(loglik)
  }
  return(loglik)
}


#---------------------------------------------
#' Generic S3 method rates
#' @export
#' @param object a fitted model object.
#' @details Method only available for ML approach.
#' @param ... further arguments passed to or from other methods.
#' @return the estimated failure rates for the PE distribution.
#'
rates <- function(object, ...) UseMethod("rates")

#' Estimated failure rates for the PE distribution
#' @aliases rates.phpe
#' @export
#' @details Method only available for ML approach.
#' @param object a fitted model object.
#' @param ... further arguments passed to or from other methods.
#' @return the estimated failure rates for the PE distribution.
#'
rates.phpe <- function(object, ...){
  if(object$approach != "mle"){
    warning("Method only available for ML approach")
  }else{
    p <- object$p
    rho <- object$rho
    m <- length(rho) - 1
    tau <- object$tau
    par <- object$fit$par
    rates <- par[(p+1):(p+m)]/tau
    return(rates)
  }
}

#' Estimated failure rates for the PE distribution
#' @aliases rates.pope
#' @export
#' @details Method only available for ML approach.
#' @param object a fitted model object.
#' @param ... further arguments passed to or from other methods.
#' @return the estimated failure rates for the PE distribution.
#'
rates.pope <- function(object, ...){
  if(object$approach != "mle"){
    warning("Method only available for ML approach")
  }else{
    p <- object$p
    rho <- object$rho
    m <- length(rho) - 1
    tau <- object$tau
    par <- object$fit$par
    rates <- par[(p+1):(p+m)]/tau
    return(rates)
  }

}

#' Estimated failure rates for the PE distribution
#' @aliases rates.yppe
#' @export
#' @details Method only available for ML approach.
#' @param object a fitted model object.
#' @param ... further arguments passed to or from other methods.
#' @return the estimated failure rates for the PE distribution.
#'
rates.yppe <- function(object, ...){
  if(object$approach != "mle"){
    warning("Method only available for ML approach")
  }else{
    p <- object$p
    q <- object$q
    rho <- object$rho
    m <- length(rho) - 1
    tau <- object$tau
    par <- object$fit$par
    rates <- par[(2*q+p+1):(2*q+p+m)]/tau
    names(rates) <- NULL
    return(rates)
  }
}

