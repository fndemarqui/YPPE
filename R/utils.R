



#---------------------------------------------
#' Generic S3 method vcov
#' @aliases vcov
#' @export
#' @param object a fitted model object
#' @param ... further arguments passed to or from other methods.
#' @return the variance-covariance matrix associated the regression coefficients.
#'
vcov <- function(object, ...) UseMethod("vcov")

#---------------------------------------------
#' Covariance of the regression coefficients
#'
#' @aliases vcov.yppe
#' @rdname vcov-methods
#' @method vcov yppe
#' @export
#' @export vcov
#' @param object an object of the class yppe
#' @param ... further arguments passed to or from other methods.
#' @return  the variance-covariance matrix associated with the regression coefficients.
#'
#'
vcov.yppe <- function(object, ...){
  p <- object$p
  q <- object$q

  if(p==0){
    std <- with(object, rep(z_sd, 2))
  }else{
    std <- with(object, c(rep(z_sd, 2), x_sd))
  }

  std <- diag(1/std)
  V <- MASS::ginv(-object$fit$hessian)[1:(2*q+p), 1:(2*q+p)]
  colnames(V) <- names(object$fit$par)[1:(2*q+p)]
  rownames(V) <- names(object$fit$par)[1:(2*q+p)]
  V <- std%*%V%*%std
  #class(V) <- "vcov.ypbp"
  return(V)
}

#---------------------------------------------
#' Generic S3 method coef
#' @aliases coef.yppe
#' @export
#' @param object a fitted model object
#' @param ... further arguments passed to or from other methods.
#' @return the estimated regression coefficients
#'
coef <- function(object, ...) UseMethod("coef")

#---------------------------------------------
#' Estimated regression coefficients
#'
#' @aliases coef.yppe
#' @rdname coef-methods
#' @method coef yppe
#' @export
#' @export coef
#' @param object an object of the class yppe
#' @param ... further arguments passed to or from other methods
#' @return  the estimated regression coefficients
#'
#'
coef.yppe <- function(object, ...){
  p <- object$p
  q <- object$q
  coeffs <- object$fit$par[1:(2*q+p)]
  return(coeffs)
}


#---------------------------------------------
#' Generic S3 method confint
#' @aliases confint.yppe
#' @export
#' @param object a fitted model object
#' @param ... further arguments passed to or from other methods.
#' @return the estimated regression coefficients
#'
confint <- function(object, ...) UseMethod("confint")

#---------------------------------------------
#' Confidence intervals for the regression coefficients
#'
#' @aliases confint.yppe
#' @rdname confint-methods
#' @method confint yppe
#' @export
#' @export confint
#' @param object an object of the class yppe
#' @param level the confidence level required
#' @param ... further arguments passed to or from other methods
#' @return  100(1-alpha) confidence intervals for the regression coefficients
#'
#'
confint.yppe <- function(object, level=0.95, ...){
  p <- object$p
  q <- object$q
  V <- vcov(object)
  par.hat <- object$fit$par[1:(2*q+p)]
  alpha <- 1-level
  d <- qnorm(1 - alpha/2)*sqrt(diag(V))
  lower <- par.hat - d
  upper <- par.hat + d
  CI <- cbind(lower, upper)
  labels <- round(100*(c(alpha/2, 1-alpha/2)),1)
  colnames(CI) <- paste0(labels, "%")
  return(CI)
}
