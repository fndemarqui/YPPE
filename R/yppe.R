

#---------------------------------------------
yppe.mle <- function(time, status, Z, n_int, rho, tau, hessian, ...) {

  n <- nrow(Z)
  q <- ncol(Z)
  idt <- as.numeric(cut(time, rho, include.lowest = TRUE))
  ttt <- matrix(nrow=n, ncol=n_int)
  for(i in 1:n){
    for(j in 1:n_int){
      ttt[i,j] <- (min(time[i], rho[j+1]) - rho[j])*(time[i] - rho[j]>0)
    }
  }

  hyper_parms = list(mu_lambda=0, sigma_lambda=2,
                     mu_psi=0, sigma_psi=3,
                     mu_phi=0, sigma_phi=3,
                     mu_beta=0, sigma_beta=4)

  stan_data <- list(status=status, Z=Z, q=q, n=n, m=n_int,
                    tau=tau, ttt=ttt, approach=0, idt=idt,
                    mu_lambda=hyper_parms$mu_lambda,
                    sigma_lambda=hyper_parms$sigma_lambda,
                    mu_psi=hyper_parms$mu_psi,
                    mu_phi=hyper_parms$sigma_phi,
                    sigma_psi=hyper_parms$sigma_psi,
                    sigma_phi=hyper_parms$sigma_psi)

  fit <- rstan::optimizing(stanmodels$yppe,data=stan_data,
                           hessian=hessian, verbose=FALSE, ...)

  fit$par <- fit$par[-grep("loglik", names(fit$par))]
  fit$theta_tilde <- fit$theta_tilde[-grep("loglik", names(fit$theta_tilde))]

  return(fit)
}

#---------------------------------------------
yppe.bayes <- function(time, status, Z, n_int, rho, tau, hyper_parms, ...) {

  n <- nrow(Z)
  q <- ncol(Z)
  idt <- as.numeric(cut(time, rho, include.lowest = TRUE))
  ttt <- matrix(nrow=n, ncol=n_int)
  for(i in 1:n){
    for(j in 1:n_int){
      ttt[i,j] <- (min(time[i], rho[j+1]) - rho[j])*(time[i] - rho[j]>0)
    }
  }
  stan_data <- list(status=status, Z=Z, q=q, n=n, m=n_int,
                    tau=tau, ttt=ttt, approach=1, idt=idt,
                    mu_lambda=hyper_parms$mu_lambda,
                    sigma_lambda=hyper_parms$sigma_lambda,
                    mu_psi=hyper_parms$mu_psi,
                    mu_phi=hyper_parms$sigma_phi,
                    sigma_psi=hyper_parms$sigma_psi,
                    sigma_phi=hyper_parms$sigma_psi)

  fit <- rstan::sampling(stanmodels$yppe, data=stan_data, ...)

  return(fit)
}

#---------------------------------------------
yppe2.mle <- function(time, status, Z, X, n_int, rho, tau, hessian, ...) {

  n <- nrow(Z)
  q <- ncol(Z)
  p <- ncol(X)
  idt <- as.numeric(cut(time, rho, include.lowest = TRUE))
  ttt <- matrix(nrow=n, ncol=n_int)
  for(i in 1:n){
    for(j in 1:n_int){
      ttt[i,j] <- (min(time[i], rho[j+1]) - rho[j])*(time[i] - rho[j]>0)
    }
  }

  hyper_parms = list(mu_lambda=0, sigma_lambda=4,
                     mu_psi=0, sigma_psi=4,
                     mu_phi=0, sigma_phi=4,
                     mu_beta=0, sigma_beta=4)

  stan_data <- list(status=status, Z=Z, X=X, p=p, q=q, n=n, m=n_int,
                    tau=tau, ttt=ttt, approach=0, idt=idt,
                    mu_lambda=hyper_parms$mu_lambda,
                    sigma_lambda=hyper_parms$sigma_lambda,
                    mu_psi=hyper_parms$mu_psi,
                    mu_phi=hyper_parms$sigma_phi,
                    sigma_psi=hyper_parms$sigma_psi,
                    sigma_phi=hyper_parms$sigma_psi,
                    mu_beta=hyper_parms$mu_beta,
                    sigma_beta=hyper_parms$sigma_beta)

  fit <- rstan::optimizing(stanmodels$yppe2,data=stan_data,
                           hessian=hessian, verbose=FALSE, ...)

  fit$par <- fit$par[-grep("loglik", names(fit$par))]
  fit$theta_tilde <- fit$theta_tilde[-grep("loglik", names(fit$theta_tilde))]

  return(fit)
}

#---------------------------------------------
yppe2.bayes <- function(time, status, Z, X, n_int, rho, tau, hyper_parms, ...) {

  n <- nrow(Z)
  q <- ncol(Z)
  p <- ncol(X)
  idt <- as.numeric(cut(time, rho, include.lowest = TRUE))
  ttt <- matrix(nrow=n, ncol=n_int)
  for(i in 1:n){
    for(j in 1:n_int){
      ttt[i,j] <- (min(time[i], rho[j+1]) - rho[j])*(time[i] - rho[j]>0)
    }
  }

  stan_data <- list(status=status, Z=Z, X=X, p=p, q=q, n=n, m=n_int,
                    tau=tau, ttt=ttt, approach=1, idt=idt,
                    mu_lambda=hyper_parms$mu_lambda,
                    sigma_lambda=hyper_parms$sigma_lambda,
                    mu_psi=hyper_parms$mu_psi,
                    mu_phi=hyper_parms$sigma_phi,
                    sigma_psi=hyper_parms$sigma_psi,
                    sigma_phi=hyper_parms$sigma_psi,
                    mu_beta=hyper_parms$mu_beta,
                    sigma_beta=hyper_parms$sigma_beta)

  fit <- rstan::sampling(stanmodels$yppe2, data=stan_data,
                         verbose=FALSE, ...)

  return(fit)
}

#---------------------------------------------

#' Fits the Yang and Prentice model with baseline distribution modelled by the piecewise exponential distribution.
#' @aliases{yppe}
#' @export
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which yppe is called.
#' @param n_int number of intervals of the PE distribution. If NULL, default value (square root of n) is used.
#' @param rho the time grid of the PE distribution. If NULL, the function timeGrid is used to compute rho.
#' @param tau the maximum time of follow-up. If NULL, tau = max(time), where time is the vector of observed survival times.
#' @param approach approach to be used to fit the model (mle: maximum likelihood; bayes: Bayesian approach).
#' @param hessian logical; If TRUE (default), the hessian matrix is returned when approach="mle".
#' @param hyper_parms a list containing the hyper-parameters of the prior distributions (when approach = "bayes"). If not specified, default values are used.
#' @param ... Arguments passed to either `rstan::optimizing` or `rstan::sampling` .
#' @return yppe returns an object of class "yppe" containing the fitted model.
#' @examples
#' \donttest{
#' # ML approach:
#' library(YPPE)
#' mle <- yppe(Surv(time, status)~arm, data=ipass, approach="mle")
#' summary(mle)
#'
#' # Bayesian approach:
#' bayes <- yppe(Surv(time, status)~arm, data=ipass, approach="bayes")
#' summary(bayes)
#' }
#'
yppe <- function(formula, data, n_int=NULL, rho=NULL, tau=NULL, hessian=TRUE,
                 approach = c("mle", "bayes"),
                 hyper_parms = list(mu_lambda=0, sigma_lambda=4,
                                    mu_psi=0, sigma_psi=4,
                                    mu_phi=0, sigma_phi=4,
                                    mu_beta=0, sigma_beta=4), ...){

  approach <- match.arg(approach)
  formula <- Formula::Formula(formula)
  mf <- stats::model.frame(formula=formula, data=data)
  Terms <- stats::terms(mf)
  resp <- stats::model.response(mf)
  time <- resp[,1]
  status <- resp[,2]
  Z <- stats::model.matrix(formula, data = mf, rhs = 1)
  X <- suppressWarnings(try( stats::model.matrix(formula, data = mf, rhs = 2), TRUE))
  labels <- colnames(Z)[-1]
  labels.ph <- colnames(X)[-1]
  Z <- matrix(Z[,-1], ncol=length(labels))
  if(ncol(X)>0){
    labels.ph <- colnames(X)[-1]
    X <- matrix(X[,-1], ncol=length(labels.ph))
  }

  n <- nrow(Z)
  q <- ncol(Z)
  p <- ncol(X)


  if(is.null(tau)){
    tau <- max(time)
  }
  time <- time/tau

  if(is.null(n_int)){
    n_int <- ceiling(sqrt(length(time)))
  }

  if(is.null(rho)){
    rho <- timeGrid(time, status, n_int)
  }

  if(approach=="mle"){
    if(p==0){
      fit <- yppe.mle(time=time, status=status, Z=Z, n_int=n_int, rho=rho,
                      tau=tau, hessian=hessian, ...)
    }else{
      fit <- yppe2.mle(time=time, status=status, Z=Z, X=X, n_int=n_int, rho=rho,
                       tau=tau, hessian=hessian, ...)
    }
  }else{
    if(p==0){
      fit <- yppe.bayes(time=time, status=status, Z=Z, n_int=n_int, rho=rho,
                        tau=tau, hyper_parms=hyper_parms, ...)
    }else{
      fit <- yppe2.bayes(time=time, status=status, Z=Z, X=X, n_int=n_int, rho=rho,
                         tau=tau, hyper_parms=hyper_parms, ...)
    }
  }

  output <- list(fit=fit)

  output$n <- n
  output$q <- q
  output$p <- p

  output$n_int <- n_int
  output$rho <- rho
  output$tau <- tau
  output$call <- match.call()
  output$formula <- formula
  output$terms <- stats::terms.formula(formula)
  output$mf <- mf
  output$labels <- labels
  output$approach <- approach

  if(p>0){
    output$labels.ph <- labels.ph
  }


  class(output) <- "yppe"
  return(output)
}


#---------------------------------------------

yppeBoot <- function(formula, data, n_int=NULL, rho=NULL, tau=NULL,
                     nboot = 4000, ...){
  formula <- Formula::Formula(formula)
  mf <- stats::model.frame(formula=formula, data=data)
  Terms <- stats::terms(mf)
  resp <- stats::model.response(mf)
  time <- resp[,1]
  status <- resp[,2]
  Z <- stats::model.matrix(formula, data = mf, rhs = 1)
  X <- suppressWarnings(try( stats::model.matrix(formula, data = mf, rhs = 2), TRUE))
  labels <- colnames(Z)[-1]
  labels.ph <- colnames(X)[-1]
  Z <- matrix(Z[,-1], ncol=length(labels))
  if(ncol(X)>0){
    labels.ph <- colnames(X)[-1]
    X <- matrix(X[,-1], ncol=length(labels.ph))
  }

  n <- nrow(Z)
  q <- ncol(Z)
  p <- ncol(X)


  if(is.null(tau)){
    tau <- max(time)
  }
  time <- time/tau

  if(is.null(n_int)){
    n_int <- ceiling(sqrt(length(time)))
  }

  if(is.null(rho)){
    rho <- timeGrid(time, status, n_int)
  }

  index <- 1:n
  index1 <- which(status==1)
  index2 <- which(status==0)
  n1 <- length(index1)
  n2 <- length(index2)
  par <- matrix(nrow=nboot, ncol=(2*q+p+n_int))

  for(step in 1:nboot){
    samp1 <- sample(index1, size=n1, replace=TRUE)
    samp2 <- sample(index2, size=n2, replace=TRUE)
    samp <- c(samp1, samp2)
    suppressWarnings({invisible(utils::capture.output(object <- yppe(formula, data=data[samp,], n_int=n_int, rho=rho, tau=tau, hessian=FALSE, approach="mle", init=0)))})
    if(class(object)!="try-error"){
      par[step, ] <- object$fit$par[-grep("log_", names(object$fit$par))]
      step <- step + 1
    }
  }

  colnames(par) <- names(object$fit$par[-grep("log_", names(object$fit$par))])

  return(par)
}
