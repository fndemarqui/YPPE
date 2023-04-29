
#---------------------------------------------
pope_fit <- function(time, status, Z, X, n_int, rho, tau,
                     hyper_parms, survreg, approach, hessian, ...) {

  n <- length(time)
  p <- ncol(X)
  q <- ncol(Z)
  idt <- as.numeric(cut(time, rho, include.lowest = TRUE))
  ttt <- TTT(time, rho)

  hyper_parms$mu_psi=0
  hyper_parms$sigma_psi=4
  hyper_parms$mu_phi=0
  hyper_parms$sigma_phi=4

  approach <- ifelse(approach=="mle", 0, 1)
  stan_data <- list(status=status, Z=Z, X=X, p=p, q=q, n=n, m=n_int,
                    tau=tau, ttt=ttt, approach=approach, idt=idt,
                    h1_gamma=hyper_parms$h1_gamma,
                    h2_gamma=hyper_parms$h2_gamma,
                    mu_psi=hyper_parms$mu_psi,
                    mu_phi=hyper_parms$sigma_phi,
                    sigma_psi=hyper_parms$sigma_psi,
                    sigma_phi=hyper_parms$sigma_psi,
                    mu_beta=hyper_parms$mu_beta,
                    sigma_beta=hyper_parms$sigma_beta,
                    survreg = survreg)

  if(approach == 0){
    fit <- rstan::optimizing(stanmodels$yppe,data=stan_data,
                             hessian=hessian, ...)
  }else{
    fit <- rstan::sampling(stanmodels$yppe, data=stan_data, ...)
  }

  return(fit)
}

#---------------------------------------------

#' pope: Fit Proportional Odds Regression Model with Piecewise Exponential baseline distribution.
#' @aliases{pope}
#' @export
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which pope is called.
#' @param n_int number of intervals of the PE distribution. If NULL, default value (square root of n) is used.
#' @param rho the time grid of the PE distribution. If NULL, the function timeGrid is used to compute rho.
#' @param tau the maximum time of follow-up. If NULL, tau = max(time), where time is the vector of observed survival times.
#' @param approach approach to be used to fit the model (mle: maximum likelihood; bayes: Bayesian approach).
#' @param hessian logical; If TRUE (default), the hessian matrix is returned when approach="mle".
#' @param hyper_parms a list containing the hyper-parameters of the prior distributions (when approach = "bayes"). If not specified, default values are used.
#' @param ... Arguments passed to either `rstan::optimizing` or `rstan::sampling` .
#' @return pope returns an object of class "pope" containing the fitted model.
#' @examples
#' \donttest{
#' # ML approach:
#' library(YPPE)
#' mle <- pope(Surv(time, status)~arm, data=ipass, n_int=10, approach="mle")
#' summary(mle)
#'
#' # Bayesian approach:
#' bayes <- pope(Surv(time, status)~arm, data=ipass, n_int=10, approach="bayes")
#' summary(bayes)
#' }
#'
pope <- function(formula, data, n_int=NULL, rho=NULL, tau=NULL, hessian=TRUE,
                 approach = c("mle", "bayes"),
                 hyper_parms = list(h1_gamma=0, h2_gamma=4,
                                    mu_beta=0, sigma_beta=4), ...){

  approach <- match.arg(approach)
  formula <- Formula::Formula(formula)
  mf <- stats::model.frame(formula=formula, data=data)
  Terms <- stats::terms(mf)
  resp <- stats::model.response(mf)
  time <- resp[,1]
  status <- resp[,2]
  X <- stats::model.matrix(formula, data = mf, rhs = 1)
  labels <- colnames(X)[-1]
  X <- matrix(X[,-1], ncol=length(labels))
  Z <- array(0, dim = c(0, 0))

  n <- length(time)
  q <- ncol(Z)
  p <- ncol(X)


  if(p > 0){
    survreg <- 4
  }else{
    survreg <- 0
  }


  if(is.null(tau)){
    tau <- max(time)
  }
  time <- time/tau

  if(is.null(rho)){
    rho <- timeGrid(time, status, n_int)
  }else{
    rho <- rho/tau
  }
  n_int <- length(rho) - 1

  fit <- pope_fit(time=time, status=status, Z=Z, X=X, n_int=n_int,
                  rho=rho, tau=tau, hyper_parms=hyper_parms,
                  survreg = survreg, approach = approach, hessian = hessian, ...)

  output <- list(fit=fit)

  output$n <- n
  output$p <- p

  output$n_int <- n_int
  output$rho <- rho*tau
  output$tau <- tau
  output$call <- match.call()
  output$formula <- formula
  output$terms <- stats::terms.formula(formula)
  output$mf <- mf
  output$labels <- labels
  output$approach <- approach
  output$survreg <- "po"
  output$npar <- p+n_int

  class(output) <- c("pope","yppe")
  return(output)
}


