


rYP <- function(u, Z, psi, phi, gamma, baseline){
  # psi = phi: PH
  # phi = 0: PO
  # psi != phi: YP
  ratio <- as.numeric(exp(Z%*%(phi-psi)))
  kappa <- exp(Z%*%phi)
  w <- as.numeric((u^(-1/kappa) - 1)*ratio)
  r <- log(1+w)
  time <- switch (baseline,
                  "weibull" = (r/gamma[2])^(1/gamma[1])  # Weibull (JAGS parametrization)
  )
  return(time)
}



#' Random generation of survival data
#' @export
#' @description Function to generate a random sample of survival data.
#' @aliases rsurv
#' @param formula formula specifying the linear predictors
#' @param covariates data frame containing the covariates used to generate the survival times
#' @param baseline baseline model (currently only the Weibull distribution is available)
#' @param gamma baseline parameters
#' @param psi short-term regression coefficients
#' @param phi long-term regression coefficients
#' @param max_fu maximum follow-up time
#'
rsurv <- function(formula, covariates,
                      baseline = "weibull",
                      gamma, psi = NULL, phi = NULL, nu = NULL, max_fu){
  formula <- Formula::Formula(formula)
  mf <- stats::model.frame(formula=formula, data=covariates)
  Z <- stats::model.matrix(formula, data = mf, rhs = 1)[,-1, drop = FALSE]
  n <- nrow(Z)
  q <- ncol(Z)
  u <- stats::runif(n)
  t <- c()
  time <- c()
  status <- c()
  c <- runif(n, 0, max_fu)
  for(i in 1:n){
    t[i] <- rYP(u = u[i], Z = Z[i,], psi = psi, phi = phi, gamma = gamma, baseline = baseline[[1]])
    time[i] <- min(t[i], c[i])
    status[i] <- as.numeric(t[i] < c[i])
  }
  data <- data.frame(time = time, status = status)
  data <- dplyr::bind_cols(data, mf)
  return(data)
}
