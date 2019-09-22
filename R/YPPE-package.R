#' The 'YPPE' package.
#'
#' @description Semiparametric modeling of lifetime data with crossing survival curves via Yang and Prentice model with piecewise exponential baseline distribution curves. Model fitting can be carried out by either likelihood-based or Bayesian approaches. The package interfaces with Stan to fit the YPPE model. It also provides a function to obtain point and interval estimates for the crossing survival times.
#'
#' @docType package
#' @name YPPE-package
#' @aliases YPPE
#' @useDynLib YPPE, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import survival
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (2019). RStan: the R interface to Stan. R package version 2.19.2. https://mc-stan.org
#'
#' Yang, S. and Prentice, R. L. (2005). Semiparametric analysis of short-term and long-term hazard ratios with two-sample survival data. Biometrika 92, 1-17.
#'
NULL
