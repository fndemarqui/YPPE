#' The 'YPPE' package.
#'
#' @description Semiparametric modeling of lifetime data with crossing survival curves via Yang and Prentice model with piecewise exponential baseline distribution curves. Details about the model can be found in Demarqui and Mayrink (2019) <arXiv:1910.02406>. Model fitting carried out via likelihood-based and Bayesian approaches. The package also provides point and interval estimation for the crossing survival times.
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
#'
#' Demarqui, F. N. and Mayrink, V. D. (2019). A fully likelihood-based approach to model survival data with crossing survival curves. <arXiv:1910.02406>
#'
#' Stan Development Team (2019). RStan: the R interface to Stan. R package version 2.19.2. https://mc-stan.org
#'
#' Yang, S. and Prentice, R. L. (2005). Semiparametric analysis of short-term and long-term hazard ratios with two-sample survival data. Biometrika 92, 1-17.
#'
NULL
