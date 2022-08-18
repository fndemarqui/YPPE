#' The 'YPPE' package.
#'
#' @description Semiparametric modeling of lifetime data with crossing survival curves via Yang and Prentice model with piecewise exponential baseline distribution. Details about the model can be found in \insertCite{2021_Demarqui}{YPPE} <doi.org/10.1214/20-BJPS471>. Model fitting carried out via likelihood-based and Bayesian approaches. The package also provides point and interval estimation for the crossing survival times.
#'
#' @docType package
#' @name YPPE-package
#' @aliases YPPE
#' @useDynLib YPPE, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom stats runif vcov
#' @import survival
#' @importFrom rstan sampling
#'
#' @references
#'
#' \insertRef{2021_Demarqui}{YPPE}
#'
#' \insertRef{2005_YangPrentice}{YPPE}
#'
#' Stan Development Team (2019). RStan: the R interface to Stan. R package version 2.19.2. https://mc-stan.org
#'
#'
NULL
