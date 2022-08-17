#' Print the summary.phpe output
#'
#' @export
#' @param x an object of the class summary.phpe.
#' @param ... further arguments passed to or from other methods.
#' @return a summary of the fitted model.
print.summary.phpe <- function(x, ...){
  if(x$approach=="mle"){
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Proportional hazards coefficients:\n")
    stats::printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
    cat("\n")
    cat("--- \n")
    cat("loglik =", x$loglik, " ", "AIC =", x$AIC,"\n")
  }else{
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Proportional hazards coefficients:\n")
    print(x$coefficients)
    cat("--- \n")
    cat("Inference for Stan model: ", x$model_name, '.\n', sep = '')
    cat(x$chains, " chains, each with iter=", x$iter,
        "; warmup=", x$warmup, "; thin=", x$thin, "; \n",
        "post-warmup draws per chain=", x$n_kept[1], ", ",
        "total post-warmup draws=", sum(x$n_kept), ".\n\n", sep = '')

  }

}


#---------------------------------------------

#' Summary for the yppe model
#'
#' @aliases summary.phpe
#' @export
#' @param object an objecto of the class 'yppe'.
#' @param ... further arguments passed to or from other methods.
#'
summary.phpe <- function(object, ...){
  p <- object$p
  if(object$approach=="mle"){
    n <- object$n
    k <- p+object$n_int

    loglik <- object$fit$value
    AIC <- -2*loglik + 2*k
    BIC <- -2*loglik + k*log(n)

    labels <- object$labels
    coefficients <- object$fit$par[1:p]
    vcov <- MASS::ginv(-object$fit$hessian)[1:p, 1:p, drop = FALSE]

    se <- sqrt(diag(vcov))
    zval <- coefficients / se
    TAB <- cbind(Estimate = coefficients,
                 StdErr = se,
                 z.value = zval,
                 p.value = 2*stats::pnorm(-abs(zval)))

    rownames(TAB) <- labels

    res <- list(call=object$call,
                coefficients=TAB,
                loglik=loglik, AIC=AIC,
                approach=object$approach)


   # Bayesiam output:
  }else{
    labels <- object$labels
    s <- rstan::summary(object$fit, pars=c("beta"))
    TAB <- round(s$summary, digits = 3)
    rownames(TAB) <- labels
    n_kept <- object$fit@sim$n_save - object$fit@sim$warmup2
    res <- list(call=object$call,
                coefficients1=TAB,
                n_kept=n_kept, model_name=object$fit@model_name,
                chains=object$fit@sim$chains, warmup=object$fit@sim$warmup,
                thin=object$fit@sim$thin, iter=object$fit@sim$iter, approach=object$approach)

  }

  class(res) <- "summary.phpe"

  return(res)
}


