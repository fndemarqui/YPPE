

extract_survregs <- function(models){
  extract_survreg <- function(fit){
    return(fit$survreg)
  }
  return(unlist(lapply(models, extract_survreg)))
}

extract_formulas <- function(object){
  return(object$formula)
}


#---------------------------------------------
#' anova method for yppe models
#'
#' @aliases anova.yppe
#' @description Compute analysis of variance (or deviance) tables for one or more fitted model objects.
#' @importFrom stats anova
#' @export
#' @param ... further arguments passed to or from other methods.
#' @return  the ANOVA table.
#'
anova.yppe <- function(...){
  models <- c(as.list(environment()), list(...))

  J <- nargs()
  labels <- paste0("Model ", 1:J, ":")
  k <- c()
  df <- c()
  k[J] <- models[[J]]$npar
  LR <- c()
  p.value <- c()
  loglik <- c()
  loglik[J] <- models[[J]]$fit$value
  for(j in 1:(J-1)){
    loglik[j] <- models[[j]]$fit$value
    LR[j] <- 2*(models[[J]]$fit$value - models[[j]]$fit$value)
    k[j] <- length(models[[j]]$fit$par)
    df[j] <- k[J]-k[j]
    p.value[j] <- stats::pchisq(LR[j], df = df[j], lower.tail = FALSE)
  }


  tab <- cbind("loglik" = loglik[-J], LR, df, 'Pr(>Chi)' = p.value)
  aux <- matrix(c(loglik[J], NA, NA, NA), nrow = 1)
  tab <- rbind(tab, aux)
  #tab <- cbind(survreg, tab)
  rownames(tab) <- labels
  formulas <- sapply(models, extract_formulas)

  survreg <- extract_survregs(models)

  cat("\n")
  for(j in 1:J){
    cat("Model", j, "(", survreg[j], "): ", deparse(formulas[[j]]), "\n")
  }
  cat("--- \n")
  stats::printCoefmat(tab, P.values=TRUE, has.Pvalue = TRUE, na.print = "-")
}
