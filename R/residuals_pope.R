


residuals.pope <- function(object, type = c("cox-snell", "martingale", "deviance")){
  type <- match.arg(type)
  mf <- object$mf
  labels <- names(mf)[-1]
  time <- sort( stats::model.response(mf)[,1])
  status <- sort( stats::model.response(mf)[,2])
  data <- data.frame(cbind(time, status, mf[,-1]))
  names(data) <- c("time", "status", names(mf)[-1])
  rho <- object$rho
  n_int <- object$n_int
  tau <- object$tau
  formula <- object$formula
  Z <- as.matrix(stats::model.matrix(formula, data = mf, rhs = 1)[,-1])
  St <- c()
  n <- length(time)

  if(object$approach=="mle"){
    par <- object$fit$par
    for(i in 1:n){
      St[i] <- popeSurv(time[i], Z[i,], par, rho, tau, n_int)
    }
  }else{ # Bayesian approach
    samp <- rstan::extract(object$fit)
    for(i in 1:n){
      aux <- apply(par, 1, popeSurv, time=time[i], z=Z[i,], rho=rho, tau=tau, n_int=n_int)
      St[i] <- apply(aux, 1, mean)
    }
  }

  r <- -log(St)
  if(type == "cox-snell"){
    return(r)
  }else if(type == "martingale"){
    m <- status - r
    return(m)
  }else{
    m <- status - r
    d <- sign(m)*( (-2*( m + status*log(status - m) ))^0.5  )
    return(d)
  }
}
