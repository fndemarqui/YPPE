
library(rstan)


model_code <-
  '
functions{

real dbeta(real x, real shape1, real shape2){
  return(exp(beta_lpdf(x|shape1, shape2)));
}

real pbeta(real x, real shape1, real shape2){
  return(exp(beta_lcdf(x|shape1, shape2)));
}



//-------------------------------------------------------------
// Likelihood function for model M1:
vector loglik1_bp(vector status, matrix Z, matrix g, matrix G, real tau,
                  vector gamma, vector psi, vector phi){

  int q = cols(Z);
  int n = rows(Z);
  vector[n] lht0;
  vector[n] Ht0;
  vector[n] Ft0;
  vector[n] lp_short;
  vector[n] lp_long;
  vector[n] log_ht;
  vector[n] log_St;
  vector[n] loglik;
  vector[n] theta;
  vector[n] logmix;
  vector[n] aux;

  lht0 = log(g*gamma) - log(tau);
  Ht0 = G*gamma;
  Ft0 = -expm1(-Ht0);
  lp_short = Z*psi;
  lp_long = Z*phi;
  theta = exp(lp_long);
  aux = lp_long - Ht0;

  for(i in 1:n){
    logmix[i] = log_mix(Ft0[i], lp_short[i], lp_long[i]);
    log_ht[i] = lp_short[i] + lp_long[i] - logmix[i] + lht0[i];
    log_St[i] = -theta[i]*(logmix[i] - aux[i]);
    loglik[i] = status[i]*log_ht[i] + log_St[i];
  }
  return loglik;
}


//-------------------------------------------------------------
// Likelihood function for model M2:
vector loglik2_bp(vector status, matrix Z, matrix g, matrix G, real tau,
                  vector gamma, vector psi, vector phi){

  int q = cols(Z);
  int n = rows(Z);
  vector[n] Rt0;
  vector[n] lrt0;
  vector[n] lp_short;
  vector[n] theta_L;
  vector[n] log_ht;
  vector[n] log_St;
  vector[n] ratio;
  vector[n] aux;
  vector[n] loglik;

  lrt0 = log(g*gamma) - log(tau);
  Rt0 = G*gamma;
  lp_short = Z*psi;
  theta_L = exp( Z*phi);
  ratio = exp( Z*(psi-phi));

  for(i in 1:n){
    aux[i] = ratio[i]*Rt0[i];
    log_ht[i] = lp_short[i] - log1p(aux[i] ) + lrt0[i];
    log_St[i] = - theta_L[i]*log1p(aux[i]);
    loglik[i] = status[i]*log_ht[i] + log_St[i];
  }

  return loglik;
}


//-------------------------------------------------------------
// Likelihood function for model M3:
vector loglik3_bp(vector status, matrix Z, matrix X, matrix g, matrix G,  real tau,
                  vector gamma, vector psi, vector phi, vector beta){

  int q = cols(Z);
  int p = cols(X);
  int n = rows(Z);

  vector[n] lht0;
  vector[n] Ht0;
  vector[n] Ft0;
  vector[n] lp_short;
  vector[n] lp_long;
  vector[n] lp_const;
  vector[n] log_ht;
  vector[n] log_St;
  vector[n] ratio;
  vector[n] loglik;
  vector[n] theta;
  vector[n] logmix;
  vector[n] aux;

  lht0 = log(g*gamma) - log(tau);
  Ht0 = G*gamma;
  Ft0 = -expm1(-Ht0);
  lp_short = Z*psi;
  lp_long = Z*phi;
  lp_const = X*beta;
  theta = exp(lp_long);
  aux = lp_long - Ht0;

  for(i in 1:n){
    logmix[i] = log_mix(Ft0[i], lp_short[i], lp_long[i]);
    log_ht[i] = lp_short[i] + lp_long[i] + lp_const[i] - logmix[i] + lht0[i];
    log_St[i] = -theta[i]*exp(lp_const[i])*(logmix[i] - aux[i]);
    loglik[i] = status[i]*log_ht[i] + log_St[i];
  }

  return loglik;
}


//-------------------------------------------------------------
// Likelihood function for model M4:
vector loglik4_bp(vector status, matrix Z, matrix X, matrix g, matrix G,  real tau,
                  vector gamma, vector psi, vector phi, vector beta){


  int q = cols(Z);
  int p = cols(X);
  int n = rows(Z);

  vector[n] Rt0;
  vector[n] lrt0;
  vector[n] lp_short;
  vector[n] lp_const;
  vector[n] theta_L;
  vector[n] log_ht;
  vector[n] log_St;
  vector[n] ratio;
  vector[n] aux;
   vector[n] loglik;

  lrt0 = log(g*gamma) - log(tau);
  Rt0 = G*gamma;
  lp_short = Z*psi;
  lp_const = X*beta;
  theta_L = exp( Z*phi);
  ratio = exp( Z*(psi-phi));

  for(i in 1:n){
    aux[i] = ratio[i]*Rt0[i];
    log_ht[i] = lp_short[i] + lp_const[i] - log1p(aux[i]) + lrt0[i];
    log_St[i] = - theta_L[i]*exp(lp_const[i])*log1p(aux[i]);
    loglik[i] = status[i]*log_ht[i] + log_St[i];
  }

  return loglik;
}

}

'

expose_stan_functions(stanc(model_code = model_code))


formula <- Surv(time, status)~trt
data <- gastric
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
tau <- max(time)
degree <- ceiling(sqrt(length(time)))

M <- 1
approahc <- 0
gamma <- rexp(degree)
psi <- rnorm(q)
phi <- rnorm(q)

gamma <- rep(1,degree)
psi <- rep(0,q)
phi <- rep(0,q)


aux <- bp(time, degree, tau)
g=aux$b
G=aux$B
m <- degree

loglik1 <- loglik1_bp(status, Z, g=g, G=G, tau,
           gamma, psi,  phi)
loglik2 <- loglik2_bp(status, Z, g=g, G=G, tau,
                     gamma, psi,  phi)

sum(loglik1)
sum(loglik2)

