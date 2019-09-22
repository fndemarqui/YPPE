

#include /chunks/loglikspe.stan

data{
  int<lower=1> n;
  int<lower=1> m;
  int<lower=1> q;
  int<lower=1> p;
  vector[n] status;
  int idt[n];
  matrix[n,q] Z;
  matrix[n,p] X;
  real<lower=0> tau;
  matrix[n,m] ttt;
  real mu_lambda;
  real mu_psi;
  real mu_phi;
  real mu_beta;
  real<lower=0> sigma_lambda;
  real<lower=0> sigma_psi;
  real<lower=0> sigma_phi;
  real<lower=0> sigma_beta;
  int<lower=0, upper=1> approach;
}


parameters{
  vector[q] psi;
  vector[q] phi;
  vector[p] beta;
  vector[m] log_lambda;
}


transformed parameters{
  vector[m] lambda;
  vector[n] loglik;
  for(k in 1:m){
    lambda[k] = exp(log_lambda[k]);
  }
  loglik = loglik2_pe(status, Z, X, tau, ttt, idt, lambda, psi, phi, beta);
}


model{
  target += loglik;
  if(approach==1){
    log_lambda ~ normal(mu_lambda, sigma_lambda);
    psi ~ normal(mu_psi, sigma_psi);
    phi ~ normal(mu_phi, sigma_phi);
    beta ~ normal(mu_beta, sigma_beta);
  }
}