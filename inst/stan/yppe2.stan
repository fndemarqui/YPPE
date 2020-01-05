

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
  vector[q] z_mean;
  vector[q] z_sd;
  vector[p] x_mean;
  vector[p] x_sd;
  real h1_gamma;
  real h2_gamma;
  real mu_psi;
  real mu_phi;
  real mu_beta;
  real<lower=0> sigma_psi;
  real<lower=0> sigma_phi;
  real<lower=0> sigma_beta;
  int<lower=0, upper=1> approach;
}

transformed data{
  matrix[n, q] Z_std;
  matrix[n, p] X_std;
  for(j in 1:q){
    for(i in 1:n){
      Z_std[i, j] = (Z[i, j] - z_mean[j])/z_sd[j];
    }
  }
  for(j in 1:p){
    for(i in 1:n){
      X_std[i, j] = (X[i, j] - x_mean[j])/x_sd[j];
    }
  }
}

parameters{
  vector[q] psi_std;
  vector[q] phi_std;
  vector[p] beta_std;
  vector<lower=0>[m] gamma_std;
}


model{
  vector[n] loglik = loglik2_pe(status, Z_std, X_std, ttt, idt, gamma_std, psi_std, phi_std, beta_std);
  target += sum(loglik);
  if(approach==1){
    gamma_std ~ lognormal(h1_gamma, h2_gamma);
    psi_std ~ normal(mu_psi, sigma_psi);
    phi_std ~ normal(mu_phi, sigma_phi);
    beta_std ~ normal(mu_beta, sigma_beta);
  }
}


generated quantities{
  vector[q] psi;
  vector[q] phi;
  vector[p] beta;
  vector<lower=0>[m] gamma;
  psi = psi_std ./ z_sd;
  phi = phi_std ./ z_sd;
  beta = beta_std ./ x_sd;
  gamma = gamma_std*exp( -sum( ((psi + phi) .* z_mean) ) - sum( beta .* x_mean) - log(tau) );
  if(approach==1){
    vector[n] loglik;
    loglik = loglik2_pe(status, Z_std, X_std, ttt, idt, gamma_std, psi_std, phi_std, beta_std);
  }
}


