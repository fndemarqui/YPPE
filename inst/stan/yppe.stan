
#include /chunks/loglikspe.stan

data{
  int<lower=1> n;
  int<lower=1> m;
  int<lower=1> q;
  vector[n] status;
  int idt[n];
  matrix[n,q] Z;
  real<lower=0> tau;
  matrix[n,m] ttt;
  vector[q] z_mean;
  vector[q] z_sd;
  real h1_gamma;
  real h2_gamma;
  real mu_psi;
  real mu_phi;
  real<lower=0> sigma_psi;
  real<lower=0> sigma_phi;
  int<lower=0, upper=1> approach;
}

transformed data{
  matrix[n, q] Z_std;
  for(j in 1:q){
    for(i in 1:n){
      Z_std[i, j] = (Z[i, j] - z_mean[j])/z_sd[j];
    }
  }
}

parameters{
  vector[q] psi_std;
  vector[q] phi_std;
  vector<lower=0>[m] gamma_std;
}

model{
  vector[n] loglik = loglik1_pe(status, Z_std, ttt, idt, gamma_std, psi_std, phi_std);
  target += sum(loglik);
  if(approach==1){
    gamma_std ~ lognormal(h1_gamma, h2_gamma);
    psi_std ~ normal(mu_psi, sigma_psi);
    phi_std ~ normal(mu_phi, sigma_phi);
  }
}

generated quantities{
  vector[q] psi;
  vector[q] phi;
  vector<lower=0>[m] gamma;
  psi = psi_std ./ z_sd;
  phi = phi_std ./ z_sd;
  gamma = gamma_std*exp( -sum( ((psi + phi) .* z_mean) ) - log(tau) );
  if(approach==1){
    vector[n] loglik;
    loglik = loglik1_pe(status, Z_std, ttt, idt, gamma_std, psi_std, phi_std);
  }
}

