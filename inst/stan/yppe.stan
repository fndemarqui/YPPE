

#include /chunks/loglikspe.stan

data{
  int<lower=1> n;
  int<lower=1> m;
  int<lower=0> q;
  int<lower=0> p;
  int survreg;
  vector[n] status;
  int idt[n];
  matrix[q == 0 ? 0 : n, q] Z;
  matrix[p == 0 ? 0 : n, p] X;
  real<lower=0> tau;
  matrix[n, m] ttt;
  real h1_gamma;
  real h2_gamma;
  real mu_beta;
  real mu_psi;
  real mu_phi;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_psi;
  real<lower=0> sigma_phi;
  int<lower=0, upper=1> approach;
}


parameters{
  vector[q] psi;
  vector[q] phi;
  vector[p == 0 ? 0 : p] beta;
  vector<lower=0>[m] gamma;
}


// transformed parameters{
//   vector[n] loglik;
//   loglik = loglik1_pe(status, Z, tau, ttt, idt, gamma, psi, phi);
// }


model{
  vector[n] loglik;
  vector[q == 0 ? 0 : n] lp_short;
  vector[q == 0 ? 0 : n] lp_long;
  vector[p == 0 ? 0 : n] lp_const;
  vector[n] ratio;
  vector[n] lht0 = log(gamma[idt]); // - log(tau);
  vector[n] Ht0 = ttt*gamma;
  if(survreg == 1){
        lp_short = Z*psi;
        lp_long = Z*phi;
        ratio = exp(Z*(psi-phi));
        loglik = loglik_yp1(status, lht0, Ht0, lp_short, lp_long, ratio, n, q);
  }else if(survreg == 2){
    lp_short = Z*psi;
    lp_long = Z*phi;
    lp_const = X*beta;
    ratio = exp(Z*(psi-phi));
    loglik = loglik_yp2(status, lht0, Ht0, lp_short, lp_long, lp_const, ratio, n, q);
  }else if(survreg == 3){
    lp_const = X*beta;
    loglik = loglik_ph(status, lht0, Ht0, lp_const, n, p);
  }else{
    lp_const = X*beta;
    loglik = loglik_po(status, lht0, Ht0, lp_const, n, p);
  }

  target += sum(loglik);
  if(approach==1){
    gamma ~ lognormal(h1_gamma, h2_gamma);
    psi ~ normal(mu_psi, sigma_psi);
    phi ~ normal(mu_phi, sigma_phi);
    beta ~ normal(mu_beta, sigma_beta);
  }
}
