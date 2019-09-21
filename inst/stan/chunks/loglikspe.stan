functions{

//-------------------------------------------------------------
// Likelihood function for model M1:
vector loglik1_pe(vector status, matrix Z, real tau, matrix ttt,  int[] idt, vector lambda, vector psi, vector phi){

  int q = cols(Z);
  int n = rows(Z);
  vector[n] lht0;
  vector[n] Ft0;
  vector[n] Ht0;
  vector[n] lp_short;
  vector[n] lp_long;
  vector[n] log_ht;
  vector[n] log_St;
  vector[n] loglik;
  vector[n] theta;
  vector[n] logmix;
  vector[n] aux;


  lht0 = log(lambda[idt]) - log(tau);
  Ht0 = ttt*lambda;
  Ft0 = -expm1(-Ht0);
  lp_short = Z*psi;
  lp_long = Z*phi;
  theta = exp(lp_long);
  aux = lp_long - Ht0;

  for(i in 1:n)
  {
    logmix[i] = log_mix(Ft0[i], lp_short[i], lp_long[i]);
    log_ht[i] = lp_short[i] + lp_long[i] - logmix[i] + lht0[i];
    log_St[i] = -theta[i]*(logmix[i] - aux[i]);
    loglik[i] = status[i]*log_ht[i] + log_St[i];
  }
  return loglik;
}



//-------------------------------------------------------------
// Likelihood function for model M3:
vector loglik2_pe(vector status, matrix Z, matrix X, real tau,
                  matrix ttt, int[] idt, vector lambda, vector psi, vector phi, vector beta){

  int q = cols(Z);
  int p = cols(X);
  int n = rows(Z);

  vector[n] lht0;
  vector[n] Ft0;
  vector[n] Ht0;
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

  lht0 = log(lambda[idt]) - log(tau);
  Ht0 = ttt*lambda;
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



}
