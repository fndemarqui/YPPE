functions{

  //-------------------------------------------------------------
  // Likelihood function for model YP1 model:

  vector loglik_yp1(vector status, vector lht0, vector Ht0, vector lp_short, vector lp_long, vector ratio, int n, int q){
    vector[n] Rt0;
    vector[n] log_ht;
    vector[n] log_St;
    vector[n] loglik;
    vector[n] theta;
    vector[n] aux;

    Rt0 = expm1(Ht0);
    theta = exp(lp_long);

    aux = ratio .* Rt0;
    log_ht = lp_short - log1p(aux) + lht0 + Ht0;
    log_St = -theta .* log1p(aux);
    loglik = status .* log_ht + log_St;

    return loglik;
  }


  //-------------------------------------------------------------
  // Likelihood function for model YP2 model:

  vector loglik_yp2(vector status, vector lht0, vector Ht0, vector lp_short, vector lp_long, vector lp_const, vector ratio, int n, int q){
    vector[n] Rt0;
    vector[n] log_ht;
    vector[n] log_St;
    vector[n] loglik;
    vector[n] theta;
    vector[n] aux;

    Rt0 = expm1(Ht0);
    theta = exp(lp_long + lp_const);

    aux = ratio .* Rt0;
    log_ht = lp_short + lp_const - log1p(aux) + lht0 + Ht0;
    log_St = -theta .* log1p(aux);
    loglik = status .* log_ht + log_St;

    return loglik;
  }

  //-------------------------------------------------------------
  // Likelihood function for model PH model:
  vector loglik_ph(vector status, vector lht0, vector Ht0, vector lp, int n, int p){
    vector[n] loglik = status .* (lht0 + lp) - Ht0 .* exp(lp);
    return loglik;
  }

  //-------------------------------------------------------------
  // Likelihood function for model PO model:
  vector loglik_po(vector status, vector lht0, vector Ht0, vector lp, int n, int p){
    vector[n] aux = exp(lp) .* expm1(Ht0);
    vector[n] loglik = status .* (lht0 + lp + Ht0) - (1+status) .* log1p(aux);
    return loglik;
  }

}
