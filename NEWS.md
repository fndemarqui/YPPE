
# YPPE Version 1.1.0

- survit function returns now a list containing two objects: time and surv associated with the observed survival times and their corresponding survival probabilities.

- inclusion of model.matrix method

- implementation of phpe() and phpo() functions to fit proportional hazards (PH) and proportional odds (PO) models.
 
- implementation of anova.yppe() function to carry out likelihood ratio tests.
 
- implemantation of AIC.yppe() function to compute AIC (Akaike information criterion).
 
- implemantation of logLik.yppe() function to extract the log-likelihood value of a fitted model.
  
- the YPPE package now requires rstan Version 2.26.  
  
- update package reference.
  

# YPPE version 1.0.1

- coef() function to extract the regression coefficientes under the MLE approach
- vcov() function to compute the variance-covariance matrix associated with the regression coefficients under the MLE approach
- confint() function to compute the 100(1-alpha)% confidence intervals under the MLE approahc

- bug fix: yppe() function now works properly if both the number of intervals and the time grid (n_int and rho) are passed as arguments; the default time grid is now computed by using all distinct observed failure times as the endpoints of the intervals



