# August 17 2022

---

## version 1.1.0

### New features

 - implementation of phpe() and phpo() functions to fit proportional hazards (PH) and proportional odds (PO) models.
 
 - implementation of anova.yppe() function to carry out likelihood ratio tests.
 
 - implemantation of AIC.yppe() function to compute AIC (Akaike information criterion).
 
  - implemantation of logLik.yppe() function to extract the log-likelihood value of a fitted model.
  
  - implementation of rsurv() function to generate survival data under the YP (and PH and PO as well) model (currently only with Weibull baseline).


# June 17 2020

---

## version 1.0.2

### New features

- survit function returns now a list containing two objects: time and surv associated with the observed survival times and their corresponding survival probabilities.

- inclusion of model.matrix method



# January 09 2020

---

## version 1.0.1

### New features

- coef() function to extract the regression coefficientes under the MLE approach
- vcov() function to compute the variance-covariance matrix associated with the regression coefficients under the MLE approach
- confint() function to compute the 100(1-alpha)% confidence intervals under the MLE approahc


### Bug fixes

- yppe() function now works properly if both the number of intervals and the time grid (n_int and rho) are passed as arguments; the default time grid is now computed by using all distinct observed failure times as the endpoints of the intervals


---

