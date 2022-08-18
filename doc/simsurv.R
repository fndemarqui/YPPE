## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(YPPE)

## ----covariates---------------------------------------------------------------

# fixing the seed for the rng:
set.seed(1234567890)

gamma <- c(1.5, 1);  # shape and scale parameters
n <- 100

# generating the set of explanatory variables:
covariates <- data.frame(
  trt = rbinom(n, 1, 0.5),
  age = rnorm(n)
)


## -----------------------------------------------------------------------------
# generating the data set:
simdata <- rsurv(
  ~ trt + age, covariates = covariates,
  baseline = "weibull",
  gamma = gamma,
  psi = c(2, -1),
  phi = c(2, -1),
  max_fu = 5
)

ph <- phpe(Surv(time, status)~trt+age, data = simdata, n_int = 10, init = 0)
po <- pope(Surv(time, status)~trt+age, data = simdata, n_int = 10, init = 0)
yp <- yppe(Surv(time, status)~trt+age, data = simdata, n_int = 10, init = 0)
anova(ph, po, yp)


## ---- po----------------------------------------------------------------------
# generating the data set:
simdata <- rsurv(
  ~ trt + age, covariates = covariates,
  baseline = "weibull",
  gamma = gamma,
  psi = c(2, -1),
  phi = c(0, 0),
  max_fu = 5
)

ph <- phpe(Surv(time, status)~trt+age, data = simdata, n_int = 10, init = 0)
po <- pope(Surv(time, status)~trt+age, data = simdata, n_int = 10, init = 0)
yp <- yppe(Surv(time, status)~trt+age, data = simdata, n_int = 10, init = 0)
anova(ph, po, yp)


## ---- yp----------------------------------------------------------------------
# generating the data set:
simdata <- rsurv(
  ~ trt + age, covariates = covariates,
  baseline = "weibull",
  gamma = gamma,
  psi = c(2, -1),
  phi = c(-1, -1),
  max_fu = 5
)

ph <- phpe(Surv(time, status)~trt+age, data = simdata, n_int = 10, init = 0)
po <- pope(Surv(time, status)~trt+age, data = simdata, n_int = 10, init = 0)
yp <- yppe(Surv(time, status)~trt+age, data = simdata, n_int = 10, init = 0)
anova(ph, po, yp)


