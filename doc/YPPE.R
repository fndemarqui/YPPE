## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(YPPE)

mle <- yppe(Surv(time, status)~trt, data=gastric, n_int = 10, approach = "mle", init = 0)
summary(mle)
coef(mle)
vcov(mle)

