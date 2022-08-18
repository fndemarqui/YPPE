## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(YPPE)
data(gastric)

ph <- phpe(Surv(time, status)~trt, data = gastric, n_int = 10, init = 0)
po <- pope(Surv(time, status)~trt, data = gastric, n_int = 10, init = 0)
yp <- yppe(Surv(time, status)~trt, data = gastric, n_int = 10, init = 0)
anova(ph, po, yp)
AIC(ph, po, yp)

