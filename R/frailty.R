

frailty <- function(x, distribution = c("gamma", "gaussian"), ...){
  distribution <- match.arg(distribution)
  id_frailty <- as.numeric(as.factor(x))
  L <- length(unique(id_frailty))
  attributes(id_frailty) <- list(L = L, distribution = distribution)
  return(id_frailty)
}


# #-----------------------------------------------------
# # testes
#
# library(survival)
# library(tidyverse)
#
#
# data(cancer, package="survival")
#
# glimpse(kidney)
#
# kidney <- kidney %>%
#   mutate(
#     sex = as.factor(sex)
#   )
#
# formula <- Surv(time, status) ~ sex  + frailty(id, distribution="gaussian")
# formula
# mf <- model.frame(formula, data = kidney)
# glimpse(mf)
# attributes(mf)
# aux <- grep("frailty", names(mf))
# if(length(aux)>0){
#   id <- mf[, aux]
#   mf <- mf[, -aux, drop = FALSE]
#   L <- attr(id, "L")
#   frailty_dist <-  attr(id, "distribution")
# }else{
#   L <- 0
#   frailty_dist <-  0
# }
#
#
# # x <- as.factor(rep(1:4, each = 2))
# # teste <- frailty(x)
# # teste
# # attr(teste, "L")
# # attr(teste, "distribution")
