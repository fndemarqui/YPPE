
library(YPBP)
library(YPPE)

mle1 <- yppe(Surv(time, status)~x1+x2, data=sim_data, init=0, n_int=1)
mle2 <- yppe(Surv(time, status)~x1|x2, data=sim_data, init=0, n_int=1)
summary(mle1)
summary(mle2)

sqrt(diag(vcov(mle1)))
sqrt(diag(vcov(mle2)))



vcov(mle1)
vcov(mle2)
coef(mle1)
coef(mle2)
confint(mle1)
confint(mle2)


bayes1 <- yppe(Surv(time, status)~x1+x2, data=sim_data, approach="bayes")
bayes2 <- yppe(Surv(time, status)~x1|x2, data=sim_data, approach="bayes")

summary(bayes1)
summary(bayes2)


data <- ipass
formula <- Surv(time, status)~arm
baseline <- "hazard"

library(YPPE)
ipass$arm <- ipass$arm
pe <- yppe(Surv(time, status)~arm, init=0, data=ipass)
summary(pe)

library(YPPE)


summary(pe2)

pe$fit$par[1:2]/sd(ipass$arm)

odds$fit$par



hazard <- yppe(Surv(time, status)~arm, baseline="hazard", init=0, data=ipass)
odds <- yppe(Surv(time, status)~arm, baseline="odds", init=0, data=ipass)
sqrt(diag(vcov(hazard)))
sqrt(diag(vcov(odds)))



summary(hazard)
summary(odds)



names(ovarian)
fit1 <- yppe(Surv(time, status)~ESA+cytored|EpoR+EphB4+age,
             baseline="hazard", data=ovarian)
fit2 <- yppe(Surv(time, status)~ESA+cytored|EpoR+EphB4+age,
             baseline="odds", data=ovarian)
summary(fit1)
summary(fit2)

labels <- fit1$labels
mf <- colnames(fit1$mf[,-1])

grepl(labels, mf)
gsub(labels, mf)

match.arg(labels, mf, several.ok = TRUE)


fit1 <- yppe(Surv(time, status)~ESA|age,
             baseline="hazard", data=ovarian)



library(YPBP)
ovarian$age <- scale(ovarian$age)
mle <- yppe(Surv(time, status)~ESA+cytored|EpoR+EphB4+age,
             baseline="hazard", data=ovarian, init=0)
summary(mle)



summary(yppe(Surv(time, status)~ESA+cytored|EpoR+EphB4+age,
            baseline="hazard", data=ovarian, init=0))

summary(yppe(Surv(time, status)~ESA+cytored|EpoR+EphB4+age,
             baseline="odds", data=ovarian, init=0))


mle <- yppe(Surv(time, status)~ESA+cytored+EpoR+EphB4+age,
            baseline="hazard", data=ovarian)
summary(mle)



bayes <- yppe(Surv(time, status)~ESA+cytored|EpoR+EphB4+age,
            baseline="hazard", approach="bayes",data=ovarian)
summary(bayes)


bayes <- yppe(Surv(time, status)~ESA+cytored|EpoR+EphB4+age,
              baseline="odds", approach="bayes",data=ovarian)
summary(bayes)


bayes <- yppe(Surv(time, status)~ESA+cytored+EpoR+EphB4+age,
              baseline="odds", approach="bayes",data=ovarian)
summary(bayes)






formula <- Surv(time, status)~ESA+cytored|EpoR+EphB4+age
data <- ovarian


data(ovarian, package="YPBP")
mle <- yppe(Surv(time, status)~ESA+cytored+EpoR+EphB4+age,
            data=ovarian, init=0, standardization = "typeI", n_int=30)

summary(mle)


bayes <- yppe(Surv(time, status)~ESA+cytored+EpoR+EphB4+age,
            data=ovarian, init=0, standardization = "typeI", n_int=14,
            approach="bayes")

summary(bayes)


library(YPPE)
fit <- yppe(Surv(time, status)~arm, data=ipass, init=0)
summary(fit)



fit1 <- yppe(Surv(time, status)~arm, data=ipass, init=0,
             standardization="typeI")
summary(fit1)

fit2 <- yppe(Surv(time, status)~arm, data=ipass, init=0,
             standardization="typeII")
summary(fit2)

cbind(fit1$fit$par, fit2$fit$par, fit1$fit$par/fit2$fit$par)



mle <- yppe(Surv(time, status)~trt, data=gastric, init=0, approach = "mle")
bayes <- yppe(Surv(time, status)~trt, data=gastric, approach = "bayes")
summary(mle)
summary(bayes)



mle <- yppe(Surv(time, status)~arm, data=ipass, init=0, approach = "mle")
bayes <- yppe(Surv(time, status)~arm, data=ipass, approach = "bayes")
summary(mle)
summary(bayes)
