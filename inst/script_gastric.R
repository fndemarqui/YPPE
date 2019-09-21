
library("crossSurv")


yppe_mle <- yppe(Surv(time, status)~trt, data=gastric, approach="mle")
yppe_bayes <- yppe(Surv(time, status)~trt, data=gastric, approach="bayes")
ypbp_mle1 <- ypbp(Surv(time, status)~trt, data=gastric, approach="mle", baseline="hazard")
ypbp_mle2 <- ypbp(Surv(time, status)~trt, data=gastric, approach="mle", baseline="odds")
ypbp_bayes1 <- ypbp(Surv(time, status)~trt, data=gastric, approach="bayes", baseline="hazard")
ypbp_bayes2 <- ypbp(Surv(time, status)~trt, data=gastric, approach="bayes", baseline="odds")


summary(mle)
summary(bayes)

newdata1 <- data.frame(trt=0)
newdata2 <- data.frame(trt=1)
t_mle <- crossTime(mle, newdata1, newdata2)
t_bayes <- crossTime(bayes, newdata1, newdata2)
t_mle
t_bayes

newdata <- data.frame(trt=as.factor(0:1))
St_mle <- survfit(mle, newdata)
St_bayes <- survfit(bayes, newdata)


ekm <- survfit(Surv(time, status)~trt, data=gastric)
plot(ekm, col=1:2)
time <- gastric$time
lines(time, St_mle[[1]])
lines(time, St_mle[[2]], col=2)
lines(time, St_bayes[[1]], lty=2)
lines(time, St_bayes[[2]], col=2, lty=2)
abline(v=t_mle, col="blue")
abline(v=t_bayes, col="green")








