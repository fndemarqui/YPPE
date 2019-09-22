
library(YPPE)

mle <- yppe(Surv(time, status)~arm, data=ipass, approach="mle", n_int=10)
bayes <- yppe(Surv(time, status)~arm, data=ipass, approach="bayes", n_int=10)

summary(mle)
summary(bayes)

newdata1 <- data.frame(arm=0)
newdata2 <- data.frame(arm=1)
t_mle <- crossTime(mle, newdata1, newdata2)
t_bayes <- crossTime(bayes, newdata1, newdata2)
t_mle
t_bayes

newdata <- data.frame(arm=as.factor(0:1))
St_mle <- survfit(mle, newdata)
St_bayes <- survfit(bayes, newdata)


ekm <- survfit(Surv(time, status)~arm, data=ipass)
plot(ekm, col=1:2)
time <- sort(ipass$time)
lines(time, St_mle[[1]])
lines(time, St_mle[[2]], col=2)
lines(time, St_bayes[[1]], lty=2)
lines(time, St_bayes[[2]], col=2, lty=2)
abline(v=t_mle, col="blue")
abline(v=t_bayes, col="green")



t_mle <- crossTime(mle, newdata1, newdata2, nboot=10)
t_mle
