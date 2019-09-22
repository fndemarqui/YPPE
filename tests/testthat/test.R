
library(YPPE)

# ML approach:
mle <- yppe(Surv(time, status)~arm, data=ipass, approach="mle")
summary(mle)
newdata1 <- data.frame(arm=0)
newdata2 <- data.frame(arm=1)
tcross <- crossTime(mle, newdata1, newdata2, nboot=10)
tcross
ekm <- survfit(Surv(time, status)~arm, data=ipass)
newdata <- data.frame(arm=0:1)
St <- survfit(mle, newdata)
time <- sort(ipass$time)
plot(ekm, col=1:2)
lines(time, St[[1]])
lines(time, St[[2]], col=2)
abline(v=tcross, col="blue")

# Bayesian approach:
bayes <- yppe(Surv(time, status)~arm, data=ipass,
              approach="bayes", chains=1, iter=10)
summary(bayes)
newdata1 <- data.frame(arm=0)
newdata2 <- data.frame(arm=1)
tcross <- crossTime(bayes, newdata1, newdata2)
tcross
ekm <- survfit(Surv(time, status)~arm, data=ipass)
newdata <- data.frame(arm=0:1)
St <- survfit(bayes, newdata)
time <- sort(ipass$time)
plot(ekm, col=1:2)
lines(time, St[[1]])
lines(time, St[[2]], col=2)
abline(v=tcross, col="blue")

