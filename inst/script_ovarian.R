

library("yppe")


data(ovarian)
ovarian <- na.exclude(ovarian)
ovarian$ESA <- as.factor(ovarian$ESA)
ovarian$EpoR <- as.factor(ovarian$EpoR)
ovarian$EphB4 <- as.factor(ovarian$EphB4)
ovarian$cytored <- as.factor(ovarian$cytored)
ovarian$age <- as.numeric(scale(ovarian$age))




mle <- yppe(Surv(time, status)~ESA+cytored|EpoR+age,
               data=ovarian, approach="mle")
summary(mle)


newdata <- data.frame(ESA=0:1, EpoR=rep(0,2))
St <- survfit(mle, newdata)
St

newdata0 <- data.frame(ESA=rep(0,2), EpoR=factor(0:1))
newdata1 <- data.frame(ESA=rep(1,2), EpoR=factor(1:0))
newdata0 <- data.frame(ESA=0, EpoR=0)
newdata1 <- data.frame(ESA=1, EpoR=0)
t <- crossSurv(mle, newdata0, newdata1)
t


newdata <- data.frame(ESA=factor(0:1), EpoR=factor(0:1))

St <- survfit(mle, newdata)
St




newdata1 <- data.frame(ESA=0, EpoR=0)
newdata2 <- data.frame(ESA=1, EpoR=0)







mle <- yppe(Surv(time, status)~ESA+cytored|EphB4+EpoR+age,
            data=ovarian, approach="mle")
summary(mle)




newdata0 <- data.frame(ESA=factor(0:1), age=c(-1,1))
newdata1 <- data.frame(ESA=factor(1:0), age=c(0,1))
newdata0
newdata1


newdata0 <- data.frame(ESA=0, cytored=0, EphB4=0, EpoR=0, age=0)
newdata1 <- data.frame(ESA=1, cytored=0, EphB4=0, EpoR=0, age=0)


newdata0 <- data.frame(ESA=0, cytored=0, EphB4=0, EpoR=0, age=0)
newdata1 <- data.frame(ESA=0, cytored=1, EphB4=0, EpoR=0, age=0)


t <- crossSurv(mle, newdata0, newdata1)
t

newdata <- data.frame(ESA=0, cytored=0, EphB4=0, EpoR=0, age=0)
St <- survfit(mle, newdata)



bayes <- yppe(Surv(time, status)~ESA+cytored|EphB4+EpoR+age,
            data=ovarian, approach="bayes")
summary(bayes)


newdata <- data.frame(ESA=0, cytored=0, EphB4=0, EpoR=0, age=0)
St <- survfit(bayes, newdata)
t <- crossSurv(bayes, newdata0, newdata1)
t
