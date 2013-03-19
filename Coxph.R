#
# Cox proportional hazard analysis of the OIF data
#
library(survival)
library(KMsurv)
library(MASS)
library(car) # for logit transformation
#
# read the data
data <- read.table(file="/Users/tselhorst/Desktop/OIF/data.txt", header=TRUE, sep="\t")
data
#
# List of Acronyms
#
# Country:____________________________________________________ Country
# Year ORV began:_____________________________________________ OIFS
# Cases at start of ORV:______________________________________ NRC
# Cases in 2010:______________________________________________ C_2010
# Last rabies case (except bats, humans and imported cases)___ LAST_CASE_Y
# Reduction of rabies:________________________________________ R_R
# ORV campaigns to eliminate rabies:__________________________ ORV_2_C100
# ORV campaigns to control* rabies (90%):_____________________ ORV_2_C90
# Size of territory (km2):____________________________________ ST
# Border length with endemic areas (km):______________________ BL
# Vaccinated area (km2):______________________________________ SVA
# Proportion of territory vaccinated:_________________________ TEA
# Mean Area Index:____________________________________________ AI
# campaigns until 2010:_______________________________________ Y <== dependent variable
# ER (ORV finished / not finished):___________________________ ER
#
#
# reshape the data so that we have a column Censor conditioned on ER which 
# tells us wheter the OIF is finishes or not. 
# When OIF is finished Censor is TRUE which tells me that the event was 
# really observed
data4surv <- transform(data, Censor=ifelse(ER=="finished", TRUE, FALSE))
data4surv$experience <- data$OIFS-1978
data4surv
attach(data4surv)
#
# create surv object with censored data
#
so <- Surv(Y, Censor, type='right')
so


#Create censoring for control
Censor2=ifelse(!is.na(ORV_2_C90), TRUE, FALSE)
YC=replace(ORV_2_C90, which(is.na(ORV_2_C90)), Y[which(is.na(ORV_2_C90))])
s2 <- Surv(Y, Censor2, type='right')
s2

#
# null model
null.model <- coxph(so~1)
summary(null.model)
# we alread know form negative binomial regression that onyl AI, TEA and log(SVA)
# are not correlated and can be used in the model. So let's start with 
# this one.
s.model <- coxph(so~AI+TEA+log(SVA))
summary(s.model)
#
# fit Coxph model in both directions. This time I want to restrict model
# parameters within interval [0,1]. This is achived with the logit transformation.
# interactions of order 2 are allowed
#
Scope = list(upper=~logit(AI)+logit(TEA)+log(SVA)+.^2,lower=~1)
phm_e=coxph(so~1, method="efron")
phm_f=stepAIC(phm_e,Scope,direction="both")
summary(phm_f)
plot(survfit(phm_f))
#
# plots for diferent figures of AI (min and max)
#
OIF.Ai <- data.frame(AI=c(min(AI),max(AI)),TEA=rep(median(TEA)))
plot(survfit(phm_f, newdata=OIF.Ai), conf.int=TRUE, lty=c(1,2),
     xlab="number of Campaigns", ylab="p(eradication not finsihed)")
legend(locator(1), legend=c('AI=0.016', 'AI=1.00'), lty=c(1,2))
#
# plots for diferent figures of TEA (min and max)
OIF.Tea <- data.frame(TEA=c(min(TEA),max(TEA)),AI=rep(median(AI),2),
                      xlab="number of Campaigns", ylab="p(eradication not finsihed)")
plot(survfit(phm_f, newdata=OIF.Tea), conf.int=TRUE, lty=c(1,2))
legend(locator(1), legend=c('TEA=0.016', 'TEA=1.00'), lty=c(1,2))
summary(survfit(phm_f, newdata=OIF.Tea))




###########################################################
#REPEAT FOR CONTROL DATA
## null model
null.model.control <- coxph(s2~1)
summary(null.model.control)
s.model.control <- coxph(s2~AI+TEA+log(SVA))
summary(s.model.control)
#
# fit Coxph model in both directions. This time I want to restrict model
# parameters within interval [0,1]. This is achived with the logit transformation.
# interactions of order 2 are allowed
#
phm_e_control=coxph(s2~1, method="efron")
phm_f_control=stepAIC(phm_e_control,Scope,direction="both")
summary(phm_f_control)
plot(survfit(phm_f_control))
#
# plots for diferent figures of AI (min and max)
#
plot(survfit(phm_f_control, newdata=OIF.Ai), conf.int=TRUE, lty=c(1,2),
     xlab="number of Campaigns", ylab="p(eradication not finsihed)")
legend(locator(1), legend=c('AI=0.016', 'AI=1.00'), lty=c(1,2))
#
# plots for diferent figures of TEA (min and max)
OIF.Tea <- data.frame(TEA=c(min(TEA),max(TEA)),AI=rep(median(AI),2),
                      xlab="number of Campaigns", ylab="p(eradication not finsihed)")
plot(survfit(phm_f_control, newdata=OIF.Tea), conf.int=TRUE, lty=c(1,2))
legend(locator(1), legend=c('TEA=0.016', 'TEA=1.00'), lty=c(1,2))
summary(survfit(phm_f_control, newdata=OIF.Tea))

#Compare elimination versus control
plot(survfit(phm_f), xlim=c(0,70), conf.int=TRUE)
lines(survfit(phm_f_control), col="red", conf.int=TRUE)

plot(survfit(phm_f, newdata=OIF.Ai), conf.int=TRUE, lty=c(1,2),
     xlab="number of Campaigns", ylab="p(rabies not eliminated")
legend(locator(1), legend=c('AI=0.016', 'AI=1.00'), lty=c(1,2))

lines(survfit(phm_f_control, newdata=OIF.Ai), conf.int=TRUE, lty=c(1,2), col="red")

plot(survfit(phm_f, newdata=OIF.Tea), conf.int=TRUE, lty=c(1,2))
legend(locator(1), legend=c('TEA=0.016', 'TEA=1.00'), lty=c(1,2))
lines(survfit(phm_f_control, newdata=OIF.Tea), conf.int=TRUE, lty=c(1,2), col="red")
