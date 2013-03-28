#
# Cox proportional hazard analysis of the OIF data
#
library(gof)
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
# Survival analysis is performed for eradication of rabies, and for control
# rabies is reduced to 10% of the rabies cases at start of OIF.
#
# Create Censoring for eraducation:
# reshape the data so that we have a column Censor conditioned on ER which 
# tells us wheter the OIF is finishes or not. 
# When OIF is finished Censor is TRUE which tells me that the event was 
# really observed
data4surv <- transform(data, Censor=ifelse(ER=="finished", TRUE, FALSE))
attach(data4surv)
#
# create surv object with censored data for eradi
#
so <- Surv(Y, Censor, type='right')
so


#Create censoring for control
Censor2=ifelse(!is.na(ORV_2_C90), TRUE, FALSE)
YC=replace(ORV_2_C90, which(is.na(ORV_2_C90)), Y[which(is.na(ORV_2_C90))])
s2 <- Surv(YC, Censor2, type='right') # <--- bug reported from Katie fixed.
s2

#
# null model
null.model <- coxph(so~1)
summary(null.model)
# we alread know form negative binomial regression that onyl AI, TEA and log(SVA)
# are not correlated and can be used in the model. So let's start with 
# this one.
s.model <- coxph(so~logit(AI)+logit(TEA)+log(SVA))
summary(s.model)
#
# The following is all for eradication.
#
# fit Coxph model in both directions. This time I want to restrict model
# parameters within interval [0,1]. This is achived with the logit transformation.
# interactions of order 2 are allowed
#
Scope = list(upper=~logit(AI)+logit(TEA)+log(SVA)+.^2,lower=~1)
phm_e=coxph(so~1, method="efron")
phm_f=stepAIC(phm_e,Scope,direction="both")
summary(phm_f)
#
# model diagnostics
#
# a) test of the proportional hazard assumption
model.diagnostics <- cox.zph(phm_f)
model.diagnostics
#
# b) plots
plot(model.diagnostics)
#
# c) influence plots
dfbeta <- residuals(phm_f, type = 'dfbeta')
dfbeta
par(mfrow=c(1,2))
for(j in 1:2){
  plot(dfbeta[,j], ylab=names(coef(phm_f))[j])
  abline(h=0,lty=2)
}
system.time(phm_f.gof <- cumres(phm_f,R=50000))
summary(phm_f.gof)
plot(phm_f.gof)
par(mfrow=c(1,1))
plot(survfit(phm_f))
#
# plots for diferent figures of AI (min and max)
# perhaps we should exclude confidence intercals to make the plot more clear
#
par(mfrow=c(1,2))
OIF.Ai <- data.frame(AI=quantile(AI, c(0.1,0.9)),TEA=rep(median(TEA)))
plot(survfit(phm_f, newdata=OIF.Ai), conf.int=FALSE, lty=c(1,2),
     xlab="number of Campaigns", ylab="p(Time to eradication)", fun="event")
legend(locator(1), legend=c(quantile(AI, c(0.1,0.9))), 
       title="AI", 
       lty=c(1,2), 
       cex=0.75, 
       box.lwd=0, 
       bg="transparent")
#
# plots for diferent figures of TEA (min and max)
OIF.Tea <- data.frame(TEA=quantile(TEA, c(0.1,0.9)),AI=rep(median(AI)))
plot(survfit(phm_f, newdata=OIF.Tea), conf.int=FALSE, lty=c(2,1),
     xlab="number of Campaigns", ylab="p(Time to eradication)", fun="event")
legend(locator(1), legend=c(quantile(TEA, c(0.1,0.9))), 
       title="TEA",
       lty=c(2,1), 
       cex=0.75, 
       box.lwd=0, 
       bg="transparent")
summary(survfit(phm_f, newdata=OIF.Tea))

survfit(phm_f, newdata=OIF.Ai)
survfit(phm_f, newdata=OIF.Tea)
survfit(phm_f)


###########################################################
# REPEAT FOR CONTROL DATA
## null model
null.model.control <- coxph(s2~1)
summary(null.model.control)
s.model.control <- coxph(s2~logit(AI)+logit(TEA)+log(SVA))
summary(s.model.control)
#
# model diagnostics
#
# a) test of the proportional hazard assumption
model.diagnostics.control <- cox.zph(s.model.control)
model.diagnostics.control
#
# b) plots
plot(model.diagnostics.control)
#
# c) influence plots
dfbeta.control <- residuals(s.model.control, type = 'dfbeta')
dfbeta.control
par(mfrow=c(1,2))
for(j in 1:2){
  plot(dfbeta.control[,j], ylab=names(coef(s.model.control))[j])
  abline(h=0,lty=2)
}
system.time(s.model.control.gof <- cumres(s.model.control,R=50000))
summary(s.model.control.gof)
plot(s.model.control.gof)
par(mfrow=c(1,1))
plot(survfit(s.model.control))
#
# fit Coxph model in both directions. This time I want to restrict model
# parameters within interval [0,1]. This is achived with the logit transformation.
# interactions of order 2 are allowed
#
phm_e_control=coxph(s2~1, method="efron")
phm_f_control=stepAIC(phm_e_control,Scope,direction="both")
summary(phm_f_control)
par(mfrow=c(1,1))
plot(survfit(phm_f_control))
survfit(phm_f_control)
#
# model diagnostics for stepwise parameter selection for control
#
# a) test of the proportional hazard assumption
model.diagnostics.control.step <- cox.zph(phm_f_control)
model.diagnostics.control.step
#
# b) plots
plot(model.diagnostics.control.step)
#
# c) influence plots
dfbeta.control.step <- residuals(s.model.control.step, type = 'dfbeta')
dfbeta.control.step
par(mfrow=c(1,2))
for(j in 1:2){
  plot(dfbeta.control.step[,j], ylab=names(coef(s.model.control.step))[j])
  abline(h=0,lty=2)
}
system.time(s.model.control.step.gof <- cumres(phm_f_step,R=50000))
summary(s.model.control.step.gof)
plot(s.model.control.step.gof)
par(mfrow=c(1,1))
plot(survfit(s.model.control.step))
###################
# plots not fixed yet
###################
OIF.Ai.control <- data.frame(AI=c(min(AI),max(AI)),TEA=rep(mean(TEA)),SVA=rep(mean(SVA)))
print(survfit(phm_f_control, newdata=OIF.Ai.control), show.rmean=TRUE)
OIF.Tea.control <- data.frame(TEA=c(min(TEA),max(TEA)),AI=rep(mean(AI)),SVA=rep(mean(SVA)))
print(survfit(phm_f_control, newdata=OIF.Tea.control), show.rmean=TRUE)
#
# plots for diferent figures of AI (min and max)
#
par(mfrow=c(1,2))
plot(survfit(phm_f_control, newdata=OIF.Ai.control), conf.int=TRUE, lty=c(1,2),
     xlab="number of Campaigns", ylab="p(eradication not finsihed)")
legend(locator(1), legend=c('AI=0.016', 'AI=1.00'), lty=c(1,2), cex=0.5)
#
# plots for diferent figures of TEA (min and max)
#OIF.Tea <- data.frame(TEA=c(min(TEA),max(TEA)),AI=rep(median(AI),2))
plot(survfit(phm_f_control, newdata=OIF.Tea.control), conf.int=TRUE, lty=c(1,2),
     xlab="number of Campaigns", ylab="p(eradication not finsihed)")
legend(locator(1), legend=c('TEA=0.016', 'TEA=1.00'), lty=c(1,2), cex=0.5)
summary(survfit(phm_f_control, newdata=OIF.Tea.control))

par(mfrow=c(1,1))
#Compare elimination versus control
plot(survfit(phm_f), xlim=c(0,70), conf.int=TRUE)
lines(survfit(phm_f_control), col="red", conf.int=TRUE)

plot(survfit(phm_f, newdata=OIF.Ai), conf.int=TRUE, lty=c(1,2),
     xlab="number of Campaigns", ylab="p(rabies not eliminated")
legend(locator(1), legend=c('AI=0.016', 'AI=1.00'), lty=c(1,2))

lines(survfit(phm_f_control, newdata=OIF.Ai.control), conf.int=TRUE, lty=c(1,2), col="red")

plot(survfit(phm_f, newdata=OIF.Tea), conf.int=TRUE, lty=c(1,2))
legend(locator(1), legend=c('TEA=0.016', 'TEA=1.00'), lty=c(1,2))
lines(survfit(phm_f_control, newdata=OIF.Tea.control), conf.int=TRUE, lty=c(1,2), col="red")