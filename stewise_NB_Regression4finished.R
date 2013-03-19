library(car)
library(MASS)
library(ggplot2)
#
# Statistical analysis of the influence of ST, TEA, AI and NRC [explanation below] on the number of
# campaigns
#
# use only countries which finished the OIF
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
data <- read.table(file="/Users/tselhorst/Desktop/OIF/data.txt", header=TRUE, sep="\t")
attach(data)
data <- subset(data, ER=="finished")
#
# Correlatons and correlation matrix plot
#
library(Hmisc)
d <- data.frame(AI, TEA, BL, log(SVA), log(ST), log(NRC))
df <- data.matrix(d)
rcorr(df, type="spearman")
#
# corr matrix plot
#
library(lattice)
# log(ST) is correlated with TEA
scatterplotMatrix(~Y+AI+TEA+BL+log(SVA)+log(ST)+log(NRC) | ER, reg.line=FALSE, smoother=FALSE)
#
# Perform stepwise main effects nb regression model
#
# Parameters of the model are:
#   Y   : number of OIF campaigns
#   TEA : territory ever affected (proportion)
#   AI  : area index
#   SVA : size of vaccination areq <== element within the calc of AI
#
#
model.nb <- glm.nb(Y ~ (TEA + AI), data = data)
model.nb.step <- stepAIC(model.nb, direction="backward")
model.nb.step$anova
summary(model.nb.step)
with(model.nb.step, cbind(res.deviance = deviance, df=df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))
#
# include second oder interactions to see whether the better fit better
model.nb.step.2 <- stepAIC(model.nb.step,~.^2, direction="both")
model.nb.step.2$anova
summary(model.nb.step.2)
with(model.nb.step.2, cbind(res.deviance = deviance, df=df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))
#
# why not use the poisson model?
#
# setup a poission regression model with is nested in the negative binomial model and assuming that the disperion parameter is
# constant.
#
model.poi <- glm(Y ~ (log(SVA) + TEA + AI), family = poisson, data = data)
summary(model.poi)
model.poi.step = stepAIC(model.poi, direction="backward")
model.poi.step$anova
summary(model.poi.step)
#
# include second order interactions to see whether the model performs better
model.poi.step.2 <- stepAIC(model.poi.step,~.^2, direction="both")
model.poi.step.2$anova
summary(model.poi.step.2)
with(model.poi.step.2, cbind(res.deviance = deviance, df=df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))
#
# compare negbin with poi
#
chi <- 2 * (logLik(model.nb.step.2) - logLik(model.poi.step.2))
chi
pchisq(chi, df = 1, lower.tail = FALSE)
#
# Confidence intervals for the coefficients of the model
#
est <- cbind(Estimate=coef(model.nb.step.2), confint(model.nb.step.2))
est
#
#######################
# Grafical representation of the reslts
#######################
#
# we need a grid for the significant interaction term TEA x AI
# sequence for AI
z1 <- seq(min(AI),max(AI),length=100)
# seauence for TEA
z2 <- seq(min(TEA),max(TEA),length=3)
# create a grid
grid.data <- expand.grid(AI=z1,TEA=z2)
grid.data$SVA <- median(SVA)
grid.data <- cbind(grid.data, predict(model.nb.step.2, grid.data, type="response", se.fit=TRUE))
grid.data <- within(grid.data, {
  LL <- fit-1.96*se.fit
  UL <- fit+1.96*se.fit
})
#
cp <- c("blue","goldenrod","indianred")
ggplot(grid.data, aes(AI, fit, color=factor(TEA)))+
  geom_ribbon(aes(ymin=LL, ymax=UL, fill=factor(TEA)),alpha=0.2)+
  geom_line(aes(color=factor(TEA)),size=1)+
  labs(x="Area Index",y="OIF campaigns")+
  scale_colour_manual(values=cp,name="TEA",guide="none")+
  scale_fill_manual(values=cp,name="TEA")+
  geom_point(data=data, aes(AI, Y), color="black")+
  annotate("text", label="SVA=60056.5",x=0.75,y=100)
#
# graphical respresentation of the effect of SVA
#
sva <- seq(min(SVA),max(SVA),length=1000)
df.sva <- data.frame(SVA=sva)
df.sva$AI <- median(AI)
df.sva$TEA <- median(TEA)
names(df.sva)
#is.data.frame(df.sva)
#pred <- predict(model.nb.step.2, df.sva, type="response", se.fit=TRUE)
df.sva <- cbind(df.sva, predict(model.nb.step.2, df.sva, type="response", se.fit=TRUE))
df.sva <- within(df.sva, {
  LL <- fit-1.96*se.fit
  UL <- fit+1.96*se.fit
})
#
ggplot(df.sva, aes(SVA, fit), color="blue")+
  geom_ribbon(aes(ymin=LL,ymax=UL), color="blue", alpha=0.2)+
  geom_line(aes(),color="blue", size=1)+
  labs(x="size of vaccinated area [SVA]", y="OIF campaigns")+
  geom_point(data=data, aes(SVA, Y), color="black")