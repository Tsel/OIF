#
# How to visualize an interaction term?
# 
# Example taken form 
#   http://stats.stackexchange.com/questions/6557/how-can-one-plot-continuous-by-continuous-interactions-in-ggplot2
#
x1 <- rnorm(100,2,10)
x2 <- rnorm(100,2,10)
y <- x1+x2+x1*x2+rnorm(100,1,2)
dat <- data.frame(y=y, x1=x1, x2=x2)
res <- lm(y~x1*x2)
summary(res)
#
z1 <- z2 <- seq(-1,1)
#
# create a grid
newdf <- expand.grid(x1=z1,x2=z2)
#
# plot the interaction
tdat <- transform(newdf, yp=predict(res, newdf))
library(ggplot2)
p <- ggplot(tdat,
            aes(y=yp, x=x1, color=factor(x2))) + stat_smooth(method=lm)
p