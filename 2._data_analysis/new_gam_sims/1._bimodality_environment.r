rm(list=ls())
source('paths.r')
source('project_functions/crib_fun.r')
library(mgcv)
library(boot)
library(betareg)

#load and prepare data.-----
d1 <- data.table(readRDS(Product_1.path))
d2 <- data.table(readRDS(Product_2.subset.path))

d1 <- d1[d1$PLT_CN %in% d2$PLT_CN,]
d1$y <- d1$relEM
d <- d1[,c('y','ndep','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')]
d <- d[complete.cases(d),]
d$y <- crib_fun(d$y)
#d$y <- ifelse(d$y == 0, min(d[d$y>0,]$y), d$y) #you can also do this, resids also look reasonable, results the same.
#d$y <- ifelse(d$y == 1, max(d[d$y<1,]$y), d$y)

#Fit gam to cribari-transformed, then logit transformed data.----
kk <- 5
m <- gam(logit(y) ~ s(ndep, k = kk) + s(PC1, k=kk) + s(PC2, k=kk) + s(PC3, k=kk) + s(PC4, k=kk) + s(PC5, k=kk) + s(PC6, k=kk) + s(PC7, k=kk) + s(PC8, k=kk) + s(PC9, k=kk) + s(PC10, k =kk),
         data = d, method = 'REML')

#plot raw and detrended relative abundance data.----
#grab mean covariates.
all.cov <- d1[,c('ndep','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')]
all.cov <- data.frame(t(colMeans(all.cov)))

par(mfrow = c(1,2))
#raw data.
hist(d$y, breaks = 20, ylim = c(0,1300))
#detrended m2.
mean.pred <- predict(m, newdata=all.cov)
detrend <- logit(d$y) - fitted(m) + rep(mean.pred, length(d$y))
hist(inv.logit(detrend), breaks = 20, ylim = c(0,1300))