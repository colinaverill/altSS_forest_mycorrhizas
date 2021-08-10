rm(list=ls())
source('paths.r')
source('project_functions/crib_fun.r')
library(mgcv)
library(boot)
library(diptest)
library(betareg)
library(data.table)

#set output path.----
output.path <- bimodality_results.path

#load and prepare data.-----
d1 <- data.table(readRDS(Product_1.path))
d2 <- data.table(readRDS(Product_2.subset.path))

d1 <- d1[d1$PLT_CN %in% d2$PLT_CN,]
d1$y <- d1$relEM
d1$county.id <- paste0(d1$STATECD,'_',d1$COUNTYCD)
d1$county.id <- as.factor(d1$county.id)
d <- d1[,c('y','relEM','ndep','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','county.id')]
d <- d[complete.cases(d),]
d$y <- crib_fun(d$y)
#d$y <- ifelse(d$y == 0, min(d[d$y>0,]$y), d$y) #you can also do this, residuls also look reasonable, results the same.
#d$y <- ifelse(d$y == 1, max(d[d$y<1,]$y), d$y)

#Fit gam to cribari-transformed data.----
kk <- 5 #number of spline knots.
nt <- 8 #number of threads on machine.
m <- bam(y ~ s(ndep, k=kk) + s(PC1, k=kk) + s(PC2, k=kk) + s(PC3, k=kk) + s(PC4, k=kk) + s(PC5, k=kk) + s(PC6, k=kk) + s(PC7, k=kk) + s(PC8, k=kk) + s(PC9, k=kk) + s(PC10, k =kk)
          + s(county.id, bs='re')
         , data = d, family=binomial(link='logit'), nthreads=nt, discrete=T)

#plot raw and detrended relative abundance data.----
#grab mean covariates.
raw.cov <- d1[,c('ndep','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','county.id')]
all.cov <- data.frame(t(colMeans(raw.cov[,1:11])))
all.cov$county.id <- '9_11'

par(mfrow = c(1,2))
#raw data.
hist(d$y, breaks = 20, ylim = c(0,1300))
#detrended m2.
mean.pred <- predict(m, newdata=all.cov, exclude = 's(county.id)', newdata.guaranteed = T)
detrend <- logit(d$y) - logit(fitted(m)) + rep(mean.pred, length(d$y)) #predictions are done on logit scale, no need for logit transform.
hist(inv.logit(detrend), breaks = 20, ylim = c(0,1300))

#hartigan's dip test.----
raw.dip <- dip.test(d$y)
detrend.dip <- dip.test(detrend)

#Save output.----
output <- cbind(d1$latitude.y, d1$longitude.y, d$y, detrend, raw.cov)
colnames(output)[1:4] <- c('latitude','longitude','relEM','relEM.detrend')
output$detrend <- inv.logit(output$relEM.detrend)
output <- list(output,raw.dip,detrend.dip)
names(output) <- c('data','raw.dip.test','detrend.dip.test')
saveRDS(output, output.path)
