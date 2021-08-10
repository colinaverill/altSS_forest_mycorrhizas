#Checking effects of stand age on relative abundance of AM vs. EM trees.
#No correlation. Nice.
rm(list=ls())
source('paths.r')
library(data.table)
library(boot)
library(mgcv)
library(car)
source('project_functions/crib_fun.r')

#set output path.----
output.path <- stand_age_results.path

#load plot data and growth/mortality/recruitment data.----
d1 <- data.table(readRDS(Product_1.path))
d2 <- data.table(readRDS(Product_2.subset.path))

#Subset plot data and rename some things.-----
d1 <- d1[d1$PLT_CN %in% d2$PLT_CN,]
d2$PLT_CN <- as.factor(d2$PLT_CN)
d1$county.id <- paste0(d1$STATECD,'_',d1$COUNTYCD)
d1$county.id <- as.factor(d1$county.id)
d2$county.ID <- as.factor(paste0(d2$STATECD,'_',d2$COUNTYCD))
d1$county.ID <- as.factor(paste0(d1$STATECD,'_',d1$COUNTYCD))
setnames(d1,'plot.BASAL','BASAL.plot')
setnames(d2,'plot.BASAL','BASAL.plot')
setnames(d1,'BASAL.ECM','BASAL.em')
setnames(d1,'BASAL.AM' ,'BASAL.am')
d1$y <- crib_fun(d1$relEM)

#gam model.----
kk <- 5
nt <- 7
m <- bam(relEM ~ s(ndep, k=kk) + s(PC1, k=kk) + s(PC2, k=kk) + s(PC3, k=kk) + s(PC4, k=kk) + s(PC5, k=kk) + s(PC6, k=kk) + s(PC7, k=kk) + s(PC8, k=kk) + s(PC9, k=kk) + s(PC10, k =kk)
         + s(county.id, bs='re') + s(STDAGE, k=kk)
         , data = d1, family=binomial(link='logit'), nthreads=nt, discrete=T) 

#detrend TWO WAYS.----
#First way leaving in the effect of stand age so it can be visualized after controlling for other covariates.
raw.cov <- d1[,c('ndep','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','county.id')]
all.cov <- data.frame(t(colMeans(raw.cov[,1:11])))
all.cov$county.id <- '9_11'
all.cov <- data.frame(d1$STDAGE, all.cov)
colnames(all.cov)[1] <- 'STDAGE'
mean.pred <- predict(m, newdata=all.cov, exclude = 's(county.id)', newdata.guaranteed = T)
detrend1 <- logit(d1$y) - logit(fitted(m)) + mean.pred #predictions are done on logit scale, no need for logit transform.

#Second way showing the bimodality remains after controlling for stand age.
raw.cov <- d1[,c('ndep','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','STDAGE','county.id')]
all.cov <- data.frame(t(colMeans(raw.cov[,1:12])))
all.cov$county.id <- '9_11'
mean.pred <- predict(m, newdata=all.cov, exclude = 's(county.id)', newdata.guaranteed = T)
detrend2 <- logit(d1$y) - logit(fitted(m)) + mean.pred #predictions are done on logit scale, no need for logit transform.

#save output.----
output <- data.frame(detrend1, detrend2, d1$STDAGE)
colnames(output) <- c('detrend1','detrend2','STDAGE')
saveRDS(output, output.path)
