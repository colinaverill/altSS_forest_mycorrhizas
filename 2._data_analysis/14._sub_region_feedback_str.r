#compare strength of feedbacks to ratio of mono-myco to mixed.
rm(list=ls())
library(data.table)
library(mgcv)
library(hexbin)
library(boot)
source('paths.r')
source('project_functions/predict_gam_well.r')

#set output path.----
output.path <- bimod_feedback_str_analysis.path

#load growth/mortality/recruitment data.----
d1 <- data.table(readRDS(Product_1.path))
d2 <- data.table(readRDS(Product_2.subset.path))

#Subset and rename some things.-----
d1 <- d1[d1$PLT_CN %in% d2$PLT_CN,]
d2$PLT_CN <- as.factor(d2$PLT_CN)
d2$county.ID <- as.factor(paste0(d2$STATECD,'_',d2$COUNTYCD))
d1$county.ID <- as.factor(paste0(d1$STATECD,'_',d1$COUNTYCD))
setnames(d1,'plot.BASAL','BASAL.plot')
setnames(d2,'plot.BASAL','BASAL.plot')
setnames(d1,'BASAL.ECM','BASAL.em')
setnames(d1,'BASAL.AM' ,'BASAL.am')

#Bin the FIA data into hexagons.----
hbin <- hexbin(d1$longitude.x, d1$latitude.x, xbins = 10, IDs = T) #xbins controls the number of hexagons.
plot(hbin) #looks cool!
d1$bin.ID <- paste0('bin_',hbin@cID)
#16 bins to fit.
check <- table(d1$bin.ID)
to_fit <- check[check > 150]
to_fit <- names(to_fit)

#grab latitude means.
lat.bin <- aggregate(latitude.x ~ bin.ID, data = d1, FUN = mean)

#Fit models to sub-regions.----
kk <- 5    #global max number of 'knots' to use in smoothing functions.
bs <- 'tp' #thin plate regression splines are default as well, but you can change things here if you like.
nt <- 8    #number of threads.

bin_fits <- list()
for(i in 1:length(to_fit)){
  #subset to hexagon.
  d1.sub <- d1[d1$bin.ID == to_fit[i],]
  d2.sub <- d2[d2$PLT_CN %in% d1.sub$PLT_CN,]
  
  #fit models.
  R.mod.em <- bam(recruit.em ~ s(relEM, k=kk, bs=bs) + s(BASAL.em, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) +
                    s(county.ID, bs = 're') +
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
                  data = d1.sub, family = 'poisson', nthreads=nt, discrete=T)
  R.mod.am <- bam(recruit.am ~ s(relEM, k=kk, bs=bs) + s(BASAL.am, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) +
                    s(county.ID, bs = 're') +
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
                  data = d1.sub, family = 'poisson', nthreads=nt, discrete=T)
  M.mod.em <- bam(mortality ~ s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) + 
                    s(county.ID, bs = 're') +
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
                  data = d2.sub[d2.sub$em == 1,], family = 'binomial', nthreads=nt, discrete=T)
  M.mod.am <- bam(mortality ~ s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) + 
                    s(county.ID, bs = 're') +
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
                  data = d2.sub[d2.sub$em == 0,], family = 'binomial', nthreads=nt, discrete=T)
  #get mean covariate set.
  all.cov <- d1.sub[,c('BASAL.plot','BASAL.am','BASAL.em','stem.density','ndep','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')]
  #mean plot-level environmental covariates.
  cov <- colMeans(all.cov)
  cov <- c(mean(d2.sub$PREVDIA.cm,na.rm=T), cov)
  names(cov)[1] <- 'PREVDIA.cm'
  
  out <- list(R.mod.em, R.mod.am, M.mod.em, M.mod.am, cov)
  names(out) <- c('R.mod.em','R.mod.am','M.mod.em','M.mod.am','mean.cov')
  bin_fits[[i]] <- out
}

#Get ratios of mono to mixed mycorrhizal stands across ecoregions.-----
mono.balance <- list()
for(i in 1:length(to_fit)){
  sub <- d1[d1$bin.ID == to_fit[i],]
  am <- nrow(sub[sub$relEM <= 0.01,])
  em <- nrow(sub[sub$relEM >= 0.99,])
  mix <- nrow(sub) - (am+em)
  mono.balance[[i]] <- c(am,em,mix,(am+em)/mix)
}
mono.balance <- data.frame(do.call(rbind, mono.balance))
rownames(mono.balance) <- to_fit
colnames(mono.balance) <- c('am','em','mix','mono_mix_ratio')
mono.balance$mono_frac <- (mono.balance$am + mono.balance$em) / (mono.balance$am + mono.balance$em + mono.balance$mix)
mono.balance$am_frac <- mono.balance$am / (mono.balance$am + mono.balance$em + mono.balance$mix)
mono.balance$em_frac <- mono.balance$em / (mono.balance$am + mono.balance$em + mono.balance$mix)

#Get an estimate of AM and EM recruitment feedback strength.----
advantage.calc <- list()
for(i in 1:length(to_fit)){
  sub.fit.am.r <- bin_fits[[i]]$R.mod.am
  sub.fit.em.r <- bin_fits[[i]]$R.mod.em
  sub.fit.am.m <- bin_fits[[i]]$M.mod.am
  sub.fit.em.m <- bin_fits[[i]]$M.mod.em
  
  #Build prediction dataset, add county level randome ffect that will be ignored.
  check1 <- sub.fit.am.r$model$county.ID
  check2 <- sub.fit.em.r$model$county.ID
  check <- check1[check1 %in% check2]
  cov <- bin_fits[[i]]$mean.cov
  relEM <- seq(0,1, by=1)
  dat <- data.frame(relEM, t(cov))
  dat$county.ID <- as.factor(check[10])
  
  #make predictions.
  am.recr <- predict_gam_well(sub.fit.am.r, newdata = dat, ranef.lab = "county.ID")$fit
  em.recr <- predict_gam_well(sub.fit.em.r, newdata = dat, ranef.lab = "county.ID")$fit
  am.mort <- predict_gam_well(sub.fit.am.m, newdata = dat, ranef.lab = "county.ID")$fit
  em.mort <- predict_gam_well(sub.fit.em.m, newdata = dat, ranef.lab = "county.ID")$fit
  
  #calculate advantange.
  em.advantage.r <- (em.recr[2]) - (em.recr[1])
  am.advantage.r <- (am.recr[1]) - (am.recr[2])
  em.advantage.m <- (em.mort[2] - em.mort[1]) * -1 #time minus 1 because less likely to die means more likely to survive.
  am.advantage.m <- (am.mort[1] - am.mort[2]) * -1
  advantage.calc[[i]] <- c(em.advantage.r, am.advantage.r,em.advantage.m,am.advantage.m)
}
advantage.calc <- data.frame(do.call(rbind, advantage.calc))
rownames(advantage.calc) <- to_fit
colnames(advantage.calc) <- c('em.recr.str','am.recr.str','em.mort.str','am.mort.str')
#calculate averages of these strengths or interactions across AM-EM groups.
advantage.calc$recr.str <- (advantage.calc$am.recr.str+advantage.calc$em.recr.str)/2
advantage.calc$mort.str <- (advantage.calc$am.mort.str+advantage.calc$em.mort.str)/2
#merge the two frames.
all.dat <- cbind(mono.balance, advantage.calc)
all.dat$bin.ID <- rownames(all.dat)

#save output.----
saveRDS(all.dat, output.path)

#fit and plot.-----
plot_it <- F
if(plot_it == T){
  par(mfrow = c(1,2))
  #mortality only fit.
  plot(logit(mono_frac) ~ mort.str, data = all.dat, 
       bty = 'l', pch = 16, cex=1.5,
       xlab = 'Survival Advantage',
       ylab = NA)
  fit.mort <- lm(logit(mono_frac) ~ mort.str, data = all.dat)
  rsq     <- round(summary(fit.mort)$adj.r.squared, 2)
  rsq.raw <- round(summary(fit.mort)$r.squared    , 2)
  rsq.lab <- bquote({R^2}[adj] == .(rsq))
  rsq.raw.lab <- bquote({R^2} == .(rsq.raw))
  abline(fit.mort, lwd = 2)
  mtext(rsq.raw.lab, side=3, line=-2.0, adj=1, at = 0.5)
  mtext(rsq.lab    , side=3, line=-3.5, adj=1, at = 0.5)
  mtext('P < 0.001', side=3, line=-4.5, adj=1, at = 0.5)
  mtext('Mycorrhizal Monodominance \n Metric',side = 2, line = 2)
  mtext('A', side = 1, line = -1.5, adj = 0.95)
  
  #recruitment only fit.
  plot(logit(mono_frac) ~ recr.str, data = all.dat, 
       bty = 'l', pch = 16, cex=1.5,
       xlab = 'Recruitment Advantage',
       ylab = NA)
  fit.recr <- lm(logit(mono_frac) ~ recr.str, data = all.dat)
  #abline(fit.recr, lwd = 1, lty=2, col = 'light grey')
  mtext('Mycorrhizal Monodominance \n Metric',side = 2, line = 2)
  mtext('N.S.', side=3, line=-2, adj=0.95)
  mtext('B', side = 1, line = -1.5, adj = 0.95)
}
