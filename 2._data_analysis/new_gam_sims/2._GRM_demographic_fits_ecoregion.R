#Fitting gams w/ separate AM-EM models and PC scores across ecoregions.
rm(list=ls())
library(data.table)
library(mgcv)
library(doParallel)
source('paths.r')
source('project_functions/tic_toc.r')

#set output path.----
output.path <- demographic_fits_gam_ecoregion.path

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

#grab ecoregions of interest.----
ecoregion.ref <- table(d1$ecoregion)
ecoregion.ref <- names(ecoregion.ref[ecoregion.ref > 1000])

#parallel environment.-----
cl <- makeCluster(28)

#set k and basis for spline.----
kk <- 5
bs <- 'cr' #thin plate regression splines too slow, and results are similar.

#Fit growth, recruitment and mortality models by ecoregion.----
output <- list()
tic()
for(i in 1:length(ecoregion.ref)){
  #subset data to ecoregion of interest.
  d1.sub <- d1[d1$ecoregion == ecoregion.ref[i],]
  d2.sub <- d2[d2$ecoregion == ecoregion.ref[i],]
  
  #fit gam models - no random effect.
  R.mod.em <- bam(recruit.em ~  s(relEM, k=kk, bs=bs) + s(BASAL.em, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) +
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d1.sub, family = 'poisson')
  R.mod.am <- bam(recruit.am ~  s(relEM, k=kk, bs=bs) + s(BASAL.am, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) +
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d1.sub, family = 'poisson')
  M.mod.em <- bam(mortality  ~  s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d2.sub[d2.sub$em == 1,], family = 'binomial', cluster = cl)
  M.mod.am <- bam(mortality  ~  s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d2.sub[d2.sub$em == 0,], family = 'binomial', cluster = cl)
  G.mod.em <- bam(DIA.cm   ~  s(relEM, k=kk, bs=bs) +  s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d2.sub[d2.sub$em == 1,], cluster = cl)
  G.mod.am <- bam(DIA.cm   ~  s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d2.sub[d2.sub$em == 0,], cluster = cl)
  
  #Grab plot environmental covariates for reference.
  #all plot-level environmental covariates. we will sample this later when simulating forests.
  all.cov <- d1.sub[,c('BASAL.plot','BASAL.am','BASAL.em','stem.density','ndep','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')]
  #mean plot-level environmental covariates.
  cov <- colMeans(all.cov)
  cov <- c(mean(d2.sub$PREVDIA.cm, na.rm=T),cov)
  names(cov)[1] <- 'PREVDIA.cm'
  
  #store in list and report.
  region.out <- list(G.mod.am, G.mod.em, M.mod.am, M.mod.em, R.mod.am, R.mod.em, cov)
  names(region.out) <- c('G.mod.am','G.mod.em','M.mod.am','M.mod.em','R.mod.am','R.mod.em','cov')
  
  #fit gam models with county level random effect.----
  M.mod.em <- bam(mortality  ~  s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(county.ID, bs = 're') + 
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d2.sub[d2.sub$em == 1,], family = 'binomial', cluster = cl)
  M.mod.am <- bam(mortality  ~  s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(county.ID, bs = 're') + 
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d2.sub[d2.sub$em == 0,], family = 'binomial', cluster = cl)
  G.mod.em <- bam(DIA.cm   ~  s(relEM, k=kk, bs=bs) +  s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(county.ID, bs = 're') + 
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d2.sub[d2.sub$em == 1,], cluster = cl)
  G.mod.am <- bam(DIA.cm   ~  s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(county.ID, bs = 're') + 
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d2.sub[d2.sub$em == 0,], cluster = cl)
  
  region.out.countyRE <- list(G.mod.am, G.mod.em, M.mod.am, M.mod.em, R.mod.am, R.mod.em, cov)
  names(region.out.countyRE) <- c('G.mod.am','G.mod.em','M.mod.am','M.mod.em','R.mod.am','R.mod.em','cov')
  
  #fit gam models with plot level random effect.----
  M.mod.em <- bam(mortality  ~  s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(PLT_CN, bs = 're') + 
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d2.sub[d2.sub$em == 1,], family = 'binomial', cluster = cl)
  M.mod.am <- bam(mortality  ~  s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(PLT_CN, bs = 're') + 
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d2.sub[d2.sub$em == 0,], family = 'binomial', cluster = cl)
  G.mod.em <- bam(DIA.cm   ~  s(relEM, k=kk, bs=bs) +  s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(PLT_CN, bs = 're') + 
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d2.sub[d2.sub$em == 1,], cluster = cl)
  G.mod.am <- bam(DIA.cm   ~  s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                    s(PLT_CN, bs = 're') + 
                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d2.sub[d2.sub$em == 0,], cluster = cl)
  
  region.out.plotRE <- list(G.mod.am, G.mod.em, M.mod.am, M.mod.em, R.mod.am, R.mod.em, cov)
  names(region.out.plotRE) <- c('G.mod.am','G.mod.em','M.mod.am','M.mod.em','R.mod.am','R.mod.em','cov')
  
  total.out <- list(region.out, region.out.countyRE, region.out.plotRE)
  names(total.out) <- c('none.re','county.re','plot.re')
  output[[i]] <- 
    msg <- paste0(ecoregion.ref[i],' fit. ',i,' of ',length(ecoregion.ref),' ecoregions complete.\n')
  cat(msg); toc()
}
names(output) <- ecoregion.ref

#Save models and size categories.----
saveRDS(output, output.path, version = 2)
cat('Output saved, script complete.\n')

#Plotting.----
plotting <- F
if(plotting == T){
  #Plot recruitment models.
  par(mfrow = c(5,2))
  output <- output$none.re
  for(i in 1:length(output)){
    new.dat <- data.frame(seq(0 ,1 ,by = 0.01), t(output[[i]]$cov))
    colnames(new.dat)[1] <- 'relEM'
    pred.r.am <- predict(output[[i]]$R.mod.am, newdata = new.dat)
    pred.r.em <- predict(output[[i]]$R.mod.em, newdata = new.dat)
    plot(pred.r.am ~ new.dat$relEM, bty = 'l', ylim = c(min(c(pred.r.am,pred.r.em)), max(c(pred.r.am,pred.r.em))))
    plot(pred.r.em ~ new.dat$relEM, bty = 'l', ylim = c(min(c(pred.r.am,pred.r.em)), max(c(pred.r.am,pred.r.em))))
  }
  
  #Plot mortality models.
  par(mfrow = c(5,2))
  for(i in 1:length(output)){
    new.dat <- data.frame(seq(0 ,1 ,by = 0.01), t(output[[i]]$cov))
    colnames(new.dat)[1] <- 'relEM'
    pred.m.am <- predict(output[[i]]$M.mod.am, newdata = new.dat)
    pred.m.em <- predict(output[[i]]$M.mod.em, newdata = new.dat)
    plot(pred.m.am ~ new.dat$relEM, bty = 'l', ylim = c(min(c(pred.m.am,pred.m.em)), max(c(pred.m.am,pred.m.em))))
    plot(pred.m.em ~ new.dat$relEM, bty = 'l', ylim = c(min(c(pred.m.am,pred.m.em)), max(c(pred.m.am,pred.m.em))))
  }
  
  #Plot Growth models.
  par(mfrow = c(5,2))
  for(i in 1:length(output)){
    new.dat <- data.frame(seq(0 ,1 ,by = 0.01), t(output[[i]]$cov))
    colnames(new.dat)[1] <- 'relEM'
    pred.g.am <- predict(output[[i]]$G.mod.am, newdata = new.dat)
    pred.g.em <- predict(output[[i]]$G.mod.em, newdata = new.dat)
    plot(pred.g.am ~ new.dat$relEM, bty = 'l', ylim = c(min(c(pred.g.am,pred.g.em)), max(c(pred.g.am,pred.g.em))))
    plot(pred.g.em ~ new.dat$relEM, bty = 'l', ylim = c(min(c(pred.g.am,pred.g.em)), max(c(pred.g.am,pred.g.em))))
  }
}
