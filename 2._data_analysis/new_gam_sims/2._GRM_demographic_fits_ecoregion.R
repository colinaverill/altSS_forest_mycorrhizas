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
setnames(d1,'plot.BASAL','BASAL.plot')
setnames(d2,'plot.BASAL','BASAL.plot')

#grab ecoregions of interest.----
ecoregion.ref <- table(d1$ecoregion)
ecoregion.ref <- names(ecoregion.ref[ecoregion.ref > 1000])

#Fit growth, recruitment and mortality models by ecoregion.----
output <- list()
tic()
for(i in 1:length(ecoregion.ref)){
  #subset data to ecoregion of interest.
  d1.sub <- d1[d1$ecoregion == ecoregion.ref[i],]
  d2.sub <- d2[d2$ecoregion == ecoregion.ref[i],]
  
  #fit gam models.
  R.mod.em <- gam(recruit.em ~  s(relEM) + s(ndep, k = 3) + s(BASAL.plot, k = 3) + s(stem.density, k = 3)+
                    s(PC1, k=3) + s(PC2, k=3) + s(PC3, k=3) + s(PC4, k=3) + s(PC5, k=3) + s(PC6, k=3) + s(PC7, k=3) + s(PC8, k=3) + s(PC9, k=3) + s(PC10, k =3), 
                  data = d1.sub, family = 'poisson')
  R.mod.am <- gam(recruit.am ~  s(relEM) + s(ndep, k = 3) + s(BASAL.plot, k = 3) + s(stem.density, k = 3)+
                    s(PC1, k=3) + s(PC2, k=3) + s(PC3, k=3) + s(PC4, k=3) + s(PC5, k=3) + s(PC6, k=3) + s(PC7, k=3) + s(PC8, k=3) + s(PC9, k=3) + s(PC10, k =3), 
                  data = d1.sub, family = 'poisson')
  M.mod.em <- gam(mortality  ~  s(relEM) + s(ndep, k = 3) + s(BASAL.plot, k = 3) + s(stem.density, k = 3) + s(PREVDIA.cm, k = 5)+
                    s(PC1, k=3) + s(PC2, k=3) + s(PC3, k=3) + s(PC4, k=3) + s(PC5, k=3) + s(PC6, k=3) + s(PC7, k=3) + s(PC8, k=3) + s(PC9, k=3) + s(PC10, k =3), 
                  data = d2.sub[d2.sub$em == 1,], family = 'binomial')
  M.mod.am <- gam(mortality  ~  s(relEM) + s(ndep, k = 3) + s(BASAL.plot, k = 3) + s(stem.density, k = 3) + s(PREVDIA.cm, k = 5)+
                    s(PC1, k=3) + s(PC2, k=3) + s(PC3, k=3) + s(PC4, k=3) + s(PC5, k=3) + s(PC6, k=3) + s(PC7, k=3) + s(PC8, k=3) + s(PC9, k=3) + s(PC10, k =3), 
                  data = d2.sub[d2.sub$em == 0,], family = 'binomial')
  G.mod.em <- gam(DIA.cm   ~  s(relEM) +  s(ndep, k = 3) + s(BASAL.plot, k = 3) + s(stem.density, k = 3) + s(PREVDIA.cm, k = 5)+
                    s(PC1, k=3) + s(PC2, k=3) + s(PC3, k=3) + s(PC4, k=3) + s(PC5, k=3) + s(PC6, k=3) + s(PC7, k=3) + s(PC8, k=3) + s(PC9, k=3) + s(PC10, k =3), 
                  data = d2.sub[d2.sub$em == 1,])
  G.mod.am <- gam(DIA.cm   ~  s(relEM) + s(ndep, k = 3) + s(BASAL.plot, k = 3) + s(stem.density, k = 3) + s(PREVDIA.cm, k = 5)+
                    s(PC1, k=3) + s(PC2, k=3) + s(PC3, k=3) + s(PC4, k=3) + s(PC5, k=3) + s(PC6, k=3) + s(PC7, k=3) + s(PC8, k=3) + s(PC9, k=3) + s(PC10, k =3), 
                  data = d2.sub[d2.sub$em == 0,])
  
  #Grab plot environmental covariates for reference.
  #all plot-level environmental covariates. we will sample this later when simulating forests.
  all.cov <- d1.sub[,c('ndep','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')]
  all.cov <- d1.sub[,c('BASAL.plot','stem.density','ndep','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')]
  #mean plot-level environmental covariates.
  cov <- colMeans(all.cov)
  cov <- c(mean(d2.sub$PREVDIA.cm, na.rm=T),cov)
  names(cov)[1] <- 'PREVDIA.cm'
  
  #store in list and report.
  region.out <- list(G.mod.am, G.mod.em, M.mod.am, M.mod.em, R.mod.am, R.mod.em, cov)
  names(region.out) <- c('G.mod.am','G.mod.em','M.mod.am','M.mod.em','R.mod.am','R.mod.em','cov')
  output[[i]] <- region.out
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
