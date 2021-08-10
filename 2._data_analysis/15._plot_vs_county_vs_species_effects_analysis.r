#Fitting Atlantic highlands ecoregion w/ plot and species random effects.
rm(list=ls())
library(data.table)
library(mgcv)
library(doParallel)
source('paths.r')
source('project_functions/tic_toc.r')

#set output path.----
output.path <- atlantic_highland_RE_test.path

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

#get genus species name
d2$genus_spp <- paste0(d2$GENUS,'_',d2$SPECIES)
d2$genus_spp <- as.factor(d2$genus_spp)

#grab Atlantic Highland ecoregion subset.----
ecoregion.ref <- table(d1$ecoregion)
ecoregion.ref <- names(ecoregion.ref[ecoregion.ref > 1000])

#set k and basis for spline.----
kk <- 5    #number of knots in smooth functions.
bs <- 'tp' #thin plate regression are default, but you can change this here.
nt <- 8    #number of cores to use while fitting.

#Fit growth, recruitment and mortality models by ecoregion.----
tic()
#subset data to ecoregion of interest (Atlantic Highlands)
d1.sub <- d1[d1$ecoregion == ecoregion.ref[1],]
d2.sub <- d2[d2$ecoregion == ecoregion.ref[1],]

#fit gam models - no random effect.
M.mod.em <- bam(mortality  ~  s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                data = d2.sub[d2.sub$em == 1,], family = 'binomial', nthreads=nt, discrete=T)
M.mod.am <- bam(mortality  ~  s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                data = d2.sub[d2.sub$em == 0,], family = 'binomial', nthreads=nt, discrete=T)

#Grab plot environmental covariates for reference.
#all plot-level environmental covariates. we will sample this later when simulating forests.
all.cov <- d1.sub[,c('BASAL.plot','BASAL.am','BASAL.em','stem.density','ndep','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')]
#mean plot-level environmental covariates.
cov <- colMeans(all.cov)
cov <- c(mean(d2.sub$PREVDIA.cm, na.rm=T),cov)
names(cov)[1] <- 'PREVDIA.cm'

#store in list and report.
region.out <- list(M.mod.am, M.mod.em, all.cov, cov)
names(region.out) <- c('M.mod.am','M.mod.em','all.cov','mean.cov')

#fit gam models with county level random effect.----
M.mod.em <- bam(mortality  ~  s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                  s(county.ID, bs = 're') + 
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                data = d2.sub[d2.sub$em == 1,], family = 'binomial', nthreads=nt, discrete=T)
M.mod.am <- bam(mortality  ~  s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                  s(county.ID, bs = 're') + 
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                data = d2.sub[d2.sub$em == 0,], family = 'binomial', nthreads=nt, discrete=T)
#  G.mod.em <- bam(DIA.cm   ~ s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) + 
#                    s(county.ID, bs = 're') +
#                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
#                  data = d2[d2$em == 1,], nthreads=nt, discrete=T)
#  G.mod.am <- bam(DIA.cm   ~ s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) + 
#                    s(county.ID, bs = 're') +
#                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
#                  data = d2[d2$em == 0,], nthreads=nt, discrete=T)

region.out.countyRE <- list(M.mod.am, M.mod.em, all.cov, cov)
names(region.out.countyRE) <- c('M.mod.am','M.mod.em','all.cov','mean.cov')
#region.out.countyRE <- list(M.mod.am, M.mod.em, G.mod.am,G.mod.em, all.cov, cov)
#names(region.out.countyRE) <- c('M.mod.am','M.mod.em','G.mod.am','G.mod.em','all.cov','mean.cov')

#fit gam models with plot level random effect.----
M.mod.em <- bam(mortality  ~  s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                  s(PLT_CN, bs = 're') + 
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                data = d2.sub[d2.sub$em == 1,], family = 'binomial', nthreads=nt, discrete=T)
M.mod.am <- bam(mortality  ~  s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                  s(PLT_CN, bs = 're') + 
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                data = d2.sub[d2.sub$em == 0,], family = 'binomial', nthreads=nt, discrete=T)
#  G.mod.em <- bam(DIA.cm   ~ s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) + 
#                    s(PLT_CN, bs = 're') +
#                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
#                  data = d2[d2$em == 1,], nthreads=nt, discrete=T)
#  G.mod.am <- bam(DIA.cm   ~ s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) + 
#                    s(PLT_CN, bs = 're') +
#                    s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k =kk, bs=bs), 
#                  data = d2[d2$em == 0,], nthreads=nt, discrete=T)

region.out.plotRE <- list(M.mod.am, M.mod.em, all.cov, cov)
names(region.out.plotRE) <- c('M.mod.am','M.mod.em','all.cov','mean.cov')
#region.out.plotRE <- list(M.mod.am, M.mod.em, G.mod.am,G.mod.em, all.cov, cov)
#names(region.out.plotRE) <- c('M.mod.am','M.mod.em','G.mod.am','G.mod.em','all.cov','mean.cov')

#fit gam models with species and county level random effects.
M.mod.em <- bam(mortality  ~  s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                  s(county.ID, bs = 're') + s(genus_spp, bs = 're') +
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                data = d2.sub[d2.sub$em == 1,], family = 'binomial', nthreads=nt, discrete=T)
M.mod.am <- bam(mortality  ~  s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs) +
                  s(county.ID, bs = 're') + s(genus_spp, bs = 're') +
                  s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                data = d2.sub[d2.sub$em == 0,], family = 'binomial', nthreads=nt, discrete=T)

region.out.countyRE.sppRE <- list(M.mod.am, M.mod.em, all.cov, cov)
names(region.out.countyRE.sppRE) <- c('M.mod.am','M.mod.em','all.cov','mean.cov')


#wrap it all.
total.out <- list(region.out, region.out.countyRE, region.out.plotRE,region.out.countyRE.sppRE)
names(total.out) <- c('none.re','county.re','plot.re','county.spp.re')
output <- total.out
toc()


#Save output.----
saveRDS(output, output.path, version = 2)
cat('Output saved, script complete.\n')

