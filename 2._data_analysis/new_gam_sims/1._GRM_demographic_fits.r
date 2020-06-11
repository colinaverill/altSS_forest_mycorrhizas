#Fitting gams w/ separate AM-EM models and PC scores.
rm(list=ls())
library(data.table)
library(mgcv)
library(doParallel)
source('paths.r')
source('project_functions/tic_toc.r')

#set output path.----
output.path <- demographic_fits_gam_separate.path

#load growth/mortality/recruitment data.----
d1 <- data.table(readRDS(Product_1.path))
d2 <- data.table(readRDS(Product_2.subset.path))

#Subset and rename some things.-----
d1 <- d1[d1$PLT_CN %in% d2$PLT_CN,]
d2$PLT_CN <- as.factor(d2$PLT_CN)
setnames(d1,'plot.BASAL','BASAL.plot')
setnames(d2,'plot.BASAL','BASAL.plot')
setnames(d1, 'BASAL.ECM','BASAL.em')
setnames(d1, 'BASAL.AM' ,'BASAL.am' )

#Set global k.----
kk <- 5

#Fit growth, recruitment and mortality models.----
#Environmental models without feedbacks.
cat('Fitting null models...\n');tic()
R.mod.em <- gam(recruit.em ~        s(BASAL.em, k = kk)  + s(ndep, k = kk) + s(BASAL.plot, k = kk) + s(stem.density, k = kk) +
                  s(PC1, k=kk) + s(PC2, k=kk) + s(PC3, k=kk) + s(PC4, k=kk) + s(PC5, k=kk) + s(PC6, k=kk) + s(PC7, k=kk) + s(PC8, k=kk) + s(PC9, k=kk) + s(PC10, k =kk), 
                data = d1, family = 'poisson', method = 'REML')
R.mod.am <- gam(recruit.am ~        s(BASAL.am, k = kk)  + s(ndep, k = kk) + s(BASAL.plot, k = kk) + s(stem.density, k = kk) +
                  s(PC1, k=kk) + s(PC2, k=kk) + s(PC3, k=kk) + s(PC4, k=kk) + s(PC5, k=kk) + s(PC6, k=kk) + s(PC7, k=kk) + s(PC8, k=kk) + s(PC9, k=kk) + s(PC10, k =kk), 
                data = d1, family = 'poisson', method = 'REML')
M.mod.em <- gam(mortality ~          s(ndep, k = kk) + s(BASAL.plot, k = kk) + s(stem.density, k = kk) + s(PREVDIA.cm, k=kk) +
                  s(PC1, k=kk) + s(PC2, k=kk) + s(PC3, k=kk) + s(PC4, k=kk) + s(PC5, k=kk) + s(PC6, k=kk) + s(PC7, k=kk) + s(PC8, k=kk) + s(PC9, k=kk) + s(PC10, k =kk), 
                data = d2[d2$em == 1,], family = 'binomial', method = 'REML')
M.mod.am <- gam(mortality ~          s(ndep, k = kk) + s(BASAL.plot, k = kk) + s(stem.density, k = kk) + s(PREVDIA.cm, k=kk) +
                  s(PC1, k=kk) + s(PC2, k=kk) + s(PC3, k=kk) + s(PC4, k=kk) + s(PC5, k=kk) + s(PC6, k=kk) + s(PC7, k=kk) + s(PC8, k=kk) + s(PC9, k=kk) + s(PC10, k =kk), 
                data = d2[d2$em == 0,], family = 'binomial', method = 'REML')
G.mod.em <- gam(DIA.cm   ~          s(ndep, k = kk) + s(BASAL.plot, k = kk) + s(stem.density, k = kk) + s(PREVDIA.cm, k=kk) +
                  s(PC1, k=kk) + s(PC2, k=kk) + s(PC3, k=kk) + s(PC4, k=kk) + s(PC5, k=kk) + s(PC6, k=kk) + s(PC7, k=kk) + s(PC8, k=kk) + s(PC9, k=kk) + s(PC10, k =kk), 
                data = d2[d2$em == 1,], method = 'REML')
G.mod.am <- gam(DIA.cm   ~          s(ndep, k = kk) + s(BASAL.plot, k = kk) + s(stem.density, k = kk) + s(PREVDIA.cm, k=kk) +
                  s(PC1, k=kk) + s(PC2, k=kk) + s(PC3, k=kk) + s(PC4, k=kk) + s(PC5, k=kk) + s(PC6, k=kk) + s(PC7, k=kk) + s(PC8, k=kk) + s(PC9, k=kk) + s(PC10, k =kk), 
                data = d2[d2$em == 0,], method = 'REML')
cat('Null models fit.\n');toc()
#wrap output and name.
n.feedback <- list(G.mod.am, G.mod.em, M.mod.am, M.mod.em, R.mod.am, R.mod.em)
names(n.feedback) <- c('G.mod.am','G.mod.em','M.mod.am','M.mod.em','R.mod.am','R.mod.em')

#Environmental models with feedbacks.
cat('Fitting feedback models...\n');tic()
  R.mod.em <- gam(recruit.em ~  s(relEM, k = kk) + s(BASAL.em, k = kk) + s(ndep, k = kk) + s(BASAL.plot, k = kk) + s(stem.density, k = kk) +
                    s(PC1, k=kk) + s(PC2, k=kk) + s(PC3, k=kk) + s(PC4, k=kk) + s(PC5, k=kk) + s(PC6, k=kk) + s(PC7, k=kk) + s(PC8, k=kk) + s(PC9, k=kk) + s(PC10, k =kk), 
                  data = d1, family = 'poisson', method = 'REML')
  R.mod.am <- gam(recruit.am ~  s(relEM, k = kk) + s(BASAL.em, k = kk)  + s(ndep, k = kk) + s(BASAL.plot, k = kk) + s(stem.density, k = kk) +
                  s(PC1, k=kk) + s(PC2, k=kk) + s(PC3, k=kk) + s(PC4, k=kk) + s(PC5, k=kk) + s(PC6, k=kk) + s(PC7, k=kk) + s(PC8, k=kk) + s(PC9, k=kk) + s(PC10, k =kk), 
                  data = d1, family = 'poisson', method = 'REML')
  M.mod.em <- gam(mortality  ~  s(relEM, k = kk) + s(ndep, k = kk) + s(BASAL.plot, k = kk) + s(stem.density, k = kk) + s(PREVDIA.cm, k=kk) +
                  s(PC1, k=kk) + s(PC2, k=kk) + s(PC3, k=kk) + s(PC4, k=kk) + s(PC5, k=kk) + s(PC6, k=kk) + s(PC7, k=kk) + s(PC8, k=kk) + s(PC9, k=kk) + s(PC10, k =kk), 
                  data = d2[d2$em == 1,], family = 'binomial', method = 'REML')
  M.mod.am <- gam(mortality  ~  s(relEM, k = kk) + s(ndep, k = kk) + s(BASAL.plot, k = kk) + s(stem.density, k = kk) + s(PREVDIA.cm, k=kk) +
                  s(PC1, k=kk) + s(PC2, k=kk) + s(PC3, k=kk) + s(PC4, k=kk) + s(PC5, k=kk) + s(PC6, k=kk) + s(PC7, k=kk) + s(PC8, k=kk) + s(PC9, k=kk) + s(PC10, k =kk), 
                  data = d2[d2$em == 0,], family = 'binomial', method = 'REML')
  G.mod.em <- gam(DIA.cm   ~  s(relEM, k = kk) +  s(ndep, k = kk) + s(BASAL.plot, k = kk) + s(stem.density, k = kk) + s(PREVDIA.cm, k=kk) +
                  s(PC1, k=kk) + s(PC2, k=kk) + s(PC3, k=kk) + s(PC4, k=kk) + s(PC5, k=kk) + s(PC6, k=kk) + s(PC7, k=kk) + s(PC8, k=kk) + s(PC9, k=kk) + s(PC10, k =kk), 
                  data = d2[d2$em == 1,], method = 'REML')
  G.mod.am <- gam(DIA.cm   ~  s(relEM, k = kk) + s(ndep, k = kk) + s(BASAL.plot, k = kk) + s(stem.density, k = kk) + s(PREVDIA.cm, k=kk) +
                  s(PC1, k=kk) + s(PC2, k=kk) + s(PC3, k=kk) + s(PC4, k=kk) + s(PC5, k=kk) + s(PC6, k=kk) + s(PC7, k=kk) + s(PC8, k=kk) + s(PC9, k=kk) + s(PC10, k =kk), 
                  data = d2[d2$em == 0,], method = 'REML')
  cat('Feedback models fit.\n');toc()


#wrap output and name.
y.feedback <- list(G.mod.am, G.mod.em, M.mod.am, M.mod.em, R.mod.am, R.mod.em)
names(y.feedback) <- c('G.mod.am','G.mod.em','M.mod.am','M.mod.em','R.mod.am','R.mod.em')

#Grab plot environmental covariates for reference.----
#all plot-level environmental covariates. we will sample this later when simulating forests.
all.cov <- d1[,c('ndep','mat','map','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')]
#mean plot-level environmental covariates.
    cov <- colMeans(all.cov)
    names(cov) <- colnames(all.cov)

#Save models and size categories.----
output <- list(n.feedback, y.feedback, cov, all.cov)
names(output) <- c('n.feedback', 'y.feedback','env.cov','all.cov')
saveRDS(output, output.path, version = 2)

#plot check.
new.dat <- data.frame(22,14151,27)
colnames(new.dat) <- c('PREVDIA.cm','BASAL.plot','BASAL.em','BASAL.am','stem.density')
new.dat <- cbind(new.dat, data.frame(t(cov)))
new.dat <- data.frame(seq(0,1,by =0.01), new.dat)
colnames(new.dat)[1] <- 'relEM'

pred.r.am <- predict(R.mod.am, newdata = new.dat)
pred.r.em <- predict(R.mod.em, newdata = new.dat)
pred.m.am <- predict(M.mod.am, newdata = new.dat)
pred.m.em <- predict(M.mod.em, newdata = new.dat)
pred.g.am <- predict(G.mod.am, newdata = new.dat)
pred.g.em <- predict(G.mod.em, newdata = new.dat)

par(mfrow = c(3,2))
plot(pred.r.am ~ new.dat$relEM, bty = 'l', ylim = c(min(c(pred.r.am,pred.r.em)), max(c(pred.r.am,pred.r.em))))
plot(pred.r.em ~ new.dat$relEM, bty = 'l', ylim = c(min(c(pred.r.am,pred.r.em)), max(c(pred.r.am,pred.r.em))))
plot(pred.m.am ~ new.dat$relEM, bty = 'l', ylim = c(min(c(pred.m.am,pred.m.em)), max(c(pred.m.am,pred.m.em))))
plot(pred.m.em ~ new.dat$relEM, bty = 'l', ylim = c(min(c(pred.m.am,pred.m.em)), max(c(pred.m.am,pred.m.em))))
plot(pred.g.am ~ new.dat$relEM, bty = 'l', ylim = c(min(c(pred.g.am,pred.g.em)), max(c(pred.g.am,pred.g.em))))
plot(pred.g.em ~ new.dat$relEM, bty = 'l', ylim = c(min(c(pred.g.am,pred.g.em)), max(c(pred.g.am,pred.g.em))))


#Check variogram.----
vario.check <- F
if(vario.check == T){
  library(sp)
  library(gstat)
  dat <- data.frame(d1$recruit.em,d1$longitude.y,d1$latitude.y)
  colnames(dat) <- c('var','y','x')
  coordinates(dat) =  ~ x + y
  sp.mod <- variogram(var~1, data=dat)
  plot(sp.mod)
  
}
test <- predict.bam(a, newdata = new.dat, exclude = 's(PLT_CN)', newdata.guaranteed = T)


a <- gam(recruit.am ~  s(relEM, k = kk) + s(ndep, k = kk) + s(BASAL.plot, k = kk) + s(stem.density, k = kk)+
                  s(PC1, k=kk) + s(PC2, k=kk) + s(PC3, k=kk) + s(PC4, k=kk) + s(PC5, k=kk) + s(PC6, k=kk) + s(PC7, k=kk) + s(PC8, k=kk) + s(PC9, k=kk) + s(PC10, k =kk), 
                data = d1, family = 'poisson', method = 'REML')

b <- gam(recruit.am ~  s(relEM, k = kk) + s(BASAL.AM, k = kk) + s(ndep, k = kk) + s(BASAL.plot, k = kk) + s(stem.density, k = kk)+
           s(PC1, k=kk) + s(PC2, k=kk) + s(PC3, k=kk) + s(PC4, k=kk) + s(PC5, k=kk) + s(PC6, k=kk) + s(PC7, k=kk) + s(PC8, k=kk) + s(PC9, k=kk) + s(PC10, k =kk), 
         data = d1, family = 'poisson', method = 'REML')

#get predictions.
all.cov <- d1[,c('ndep','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')]
#mean plot-level environmental covariates.
cov <- colMeans(all.cov)
names(cov) <- colnames(all.cov)
#plot check.
new.dat <- data.frame(22,14151,27,8210)
colnames(new.dat) <- c('PREVDIA.cm','BASAL.plot','stem.density','BASAL.AM')
new.dat <- cbind(new.dat, data.frame(t(cov)))
new.dat <- data.frame(seq(0,1,by =0.01), new.dat)
colnames(new.dat)[1] <- 'relEM'


pred.a <- predict(a, newdata = new.dat)
pred.b <- predict(b, newdata = new.dat)
plot(pred.a ~ new.dat$relEM, bty = 'l', ylim = c(min(c(pred.a,pred.b)), max(c(pred.a,pred.b))), cex = 0)
lines(smooth.spline(pred.a ~ new.dat$relEM), lwd = 2)
lines(smooth.spline(pred.b ~ new.dat$relEM), lwd = 2, col = 'purple')
