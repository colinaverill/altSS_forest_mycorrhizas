#compare spatial autocorrelation in AM vs. EM recruitment among models.
rm(list=ls())
source('paths.r')
library(gstat)
library(sp)
library(mgcv)

#set output path.----
output.path <- variogram_data.path

#load data and models.----
d1 <- (readRDS(Product_1.path))
d2 <- (readRDS(Product_2.subset.path))
d1 <- d1[d1$PLT_CN %in% d2$PLT_CN,]
fit <- readRDS(demographic_fits_gam_separate_plus_re_county.path)

#Get AM and EM recruitment residuals of an intercept only model, and of the gams.----
raw.am.resid <- residuals(glm(d1$recruit.am ~ 1, family = poisson))
raw.em.resid <- residuals(glm(d1$recruit.em ~ 1, family = poisson))
fit.am.resid <- residuals(fit$y.feedback$R.mod.am)
fit.em.resid <- residuals(fit$y.feedback$R.mod.em)
resid.list <- list(raw.am.resid, raw.em.resid, fit.am.resid, fit.em.resid)



#fit variograms.----
#define spatial coordinates.
lat <- d1$latitude.x
lon <- d1$longitude.x

#fit you variograms.
vario.list <- list()
vario.model.list <- list()
for(i in 1:length(resid.list)){
  test.dat <- data.frame(resid.list[[i]], lat, lon)
  colnames(test.dat)[1] <- 'resid'
  coordinates(test.dat) <- ~ lon + lat
  vario       <-     variogram(resid~ 1, data = test.dat)           #get semivariance values
  vario.model <- fit.variogram(vario, vgm(model='Sph', nugget=T)) #fit a variogram model with a Spherical error structure and a "nugget" (nugget = intercept)
  vario.list[[i]] <- vario
  vario.model.list[[i]] <- vario.model
}

#Save output.----
lab <- c('raw.am','raw.em','fit.am','fit.em')
names(vario.list)       <- lab
names(vario.model.list) <- lab
output <- list(vario.list, vario.model.list)
names(output) <- c('vario.list','vario.model.list')
saveRDS(output, output.path)

#plot results.----
lab <- c('Raw AM Residuals','Raw EM Residuals','Fitted AM Residuals','Fitted EM Residuals')
par(mfrow = c(2,2))
for(i in 1:length(vario.list)){
  #plot semivariance as a function of distance.
  plot(vario.list[[i]]$gamma ~ vario.list[[i]]$dist, 
       bty = 'l', 
       ylim = c(0, max(vario.list[[i]]$gamma)*1.05), 
       pch = 16, 
       ylab = 'semivariance', xlab = 'distance (km)',
       main = lab[i])
  #calculate spherical fit line.
  par.mat <- vario.model.list[[i]]
  fit <- par.mat[1,2] + par.mat[2,2] * ((1.5 * (vario.list[[i]]$dist / par.mat[2,3])) - (0.5*(vario.list[[i]]$dist^3 / par.mat[2,3]^3)))
  fit <- ifelse(vario.list[[i]]$dist > par.mat[2,3], 
                par.mat[1,2] + par.mat[2,2],
                fit)
  #plot spherical fit line.
  lines(smooth.spline(fit ~ vario.list[[i]]$dist))
}

