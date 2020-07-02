rm(list = ls())
library(gstat)
library(sp)
source('paths.r')
d <- readRDS(demographic_fits_gam_separate_plus_re_county.path)
d <- readRDS(demographic_fits_gam_separate.path)
d1 <- (readRDS(Product_1.path))
d2 <- (readRDS(Product_2.subset.path))
d1 <- d1[d1$PLT_CN %in% d2$PLT_CN,]

#grab eresiduals, as well as lat/lon
resid <- residuals(d$y.feedback$R.mod.em)
lat <- d1$latitude.x
lon <- d1$longitude.x
test.dat <- data.frame(resid, lat, lon)

#get semivariances and
coordinates(test.dat) <- ~ lon + lat
vario.1 <- variogram(resid~ 1, data = test.dat) #get semivariance values
vario.model.1 <- fit.variogram(vario.1, vgm(model='Sph', nugget=T)) #fit a variogram model with a Spherical error structure and a "nugget" (nugget = intercept)
plot(vario.1, model=vario.model.1, xlab = 'distance (km)')
