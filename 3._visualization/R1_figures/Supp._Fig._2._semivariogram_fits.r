#plot variograms of raw and fitted residuals of AM and EM recruitment.
rm(list=ls())
source('paths.r')

#output path.----
output.path <- 'figures/Supplemental_Fig._X._semivariograms_fits.png'

#load data.----
d <- readRDS(variogram_data.path)
vario.list       <- d$vario.list
vario.model.list <- d$vario.model.list

#plot results.----
#png save line.
png(output.path, width = 8, height = 10, units = 'in', res = 300)

#plot
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

#end plot.----
dev.off()
