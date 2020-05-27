rm(list=ls())
source('paths.r')

#load ecoregion data.
d <- readRDS(demographic_fits_gam_ecoregion.path)
names(d) <- c('atlantic_highlands','mixed_wood_plains','mixed_wood_shield','appalachian_forests','southeastern_plains')
lab <- names(d)


#setup save paths.
out.dir <- 'figures/ecoregion_plots/'
cmd <- paste0('mkdir -p ',out.dir)
system(cmd)
output.paths <- paste0(out.dir,paste0(0,1:5),'_',lab,'.png')

for(i in 1:length(output.paths)){
  #grab ecoregion results of interest.----
         z <- d[[i]] 
  path.out <- output.paths[i]
  
  #setup covariates, generate predictions.----
  am.dat <- c(0,z$cov)
  names(am.dat)[1] <- 'relEM'
  em.dat <- am.dat
  em.dat[1] <- 1
  dat <- data.frame(rbind(am.dat, em.dat))
  #make predictions on linear scale.
  am.mort <- predict(z$M.mod.am, newdata = dat, se.fit = T)
  em.mort <- predict(z$M.mod.em, newdata = dat, se.fit = T)
  am.recr <- predict(z$R.mod.am, newdata = dat, se.fit = T)
  em.recr <- predict(z$R.mod.em, newdata = dat, se.fit = T)
  N <- 300
  y.am.am.M <- rnorm(N, am.mort$fit[1], am.mort$se.fit[1])
  y.am.em.M <- rnorm(N, am.mort$fit[2], am.mort$se.fit[2])
  y.em.am.M <- rnorm(N, em.mort$fit[1], em.mort$se.fit[1])
  y.em.em.M <- rnorm(N, em.mort$fit[2], em.mort$se.fit[2])
  y.am.am.R <- rnorm(N, am.recr$fit[1], am.recr$se.fit[1])
  y.am.em.R <- rnorm(N, am.recr$fit[2], am.recr$se.fit[2])
  y.em.am.R <- rnorm(N, em.recr$fit[1], em.recr$se.fit[1])
  y.em.em.R <- rnorm(N, em.recr$fit[2], em.recr$se.fit[2])
  
  #back transform mortality and recruitment data.
  y.am.M <- boot::inv.logit(c(y.am.am.M,y.am.em.M))
  y.em.M <- boot::inv.logit(c(y.em.am.M,y.em.em.M))
  y.am.R <-             exp(c(y.am.am.R,y.am.em.R))
  y.em.R <-             exp(c(y.em.am.R,y.em.em.R))
  
  #setup x positions.
  jit <- 0.05
  x.pos <- c(0.5,1,1.5,2)
  x1 <- rnorm(N,x.pos[1],jit)
  x2 <- rnorm(N,x.pos[2],jit)
  x3 <- rnorm(N,x.pos[3],jit)
  x4 <- rnorm(N,x.pos[4],jit)
  x.am <- c(x1,x2)
  x.em <- c(x3,x4)
  
  #Start plots.----
  png(path.out, width = 8, height = 3.5, units = 'in', res = 300)
  par(mfrow = c(1,2),
      mar = c(1,3.5,1,1))
  
  #Mortality plot.
  cols <- c('green','purple')
  col.1 <- c(rep(cols[1],N),rep(cols[2],N))
  limx <- c(min(c(x.am,x.em)), max(c(x.am,x.em)))
  limy <- c(0, max(c(y.am.M,y.em.M))*1.05)
  #plot.
  plot(y.am.M ~ x.am, pch = 16, xlim = limx, ylim = limy, cex = 0.3, bty = 'l', ylab = NA, xlab = NA, xaxt = 'n', col = col.1)
  points(y.em.M ~ x.em, pch = 16, cex = 0.3, col = col.1)
  abline(v = 1.25, lwd=2, lty = 3, col = 'light gray')
  #labels.
  mtext('Mortality Probability', side = 2, cex = 1.0, line = 2.5)
  legend('bottomleft',legend = c('ectomycorrhizal','arbuscular mycorrhizal'), 
         pch = 16, col = rev(cols), bty = 'n', cex = 0.9,
         x.intersp = .75, xpd = T, 
         horiz = F)
  mtext('AM Forest',side = 3, line = -1, adj = 0.125, cex = 1)
  mtext('EM Forest',side = 3, line = -1, adj = 0.9, cex = 1)
  
  #Recruitment plot.
  limy <- c(0, max(c(y.am.R,y.em.R))*1.05)
  #plot.
  plot(y.am.R ~ x.am, pch = 16, xlim = limx, ylim = limy, cex = 0.3, bty = 'l', ylab = NA, xlab = NA, xaxt = 'n', col = col.1)
  points(y.em.R ~ x.em, pch = 16, cex = 0.3, col = col.1)
  abline(v = 1.25, lwd=2, lty = 3, col = 'light gray')
  #labels.
  lab <- expression(paste('Recruitment - Poisson ',lambda))
  mtext(lab, side = 2, cex = 1.0, line = 2.5)
  mtext('AM Forest',side = 3, line = -1, adj = 0.125, cex = 1)
  mtext('EM Forest',side = 3, line = -1, adj = 0.9, cex = 1)
  
  #end plots.-----
  dev.off()
  
}
