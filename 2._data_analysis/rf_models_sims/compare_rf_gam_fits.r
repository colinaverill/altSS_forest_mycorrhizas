rm(list=ls())
source('paths.r')
library(mgcv)
library(randomForest)

d.gam <- readRDS(demographic_fits_gam_separate.path)
d.for <- readRDS(rf_demographic_fits.path)

#get predictions across range.
new.dat <- data.frame(22,14151,27)
colnames(new.dat) <- c('PREVDIA.cm','BASAL.plot','stem.density')
new.dat <- cbind(new.dat, data.frame(t(d.gam$env.cov)))
new.dat <- data.frame(seq(0.01,0.99,by =0.01), new.dat)
colnames(new.dat)[1] <- 'relEM'


#mortality
m.am.gam <- predict(d.gam$y.feedback$M.mod.am, newdata = new.dat)
m.em.gam <- predict(d.gam$y.feedback$M.mod.em, newdata = new.dat)
m.am.rf  <- predict(d.for$models$mort.mod.am , newdata = new.dat)
m.em.rf  <- predict(d.for$models$mort.mod.em , newdata = new.dat)
m.am.rf <- boot::logit(m.am.rf)
m.em.rf <- boot::logit(m.em.rf)
limy.m.am <- c(min(c(m.am.gam,m.am.rf)), max(c(m.am.gam,m.am.rf)))
limy.m.em <- c(min(c(m.em.gam,m.em.rf)), max(c(m.em.gam,m.em.rf)))
r.am.gam <- predict(d.gam$y.feedback$R.mod.am, newdata = new.dat)
r.em.gam <- predict(d.gam$y.feedback$R.mod.em, newdata = new.dat)
r.am.rf  <- predict(d.for$models$recr.mod.am , newdata = new.dat)
r.em.rf  <- predict(d.for$models$recr.mod.em , newdata = new.dat)
r.am.rf <- log(r.am.rf)
r.em.rf <- log(r.em.rf)
limy.r.am <- c(min(c(r.am.gam,r.am.rf)), max(c(r.am.gam,r.am.rf)))
limy.r.em <- c(min(c(r.em.gam,r.em.rf)), max(c(r.em.gam,r.em.rf)))


par(mfrow = c(2,2))
#AM mortality
plot(m.am.gam ~ new.dat$relEM, bty = 'l', cex = 0, ylim = limy.m.am)
lines(smooth.spline(m.am.gam ~ new.dat$relEM), lwd = 2, col = 'purple')
lines(smooth.spline(m.am.rf  ~ new.dat$relEM), lwd = 2, col = 'green')
#EM mortality
plot(m.em.gam ~ new.dat$relEM, bty = 'l', cex = 0, ylim = limy.m.em)
lines(smooth.spline(m.em.gam ~ new.dat$relEM), lwd = 2, col = 'purple')
lines(smooth.spline(m.em.rf  ~ new.dat$relEM), lwd = 2, col = 'green')
#AM recruitment
plot(r.am.gam ~ new.dat$relEM, bty = 'l', cex = 0, ylim = limy.r.am)
lines(smooth.spline(r.am.gam ~ new.dat$relEM), lwd = 2, col = 'purple')
lines(smooth.spline(r.am.rf  ~ new.dat$relEM), lwd = 2, col = 'green')
#EM recruitment
plot(r.em.gam ~ new.dat$relEM, bty = 'l', cex = 0, ylim = limy.r.em)
lines(smooth.spline(r.em.gam ~ new.dat$relEM), lwd = 2, col = 'purple')
lines(smooth.spline(r.em.rf  ~ new.dat$relEM), lwd = 2, col = 'green')




