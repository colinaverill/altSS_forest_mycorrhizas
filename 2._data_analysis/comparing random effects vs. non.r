#compare old vs. new feedbacks.
rm(list = ls())
library(mgcv)
source('paths.r')

old <- readRDS(demographic_fits_gam_separate.path)
new <- readRDS(demographic_fits_gam_separate_plus_re_county.path)

cov <- old$env.cov
new.dat <- data.frame(22,14151,27,6500,6500)
colnames(new.dat) <- c('PREVDIA.cm','BASAL.plot','stem.density','BASAL.am','BASAL.em')
new.dat <- cbind(new.dat, data.frame(t(cov)))
new.dat <- data.frame(seq(0,1,by =0.01), new.dat)
colnames(new.dat)[1] <- 'relEM'

pred.r.am.old <- predict(old$y.feedback$R.mod.am, newdata = new.dat)
pred.r.em.old <- predict(old$y.feedback$R.mod.em, newdata = new.dat)
pred.m.am.old <- predict(old$y.feedback$M.mod.am, newdata = new.dat, exclude = 's(county.ID)', newdata.guaranteed = T)
pred.m.em.old <- predict(old$y.feedback$M.mod.em, newdata = new.dat, exclude = 's(county.ID)', newdata.guaranteed = T)
pred.g.am.old <- predict(old$y.feedback$G.mod.am, newdata = new.dat, exclude = 's(county.ID)', newdata.guaranteed = T)
pred.g.em.old <- predict(old$y.feedback$G.mod.em, newdata = new.dat, exclude = 's(county.ID)', newdata.guaranteed = T)

pred.r.am.new <- predict(new$y.feedback$R.mod.am, newdata = new.dat)
pred.r.em.new <- predict(new$y.feedback$R.mod.em, newdata = new.dat)
pred.m.am.new <- predict(new$y.feedback$M.mod.am, newdata = new.dat, exclude = 's(county.ID)', newdata.guaranteed = T)
pred.m.em.new <- predict(new$y.feedback$M.mod.em, newdata = new.dat, exclude = 's(county.ID)', newdata.guaranteed = T)
pred.g.am.new <- predict(new$y.feedback$G.mod.am, newdata = new.dat, exclude = 's(county.ID)', newdata.guaranteed = T)
pred.g.em.new <- predict(new$y.feedback$G.mod.em, newdata = new.dat, exclude = 's(county.ID)', newdata.guaranteed = T)

par(mfrow = c(3,2))
#mortality.
plot(pred.m.am.old ~ new.dat$relEM, bty = 'l', cex = 0,
     ylim = c(min(c(pred.m.am.old,pred.m.am.new,pred.m.em.old,pred.m.em.new)), max(c(pred.m.am.old,pred.m.em.old,pred.m.am.new,pred.m.em.new))))
lines(smooth.spline(pred.m.am.old ~ new.dat$relEM), lwd = 2)
lines(smooth.spline(pred.m.am.new ~ new.dat$relEM), lwd = 2, col = 'purple')
plot(pred.m.em.old ~ new.dat$relEM, bty = 'l', cex = 0,
     ylim = c(min(c(pred.m.am.old,pred.m.am.new,pred.m.em.old,pred.m.em.new)), max(c(pred.m.am.old,pred.m.em.old,pred.m.am.new,pred.m.em.new))))
lines(smooth.spline(pred.m.em.old ~ new.dat$relEM), lwd = 2)
lines(smooth.spline(pred.m.em.new ~ new.dat$relEM), lwd = 2, col = 'purple')
#growth.
plot(pred.g.am.old ~ new.dat$relEM, bty = 'l', cex = 0,
     ylim = c(min(c(pred.g.am.old,pred.g.am.new,pred.g.em.old,pred.g.em.new)), max(c(pred.g.am.old,pred.g.em.old,pred.g.am.new,pred.g.em.new))))
lines(smooth.spline(pred.g.am.old ~ new.dat$relEM), lwd = 2)
lines(smooth.spline(pred.g.am.new ~ new.dat$relEM), lwd = 2, col = 'purple')
plot(pred.g.em.old ~ new.dat$relEM, bty = 'l', cex = 0,
     ylim = c(min(c(pred.g.am.old,pred.g.am.new,pred.g.em.old,pred.g.em.new)), max(c(pred.g.am.old,pred.g.em.old,pred.g.am.new,pred.g.em.new))))
lines(smooth.spline(pred.g.em.old ~ new.dat$relEM), lwd = 2)
lines(smooth.spline(pred.g.em.new ~ new.dat$relEM), lwd = 2, col = 'purple')
#recruitment. should be nearly identical. It is.
plot(pred.r.am.old ~ new.dat$relEM, bty = 'l', cex = 0,
     ylim = c(min(c(pred.r.am.old,pred.r.am.new,pred.r.em.old,pred.r.em.new)), max(c(pred.r.am.old,pred.r.em.old,pred.r.am.new,pred.r.em.new))))
lines(smooth.spline(pred.r.am.old ~ new.dat$relEM), lwd = 2)
lines(smooth.spline(pred.r.am.new ~ new.dat$relEM), lwd = 2, col = 'purple')
plot(pred.r.em.old ~ new.dat$relEM, bty = 'l', cex = 0,
     ylim = c(min(c(pred.r.am.old,pred.r.am.new,pred.r.em.old,pred.r.em.new)), max(c(pred.r.am.old,pred.r.em.old,pred.r.am.new,pred.r.em.new))))
lines(smooth.spline(pred.r.em.old ~ new.dat$relEM), lwd = 2)
lines(smooth.spline(pred.r.em.new ~ new.dat$relEM), lwd = 2, col = 'purple')


