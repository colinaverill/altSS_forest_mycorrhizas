#Showing influenbce of each PC on the 6 models.
#Use county RE models.
rm(list=ls())
source('paths.r')
library(mgcv)

#load data.----
d <- readRDS(demographic_fits_gam_separate_plus_re_county.path)

#save line.----
jpeg(filename='figures/Supp.Fig._X._AM_PC_GRM_relationships.jpeg',width=14,height=7,units='in',res=300)

#global plot settings.----
par(mfrow = c(3, 10))
par(mar = c(4,3,2,1))
par(oma = c(1,4,1,1))

#Growth plots.----
plot(d$y.feedback$G.mod.am, select =  7, ylim = c(-1, 1), bty = 'l', ylab = NA)
mtext('growth rate', side = 2, line = 4)
lab <- expression(paste('cm * (5 years)'^-1))
mtext(lab, side = 2, cex = 0.7, line = 2.5)
plot(d$y.feedback$G.mod.am, select =  8, ylim = c(-1, 1), bty = 'l', ylab = NA)
plot(d$y.feedback$G.mod.am, select =  9, ylim = c(-1, 1), bty = 'l', ylab = NA)
plot(d$y.feedback$G.mod.am, select = 10, ylim = c(-1, 1), bty = 'l', ylab = NA)
plot(d$y.feedback$G.mod.am, select = 11, ylim = c(-1, 1), bty = 'l', ylab = NA)
plot(d$y.feedback$G.mod.am, select = 12, ylim = c(-1, 1), bty = 'l', ylab = NA)
plot(d$y.feedback$G.mod.am, select = 13, ylim = c(-1, 1), bty = 'l', ylab = NA)
plot(d$y.feedback$G.mod.am, select = 14, ylim = c(-1, 1), bty = 'l', ylab = NA)
plot(d$y.feedback$G.mod.am, select = 15, ylim = c(-1, 1), bty = 'l', ylab = NA)
plot(d$y.feedback$G.mod.am, select = 16, ylim = c(-1, 1), bty = 'l', ylab = NA)

#Recruitment plots.----
plot(d$y.feedback$R.mod.am, select =  7, ylim = c(-2, 2), bty = 'l', ylab = 'recruitment')
mtext('recruitment', side = 2, line = 4)
lab <- expression(paste('Poisson ',lambda))
mtext(lab, side = 2, cex = 0.7, line = 2.5)
plot(d$y.feedback$R.mod.am, select =  8, ylim = c(-2, 2), bty = 'l', ylab = NA)
plot(d$y.feedback$R.mod.am, select =  9, ylim = c(-2, 2), bty = 'l', ylab = NA)
plot(d$y.feedback$R.mod.am, select = 10, ylim = c(-2, 2), bty = 'l', ylab = NA)
plot(d$y.feedback$R.mod.am, select = 11, ylim = c(-2, 2), bty = 'l', ylab = NA)
plot(d$y.feedback$R.mod.am, select = 12, ylim = c(-2, 2), bty = 'l', ylab = NA)
plot(d$y.feedback$R.mod.am, select = 13, ylim = c(-2, 2), bty = 'l', ylab = NA)
plot(d$y.feedback$R.mod.am, select = 14, ylim = c(-2, 2), bty = 'l', ylab = NA)
plot(d$y.feedback$R.mod.am, select = 15, ylim = c(-2, 2), bty = 'l', ylab = NA)
plot(d$y.feedback$R.mod.am, select = 16, ylim = c(-2, 2), bty = 'l', ylab = NA)

#Mortality plots.----
plot(d$y.feedback$M.mod.am, select =  7, ylim = c(-2, 2), bty = 'l', ylab = 'mortality')
mtext('mortality', side = 2, line = 4)
mtext('logit(probability)', side = 2, cex = 0.7, line = 2.5)
plot(d$y.feedback$M.mod.am, select =  8, ylim = c(-2, 2), bty = 'l', ylab = NA)
plot(d$y.feedback$M.mod.am, select =  9, ylim = c(-2, 2), bty = 'l', ylab = NA)
plot(d$y.feedback$M.mod.am, select = 10, ylim = c(-2, 2), bty = 'l', ylab = NA)
plot(d$y.feedback$M.mod.am, select = 11, ylim = c(-2, 2), bty = 'l', ylab = NA)
plot(d$y.feedback$M.mod.am, select = 12, ylim = c(-2, 2), bty = 'l', ylab = NA)
plot(d$y.feedback$M.mod.am, select = 13, ylim = c(-2, 2), bty = 'l', ylab = NA)
plot(d$y.feedback$M.mod.am, select = 14, ylim = c(-2, 2), bty = 'l', ylab = NA)
plot(d$y.feedback$M.mod.am, select = 15, ylim = c(-2, 2), bty = 'l', ylab = NA)
plot(d$y.feedback$M.mod.am, select = 16, ylim = c(-2, 2), bty = 'l', ylab = NA)

#outer labels.----
mtext('Arbuscular mycorrhizal fits', side = 3, line = -1, outer = T, adj = 0.03)

#end plot.----
dev.off()

