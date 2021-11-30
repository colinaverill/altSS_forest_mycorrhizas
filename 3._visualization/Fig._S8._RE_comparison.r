rm(list=ls())
source('paths.r')
library(mgcv)

#set output path.----
output.path <- 'figures/Extended_Data_Fig._8._RE_comparison.pdf'

#load ecoregion data.
d <- readRDS(atlantic_highland_RE_test.path)

#save line.----
#jpeg(filename=output.path, width=8, height=6, units='in', res=300)
pdf(file=output.path, width=8, height=6)

#plot two plots separately and then on top of each other.----
par(mfrow = c(2,3))
#COUNTY VS PLOT RANDOM EFFECT.
plot(d$county.re$M.mod.em, select = 1, ylim = c(-0.5, 1.5), 
     bty = 'l', main = 'county \nrandom effect', col = 'red',
     xlab = 'relative abundance  \nEM trees', ylab = 'mortality index')
mtext('A',side=3, adj=0.95, line=-2)
plot(d$plot.re$M.mod.em  , select = 1, ylim = c(-0.5, 1.5), 
     bty = 'l', main = 'plot \nrandom effect', col = 'blue',
     xlab = 'relative abundance \nEM trees', ylab = 'mortality index')
mtext('B',side=3, adj=0.95, line=-2)
#together
plot(d$county.re$M.mod.em, select = 1, ylim = c(-0.5, 1.5), 
     bty = 'l', main = 'combined', col = 'blue',
     xlab = 'relative abundance \nEM trees', ylab = 'mortality index')
par(new = T)
plot(d$plot.re$M.mod.em  , select = 1, ylim = c(-0.5, 1.5), 
     bty = 'l', col = 'red', xlab = NA, ylab = NA, axes = F)
mtext('C',side=3, adj=0.95, line=-2)

#WITHOUT / WITH SPECIES RANDOM EFFECT.
plot(d$county.re$M.mod.em, select = 1, ylim = c(-0.5, 1.5), 
     bty = 'l', main = 'without species \nrandom effect', col = 'red',
     xlab = 'relative abundance  \nEM trees', ylab = 'mortality index')
mtext('D',side=3, adj=0.95, line=-2)
plot(d$county.spp.re$M.mod.em  , select = 1, ylim = c(-0.5, 1.5), 
     bty = 'l', main = 'with species \nrandom effect', col = 'blue',
     xlab = 'relative abundance \nEM trees', ylab = 'mortality index')
mtext('E',side=3, adj=0.95, line=-2)
#together
plot(d$county.re$M.mod.em, select = 1, ylim = c(-0.5, 1.5), 
     bty = 'l', main = 'combined', col = 'blue',
     xlab = 'relative abundance \nEM trees', ylab = 'mortality index')
par(new = T)
plot(d$county.spp.re$M.mod.em  , select = 1, ylim = c(-0.5, 1.5), 
     bty = 'l', col = 'red', xlab = NA, ylab = NA, axes = F)
mtext('F',side=3, adj=0.95, line=-2)

#end plot.----
dev.off()
