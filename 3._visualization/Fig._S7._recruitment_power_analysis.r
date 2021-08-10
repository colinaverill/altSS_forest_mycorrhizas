#plot recruitment power analysis, and correlations between total abundance and relative abundance of AM/EM trees.
rm(list=ls())
library(data.table)
library(car)
source('paths.r')

#set output path.----
output.path <- 'figures/Fig._S7._recruitment_power_analysis.jpeg'

#load plot data and growth/mortality/recruitment data.----
d1 <- data.table(readRDS(Product_1.path))
d2 <- data.table(readRDS(Product_2.subset.path))
results <- readRDS(recruit_power_analysis.path)

#Subset plot data and rename some things.-----
d1 <- d1[d1$PLT_CN %in% d2$PLT_CN,]
d2$PLT_CN <- as.factor(d2$PLT_CN)
d2$county.ID <- as.factor(paste0(d2$STATECD,'_',d2$COUNTYCD))
d1$county.ID <- as.factor(paste0(d1$STATECD,'_',d1$COUNTYCD))
setnames(d1,'plot.BASAL','BASAL.plot')
setnames(d2,'plot.BASAL','BASAL.plot')
setnames(d1,'BASAL.ECM','BASAL.em')
setnames(d1,'BASAL.AM' ,'BASAL.am')

#save line.----
jpeg(output.path, width=8, height=4, units='in', res=300)

#global plot settings.----
par(mfrow=c(1,3))
#EM corrleation.----
plot(BASAL.em ~ relEM, data = d1, bty = 'l', cex = 0.3)
fit.em <- lm(BASAL.em ~ relEM, data= d1)
x <- data.frame(seq(0, 1, by = 0.01))
colnames(x) <- 'relEM'
pred.em <- predict(fit.em, newdata = x)
#plot spline fit.
lines(smooth.spline(pred.em ~ x[,1]), lwd = 2, col = 'green')
#report rsq and VIF.
rsq <- round(summary(fit.em)$r.sq, 2)
vif <- round(vif(lm(mat ~ relEM + BASAL.em, data = d1))[1], 2)
rsq.lab <- bquote({R^2} == .(rsq))
vif.lab <- paste0('VIF = ',vif)
mtext(rsq.lab, side = 3, line = -2.0, adj = 0.05)
mtext(vif.lab, side = 3, line = -3.4, adj = 0.05)

#AM plot.----
plot(BASAL.am ~ relEM, data = d1, bty = 'l', cex = 0.3,
     ylab = 'Basal Area AM Trees',
     xlab = 'relative abundance EM trees')
fit.am <- lm(BASAL.am ~ relEM, data = d1)
x <- data.frame(seq(0, 1, by = 0.01))
colnames(x) <- 'relEM'
pred.am <- predict(fit.am, newdata = x)
#plot spline fit.
lines(smooth.spline(pred.am ~ x[,1]), lwd = 2, col = 'green')
#report rsq and VIF.
rsq <- round(summary(fit.am)$r.sq, 2)
formatC(summary(fit.am)$r.sq, digits = 2, drop0trailing = F)
vif <- round(vif(lm(mat ~ relEM + BASAL.am, data = d1))[1], 2)
rsq.lab <- bquote({R^2} == .(rsq))
vif.lab <- paste0('VIF = ',vif)
mtext(rsq.lab, side = 3, line = -2.0, adj = 0.05)
mtext(vif.lab, side = 3, line = -3.4, adj = 0.05)


#panel 3. Power analysis.----
x <- c(1:nrow(results))
plot(results$mu ~ x, pch = 16, bty='l',cex =2,
     ylim = c(min(results$lo95)*1.05, max(results$hi95)*1.05),
     xaxt='n', ylab=NA, xlab=NA)
abline(h=2, lwd = 2, lty = 2)
arrows(x, results$lo95, x, results$hi95, length=0.04, angle=90, code=3)
mtext('parameter estimate',side = 2, line = 2.3)
mtext('sampling effort'  , side = 1, line = 3)
text(cex=1, x= c(1:nrow(results)), y=0.2, results$N, srt = 45, xpd=TRUE, adj=1)

#end plot.----
dev.off()
