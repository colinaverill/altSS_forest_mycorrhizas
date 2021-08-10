#Visualized biomdality - feedback strength analysis.
rm(list=ls())
source('paths.r')
library(boot)

#set output path.----
output.path <- 'figures/Fig._S3._bimod_feedback_strength.jpeg'

#load data.----
all.dat <- readRDS(bimod_feedback_str_analysis.path)

#jpeg save line.----
jpeg(output.path, width=8, height=5, units='in', res=300)

#plot results.----
par(mfrow = c(1,2))
#mortality only fit.
plot(logit(mono_frac) ~ mort.str, data = all.dat, 
     bty = 'l', pch = 16, cex=1.5,
     xlab = 'Survival Advantage',
     ylab = NA)
fit.mort <- lm(logit(mono_frac) ~ mort.str, data = all.dat)
rsq     <- round(summary(fit.mort)$adj.r.squared, 2)
rsq.raw <- round(summary(fit.mort)$r.squared    , 2)
rsq.lab <- bquote({R^2}[adj] == .(rsq))
rsq.raw.lab <- bquote({R^2} == .(rsq.raw))
abline(fit.mort, lwd = 2)
mtext(rsq.raw.lab, side=3, line=-2.0, adj=1, at = 0.5)
mtext(rsq.lab    , side=3, line=-3.5, adj=1, at = 0.5)
mtext('P < 0.001', side=3, line=-4.5, adj=1, at = 0.5)
mtext('Mycorrhizal Monodominance \n Metric',side = 2, line = 2)
mtext('A', side = 1, line = -1.5, adj = 0.95)

#recruitment only fit.
plot(logit(mono_frac) ~ recr.str, data = all.dat, 
     bty = 'l', pch = 16, cex=1.5,
     xlab = 'Recruitment Advantage',
     ylab = NA)
fit.recr <- lm(logit(mono_frac) ~ recr.str, data = all.dat)
#abline(fit.recr, lwd = 1, lty=2, col = 'light grey')
mtext('Mycorrhizal Monodominance \n Metric',side = 2, line = 2)
mtext('N.S.', side=3, line=-2, adj=0.95)
mtext('B', side = 1, line = -1.5, adj = 0.95)

#end plot.----
dev.off()