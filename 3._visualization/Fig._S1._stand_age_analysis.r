#plot stand age results.
rm(list=ls())
source('paths.r')
library(boot)

#set output path.----
#output.path <- 'figures/Fig._S1._stand_age_analysis.jpeg'
output.path <- 'figures/Extended_Data_Fig._1._stand_age_analysis.pdf'

#load data.----
d <- readRDS(stand_age_results.path)

#save line.----
#jpeg(output.path, width = 8, height = 5, units = 'in', res =300)
pdf(output.path, width = 8, height = 5)

#plot age relationship after detrending.----
par(mfrow = c(1,2))
#panel 1.
plot(inv.logit(d$detrend1) ~ d$STDAGE, cex = 0.4, pch = 16, bty = 'l',
     xlab = 'Stand Age', ylab = 'Relative Abundance EM Trees')
fit <- lm(inv.logit(d$detrend1) ~ d$STDAGE)
abline(fit, col = 'green', lwd = 2)
mtext('A',side = 3, line = 0, adj=0.05)

#panel 2.
hist(inv.logit(d$detrend2), breaks = 30, ylab = 'Number of Forests',xlab='Relative Abundance EM Trees', main = NA)
mtext('B',side = 3, line = 0, adj=0.05)

#end plot.----
dev.off()