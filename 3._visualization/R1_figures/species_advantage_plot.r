#Plot species level effects. Differences always performance under EM relative to performance under AM.
rm(list = ls())
source('paths.r')
library(mgcv)

#load data.----
d <- readRDS(demographic_fits_gam_species.path)
spp.key <- d$spp.key
d$spp.key <- NULL

#aggregate species specific contrasts under EM vs. AM condition.----
spp.result <- list()
for(i in 1:length(d)){
  spp.result[[i]] <- unlist(d[[i]]$post.contrast)
}
spp.result <- data.frame(do.call(rbind, spp.result))
spp.result$name <- names(d)
spp.key <- merge(spp.key, spp.result)

#setup pch key.----
spp.key$pch <- 16
spp.key$pch <- ifelse(spp.key$MYCO_ASSO == 'ECM' & spp.key$gymno == 1, 17, spp.key$pch)
spp.key$pch <- ifelse(spp.key$MYCO_ASSO ==  'AM'                     ,  1, spp.key$pch)
spp.key$pch <- ifelse(spp.key$MYCO_ASSO ==  'AM' & spp.key$gymno == 1,  2, spp.key$pch)

#Plot panels and save line.----
png('test.png', width = 6, height = 8, units = 'in', res = 300)
par(mfrow = c(2,1),mar = c(0.1,6,6,0.1), oma = c(1,1,0,1))

#recruitment plot.----
spp.key <- spp.key[order(spp.key$r.mu),]
x    <- spp.key$r.mu
x.sd <- spp.key$r.sd
y <- seq(1:nrow(spp.key))
limx <- max(abs(c((x - 1.06*x.sd),(x+1.06*x.sd))))
limx <- c(-limx, limx)

plot(y ~ x, bty = 'n', xaxt = 'n', yaxt = 'n', cex = NA, ylab = NA, xlab = NA, xlim = limx)
axis(3)
abline (v = 0, lty = 2, col = 'light gray')
#arrows(x - x.sd*1.96, y, x + x.sd*1.96, y, code = 3, angle=90, length=0.0)
arrows(x - x.sd*1.00, y, x + x.sd*1.00, y, code = 3, angle=90, length=0.0)
points(y ~ x, pch = spp.key$pch, cex = 1.2)
#y-axis species labels.
text(spp.key$name,cex = 0.8, y = y, x = min(limx), adj = 1, xpd = T)
mtext('Recruitment Advantage\namong AM Trees' , side = 3, line = 2, at = limx[1]*0.5, adj = 0.5, cex = 0.8)
mtext('Recruitment Advantage\namong ECM Trees', side = 3, line = 2, at = limx[2]*0.5, adj = 0.5, cex = 0.8)
#legend
legend(x = 1.5, y = 8, legend = c('AM hardwood','AM conifer', 'ECM hardwood','ECM confier'), pch = c(1,2,16,17) , ncol=, bg = NA, bty = 'n', cex = 0.8)

#mortality plot.----
spp.key <- spp.key[order(spp.key$m.mu),]
x    <- spp.key$m.mu
x.sd <- spp.key$m.sd
y <- seq(1:nrow(spp.key))
limx <- max(abs(c((x - 1.06*x.sd),(x+1.06*x.sd))))
limx <- c(-limx, limx)

plot(y ~ x, bty = 'n', xaxt = 'n', yaxt = 'n', cex = NA, ylab = NA, xlab = NA, xlim = limx)
axis(3)
abline (v = 0, lty = 2, col = 'light gray')
#arrows(x - x.sd*1.96, y, x + x.sd*1.96, y, code = 3, angle=90, length=0.0)
arrows(x - x.sd*1.00, y, x + x.sd*1.00, y, code = 3, angle=90, length=0.0)
points(y ~ x, pch = spp.key$pch, cex = 1.2)
#y-axis species labels.
text(spp.key$name,cex = 0.8, y = y, x = min(limx), adj = 1, xpd = T)
mtext('Survival Advantage\namong AM Trees' , side = 3, line = 2, at = limx[1]*0.5, adj = 0.5, cex = 0.8)
mtext('Survival Advantage\namong ECM Trees', side = 3, line = 2, at = limx[2]*0.5, adj = 0.5, cex = 0.8)


#end plot.----
dev.off()

