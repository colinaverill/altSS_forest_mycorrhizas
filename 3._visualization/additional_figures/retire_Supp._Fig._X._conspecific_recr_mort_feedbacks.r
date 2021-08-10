#Plot conspecific recruitment and survival effects at thhe species level, holding relEM constant.
#Convert mortality to survival.
rm(list=ls())
library(mgcv)
library(data.table)
source('paths.r')
source('project_functions/predict_gam_well.r')
#function to draw randomly from multivariate normal distribution.-----
rmvn <- function(n,mu,sig) { ## MVN random deviates
  L <- mroot(sig);m <- ncol(L);
  t(mu + L%*%matrix(rnorm(m*n),m,n)) 
}
#define number of posterior draws for plotting.----
N <- 300
#load species level recruitment and mortality data.----
d <- readRDS(demographic_fits_gam_species.path)
spp.key <- d$spp.key
d$spp.key <- NULL
g.mu <- rep(0, length(d))
g.sd <- rep(0, length(d))
r.mu <- rep(0, length(d))
r.sd <- rep(0, length(d))
m.mu <- rep(0, length(d))
m.sd <- rep(0, length(d))

#Grab responses.
out <- list()
for(i in 1:length(d)){
  recr.dat  <- plot(d[[i]]$R.mod, select = 2)
  recr.dat  <- recr.dat[[2]]$fit
  recr.diff <- recr.dat[length(recr.dat)] - recr.dat[1]
  mort.dat  <- plot(d[[i]]$M.mod, select = 2)
  mort.dat  <- mort.dat[[2]]$fit
  mort.diff <- mort.dat[length(mort.dat)] - mort.dat[1]
  out[[i]] <- c(recr.diff, mort.diff)
}
out <- data.frame(do.call(rbind, out))
colnames(out) <- c('recr.diff','mort.diff')
out$species <- names(d)
out$MYCO_ASSO <- spp.key$MYCO_ASSO
out$gymno <- spp.key$gymno
spp.key <- out

#setup pch key.
spp.key$pch <- 16
spp.key$pch <- ifelse(spp.key$MYCO_ASSO == 'ECM' & spp.key$gymno == 1, 17, spp.key$pch)
spp.key$pch <- ifelse(spp.key$MYCO_ASSO ==  'AM'                     ,  1, spp.key$pch)
spp.key$pch <- ifelse(spp.key$MYCO_ASSO ==  'AM' & spp.key$gymno == 1,  2, spp.key$pch)

par(mfrow = c(1,2))

#recruitment effects conspecific.----
spp.key <- spp.key[order(spp.key$recr.diff),]
x    <- spp.key$recr.diff
y <- seq(1:nrow(spp.key))
limx <- c(min(x)*1.05, max(x)*1.05)

plot(y ~ x, bty = 'n', xaxt = 'n', yaxt = 'n', cex = NA, ylab = NA, xlab = NA, xlim = limx)
axis(3)
abline (v = 0, lty = 2, col = 'light gray')
#arrows(x - x.sd*1.96, y, x + x.sd*1.96, y, code = 3, angle=90, length=0.0)
points(y ~ x, pch = spp.key$pch, cex = 1.2)
#y-axis species labels.
text(spp.key$species,cex = 0.8, y = y, x = min(limx), adj = 1, xpd = T)
mtext('Recruitment Advantage\namong Heterspecifics' , side = 3, line = 2, at = limx[1]*0.5, adj = 0.5, cex = 0.7)
mtext('Recruitment Advantage\namong Conspecifics', side = 3, line = 2, at = limx[2]*0.5, adj = 0.5, cex = 0.7)
mtext('a.', side = 1, line = -1.65, adj = 0.98, cex = 0.9, font = 1)
#legend
legend(x = 1.5, y = 8, legend = c('AM hardwood','AM conifer', 'EM hardwood','EM confier'), pch = c(1,2,16,17) , ncol=, bg = NA, bty = 'n', cex = 0.8)

#mortality effects conspecific.----
spp.key <- spp.key[order(spp.key$mort.diff),]
x    <- spp.key$mort.diff
y <- seq(1:nrow(spp.key))
limx <- c(min(x)*1.05, max(x)*1.05)

plot(y ~ x, bty = 'n', xaxt = 'n', yaxt = 'n', cex = NA, ylab = NA, xlab = NA, xlim = limx)
axis(3)
abline (v = 0, lty = 2, col = 'light gray')
#arrows(x - x.sd*1.96, y, x + x.sd*1.96, y, code = 3, angle=90, length=0.0)
points(y ~ x, pch = spp.key$pch, cex = 1.2)
#y-axis species labels.
text(spp.key$species,cex = 0.8, y = y, x = min(limx), adj = 1, xpd = T)
mtext('Mortality Greater \namong Heterospecifics' , side = 3, line = 2, at = limx[1]*0.5, adj = 0.5, cex = 0.7)
mtext('Mortality Greater \namong Conspecifics', side = 3, line = 2, at = limx[2]*0.5, adj = 0.5, cex = 0.7)
mtext('a.', side = 1, line = -1.65, adj = 0.98, cex = 0.9, font = 1)
#legend
