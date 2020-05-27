#Plot species level effects. Differences always performance under EM relative to performance under AM.
rm(list = ls())
source('paths.r')
library(mgcv)
#function to draw randomly from multivariate normal distribution.
rmvn <- function(n,mu,sig) { ## MVN random deviates
  L <- mroot(sig);m <- ncol(L);
  t(mu + L%*%matrix(rnorm(m*n),m,n)) 
}

#load data.----
d <- readRDS(demographic_fits_gam_species.path)
spp.key <- d$spp.key
d$spp.key <- NULL
cov <- readRDS(demographic_fits_gam_ecoregion.path)
cov <- cov$`ATLANTIC HIGHLANDS`$cov

#General covariate data.----
#covariate data.
newdat.am <- c(0,cov)
newdat.em <- c(1,cov)
newdat <- data.frame(rbind(newdat.am, newdat.em))
colnames(newdat)[1] <- 'relEM'


#draw from posterior.----
spp.result <- list()
for(k in 1:length(d)){
  spp <- d[[k]]
  g.mod <- spp$G.mod
  m.mod <- spp$M.mod
  r.mod <- spp$R.mod
  
  #covert X predictors to match knots.
  Xp.g <- predict(g.mod, newdata = newdat, type = 'lpmatrix')
  Xp.m <- predict(m.mod, newdata = newdat, type = 'lpmatrix')
  Xp.r <- predict(r.mod, newdata = newdat, type = 'lpmatrix')
  
  #Draw matrix of correlated predictors from posterior.
  g.par <- rmvn(1000,coef(g.mod),g.mod$Vp)
  m.par <- rmvn(1000,coef(m.mod),m.mod$Vp)
  r.par <- rmvn(1000,coef(r.mod),r.mod$Vp)
  #Multiply predictors by each posterior draw of parameters.
  g.pred <- list()
  m.pred <- list()
  r.pred <- list()
  for(i in 1:1000){
    g.pred[[i]] <- t(Xp.g %*% g.par[i,])
    m.pred[[i]] <- t(Xp.m %*% m.par[i,])
    r.pred[[i]] <- t(Xp.r %*% r.par[i,])
  }
  g.pred <- do.call(rbind, g.pred)
  m.pred <- do.call(rbind, m.pred)
  r.pred <- do.call(rbind, r.pred)
  g.mu <- mean(g.pred[,2] - g.pred[,1])
  g.sd <-   sd(g.pred[,2] - g.pred[,1])
  m.mu <- mean(m.pred[,2] - m.pred[,1])
  m.sd <-   sd(m.pred[,2] - m.pred[,1])
  r.mu <- mean(r.pred[,2] - r.pred[,1])
  r.sd <-   sd(r.pred[,2] - r.pred[,1])
  out <- list(g.mu, g.sd, m.mu, m.sd, r.mu, r.sd)
  names(out) <- c('g.mu','g.sd','m.mu','m.sd','r.mu','r.sd')
  spp.result[[k]] <- out
}
names(spp.result) <- names(d)

#Prepare data for recruitment and mortality plots.----
g.mu <- rep(0, length(spp.result))
g.sd <- rep(0, length(spp.result))
r.mu <- rep(0, length(spp.result))
r.sd <- rep(0, length(spp.result))
m.mu <- rep(0, length(spp.result))
m.sd <- rep(0, length(spp.result))
for(i in 1:length(spp.result)){
  g.mu[i] <- spp.result[[i]]$g.mu
  g.sd[i] <- spp.result[[i]]$g.sd
  m.mu[i] <- spp.result[[i]]$m.mu
  m.sd[i] <- spp.result[[i]]$m.sd
  r.mu[i] <- spp.result[[i]]$r.mu
  r.sd[i] <- spp.result[[i]]$r.sd
}
spp.key$g.mu <- g.mu
spp.key$g.sd <- g.sd
spp.key$m.mu <- m.mu
spp.key$m.sd <- m.sd
spp.key$r.mu <- r.mu
spp.key$r.sd <- r.sd

#setup pch key.
spp.key$pch <- 16
spp.key$pch <- ifelse(spp.key$MYCO_ASSO == 'ECM' & spp.key$gymno == 1, 17, spp.key$pch)
spp.key$pch <- ifelse(spp.key$MYCO_ASSO ==  'AM'                     ,  1, spp.key$pch)
spp.key$pch <- ifelse(spp.key$MYCO_ASSO ==  'AM' & spp.key$gymno == 1,  2, spp.key$pch)

#Begin plots.----
png('test.png', width = 6, height = 8, units = 'in', res = 300)
par(mfrow = c(2,1),mar = c(0.1,6,6,0.1), oma = c(1,1,0,1))
#recruitment.
spp.key <- spp.key[order(spp.key$r.mu),]
x    <- spp.key$r.mu
x.sd <- spp.key$r.sd
y <- seq(1:nrow(spp.key))
limx <- max(abs(c((x - 1.96*x.sd),(x+1.96*x.sd))))
limx <- c(-limx, limx)

plot(y ~ x, bty = 'n', xaxt = 'n', yaxt = 'n', cex = NA, ylab = NA, xlab = NA, xlim = limx)
axis(3)
abline (v = 0, lty = 2, color = 'light gray')
arrows(x - x.sd*1.96, y, x + x.sd*1.96, y, code = 3, angle=90, length=0.0)
points(y ~ x, pch = spp.key$pch, cex = 1.2)
#y-axis species labels.
text(spp.key$name,cex = 0.8, y = y, x = min(limx), adj = 1, xpd = T)
mtext('Recruitment Advantage\namong AM Trees' , side = 3, line = 2, at = limx[1]*0.5, adj = 0.5, cex = 0.8)
mtext('Recruitment Advantage\namong ECM Trees', side = 3, line = 2, at = limx[2]*0.5, adj = 0.5, cex = 0.8)
#legend
legend(x = 1.5, y = 8, legend = c('AM hardwood','AM conifer', 'ECM hardwood','ECM confier'), pch = c(1,2,16,17) , ncol=, bg = NA, bty = 'n', cex = 0.8)



#mortality.
spp.key <- spp.key[order(spp.key$m.mu),]
x    <- spp.key$m.mu
x.sd <- spp.key$m.sd
y <- seq(1:nrow(spp.key))
limx <- max(abs(c((x - 1.96*x.sd),(x+1.96*x.sd))))
limx <- c(-limx, limx)

plot(y ~ x, bty = 'n', xaxt = 'n', yaxt = 'n', cex = NA, ylab = NA, xlab = NA, xlim = limx)
axis(3)
abline (v = 0, lty = 2, color = 'light gray')
arrows(x - x.sd*1.96, y, x + x.sd*1.96, y, code = 3, angle=90, length=0.0)
points(y ~ x, pch = spp.key$pch, cex = 1.2)
#y-axis species labels.
text(spp.key$name,cex = 0.8, y = y, x = min(limx), adj = 1, xpd = T)
mtext('Survival Advantage\namong ECM Trees', side = 3, line = 2, at = limx[1]*0.5, adj = 0.5, cex = 0.8)
mtext('Survival Advantage\namong AM Trees' , side = 3, line = 2, at = limx[2]*0.5, adj = 0.5, cex = 0.8)


#end plot.
dev.off()

