#Plot recruitment and survival of all species together as AM vs. EM and at species level.
#Convert mortality to survival.
rm(list=ls())
library(mgcv)
library(data.table)
source('paths.r')
#function to draw randomly from multivariate normal distribution.-----
rmvn <- function(n,mu,sig) { ## MVN random deviates
  L <- mroot(sig);m <- ncol(L);
  t(mu + L%*%matrix(rnorm(m*n),m,n)) 
}
#define number of posterior draws for plotting.----
N <- 300

#load, workup mortality (survival) data.----
#Load fits and data.
fit <- readRDS(demographic_fits.path)
env <- fit$env.cov
fit <- fit$y.feedback$M.mod

#Estimate impact of relEM on AM vs. EM mortality.
am.am.dat <- c(0,0,25000,30,22,env)
names(am.am.dat)[1:5] <- c('em','relEM','BASAL.plot','stem.density','PREVDIA.cm')
am.em.dat <- am.am.dat
am.em.dat['relEM'] <- 1
em.em.dat <- am.em.dat
em.em.dat['em'] <- 1
em.am.dat <- em.em.dat
em.am.dat['relEM'] <- 0
am.am <- predict(fit, newdata = data.frame(t(am.am.dat)), se.fit = T)
am.em <- predict(fit, newdata = data.frame(t(am.em.dat)), se.fit = T)
em.am <- predict(fit, newdata = data.frame(t(em.am.dat)), se.fit = T)
em.em <- predict(fit, newdata = data.frame(t(em.em.dat)), se.fit = T)
y.am.am <- rnorm(N, am.am$fit, am.am$se.fit)
y.am.em <- rnorm(N, am.em$fit, am.em$se.fit)
y.em.am <- rnorm(N, em.am$fit, em.am$se.fit)
y.em.em <- rnorm(N, em.em$fit, em.em$se.fit)
y.am <- exp(c(y.am.am, y.em.am)) #undo log link to get to scale of recruitment.
y.em <- exp(c(y.am.em, y.em.em))

#convert to survival, rather than mortality.
y.am.surv <- 1 - y.am
y.em.surv <- 1 - y.em

#load, workup recruitment data.----
#Load fits and data.
fit <- readRDS(demographic_fits.path)
env <- fit$env.cov
fit <- fit$y.feedback
d <- data.table(readRDS(Product_2.subset.path))

#get mean and se.
am.am.dat <- c(0,0,25000,30,env)
names(am.am.dat)[1:4] <- c('am.density','relEM','BASAL.plot','stem.density')
am.em.dat <- am.am.dat
am.em.dat['relEM'] <- 1
em.em.dat <- am.em.dat
names(em.em.dat)[1] <- 'em.density'
em.am.dat <- em.em.dat
em.am.dat['relEM'] <- 0
am.am <- predict(fit$R.mod.am, newdata = data.frame(t(am.am.dat)), se.fit = T)
am.em <- predict(fit$R.mod.am, newdata = data.frame(t(am.em.dat)), se.fit = T)
em.am <- predict(fit$R.mod.em, newdata = data.frame(t(em.am.dat)), se.fit = T)
em.em <- predict(fit$R.mod.em, newdata = data.frame(t(em.em.dat)), se.fit = T)
N <- 300
y.am.am <- rnorm(N, am.am$fit, am.am$se.fit)
y.am.em <- rnorm(N, am.em$fit, am.em$se.fit)
y.em.am <- rnorm(N, em.am$fit, em.am$se.fit)
y.em.em <- rnorm(N, em.em$fit, em.em$se.fit)
y.am.recr <- exp(c(y.am.am, y.em.am)) #undo log link to get to scale of recruitment.
y.em.recr <- exp(c(y.am.em, y.em.em))

#load species level recruitment and mortality data.----
d <- readRDS(demographic_fits_gam_species.path)
spp.key <- d$spp.key
d$spp.key <- NULL
cov <- readRDS(demographic_fits_gam_ecoregion.path)
cov <- cov$`ATLANTIC HIGHLANDS`$cov

#General spp-level covariate data.
newdat.am <- c(0,cov)
newdat.em <- c(1,cov)
newdat <- data.frame(rbind(newdat.am, newdat.em))
colnames(newdat)[1] <- 'relEM'


#draw from species-level posteriors.----
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
  #convert mortality to survival.
  m.pred <- boot::logit(1 - boot::inv.logit(m.pred))
  
  g.mu <-     mean(g.pred[,2] - g.pred[,1])
  g.sd <-       sd(g.pred[,2] - g.pred[,1])
  g.95 <- quantile(g.pred[,2] - g.pred[,1], probs = c(0.025, 0.975))
  m.mu <-     mean(m.pred[,2] - m.pred[,1])
  m.sd <-       sd(m.pred[,2] - m.pred[,1])
  m.95 <- quantile(m.pred[,2] - m.pred[,1], probs = c(0.025, 0.975))
  r.mu <- mean(r.pred[,2] - r.pred[,1])
  r.sd <-   sd(r.pred[,2] - r.pred[,1])
  r.95 <- quantile(r.pred[,2] - r.pred[,1], probs = c(0.025, 0.975))
  out <- list(g.mu, g.sd, g.95, m.mu, m.sd, m.95, r.mu, r.sd, r.95)
  names(out) <- c('g.mu','g.sd','g.95','m.mu','m.sd','m.95','r.mu','r.sd','r.95')
  spp.result[[k]] <- out
}
names(spp.result) <- names(d)

#Organize species level data for recruitment and survival plots.----
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


#setup x positions for first two panels.-----
jit <- 0.05
x.pos <- c(0.5,1,1.5,2)
x1 <- rnorm(N,x.pos[1],jit)
x2 <- rnorm(N,x.pos[2],jit)
x3 <- rnorm(N,x.pos[3],jit)
x4 <- rnorm(N,x.pos[4],jit)
x.am <- c(x1,x2)
x.em <- c(x3,x4)

#Setup panels.----
par(mfrow = c(2,2),
    mar = c(1,6,5,1),
    oma = c(1, 5, 1,1))

#plot survival effects.----
cols <- c('green','purple')
col.1 <- c(rep(cols[1],N),rep(cols[2],N))
limx <- c(min(c(x.am,x.em)), max(c(x.am,x.em)))
limy <- c(min(c(y.am.surv,y.em.surv)), 1)
#plot.
plot(y.am.surv ~ x.am, pch = 16, xlim = limx, ylim = limy, cex = 0.3, bty = 'l', ylab = NA, xlab = NA, xaxt = 'n', col = col.1)
points(y.em.surv ~ x.em, pch = 16, cex = 0.3, col = col.1)
abline(v = 1.25, lwd=2, lty = 3, col = 'light gray')
#labels.
mtext('Survival Probability', side = 2, cex = 1.0, line = 2.5)
legend('topleft',legend = c('EM tree','AM tree'), 
       pch = 16, col = rev(cols), bty = 'n', cex = 0.9,
       x.intersp = .75, xpd = T, 
       horiz = F)
mtext('AM Forest',side = 1, line = 1, adj = 0.125, cex = 1)
mtext('EM Forest',side = 1, line = 1, adj = 0.9, cex = 1)
mtext('A', side = 1, line = -1.75, adj = 0.95, cex = 1.0, font = 2)


#plot recruitment effects.----
cols <- c('green','purple')
col.1 <- c(rep(cols[1],N),rep(cols[2],N))
labx <- c('AM in EM Forest','AM in AM Forest','EM in AM Forest','EM in EM Forest')
limx <- c(min(c(x.am,x.em)), max(c(x.am,x.em)))
limy <- c(0, max(c(y.am.recr,y.em.recr)))
#plot.
  plot(y.am.recr ~ x.am, pch = 16, xlim = limx, ylim = limy, cex = 0.3, bty = 'l', ylab = NA, xlab = NA, xaxt = 'n', col = col.1)
points(y.em.recr ~ x.em, pch = 16, cex = 0.3, col = col.1)
abline(v = 1.25, lwd=2, lty = 3, col = 'light gray')
#labels.
lab <- expression(paste('Recruitment - Poisson ',lambda))
mtext(lab, side = 2, cex = 1.0, line = 2.5)
mtext('AM Forest',side = 1, line = 1, adj = 0.125, cex = 1)
mtext('EM Forest',side = 1, line = 1, adj = 0.9, cex = 1)
mtext('B', side = 1, line = -1.75, adj = 0.95, cex = 1.0, font = 2)

#species level mortality.-----
spp.key <- spp.key[order(spp.key$m.mu),]
x    <- spp.key$m.mu
x.sd <- spp.key$m.sd
y <- seq(1:nrow(spp.key))
limx <- max(abs(c((x - 1.96*x.sd),(x+1.96*x.sd))))
limx <- c(-limx, limx)

plot(y ~ x, bty = 'n', xaxt = 'n', yaxt = 'n', cex = NA, ylab = NA, xlab = NA, xlim = limx)
axis(3)
abline (v = 0, lty = 2, col = 'light gray')
arrows(x - x.sd, y, x + x.sd, y, code = 3, angle=90, length=0.0)
points(y ~ x, pch = spp.key$pch, cex = 1.2)
#y-axis species labels.
text(spp.key$name,cex = 0.8, y = y, x = min(limx), adj = 1, xpd = T)
mtext('Survival Advantage\namong AM Trees' , side = 3, line = 2, at = limx[1]*0.5, adj = 0.5, cex = 0.8)
mtext('Survival Advantage\namong ECM Trees', side = 3, line = 2, at = limx[2]*0.5, adj = 0.5, cex = 0.8)

#species level recruitment.----
spp.key <- spp.key[order(spp.key$r.mu),]
x    <- spp.key$r.mu
x.sd <- spp.key$r.sd
y <- seq(1:nrow(spp.key))
limx <- max(abs(c((x - 1.96*x.sd),(x+1.96*x.sd))))
limx <- c(-limx, limx)

plot(y ~ x, bty = 'n', xaxt = 'n', yaxt = 'n', cex = NA, ylab = NA, xlab = NA, xlim = limx)
axis(3)
abline (v = 0, lty = 2, col = 'light gray')
arrows(x - x.sd*1.96, y, x + x.sd*1.96, y, code = 3, angle=90, length=0.0)
points(y ~ x, pch = spp.key$pch, cex = 1.2)
#y-axis species labels.
text(spp.key$name,cex = 0.8, y = y, x = min(limx), adj = 1, xpd = T)
mtext('Recruitment Advantage\namong AM Trees' , side = 3, line = 2, at = limx[1]*0.5, adj = 0.5, cex = 0.8)
mtext('Recruitment Advantage\namong ECM Trees', side = 3, line = 2, at = limx[2]*0.5, adj = 0.5, cex = 0.8)
#legend
legend(x = 1.5, y = 8, legend = c('AM hardwood','AM conifer', 'ECM hardwood','ECM confier'), pch = c(1,2,16,17) , ncol=, bg = NA, bty = 'n', cex = 0.8)



#end plots.----
#dev.off()
