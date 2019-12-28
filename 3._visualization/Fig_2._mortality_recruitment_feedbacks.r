#5-panel to show recruitment and mortality positive feedbacks, demo simulation altnerative stable states, and hysteresis.
rm(list=ls())
source('paths.r')
library(mgcv)

#set output path.----
output.path <- Fig_2.path

#setup 2-panel and save line.----
png(output.path, width = 8, height = 4, units = 'in', res = 300)
#Setup 2-panel 
par(mfrow = c(1,2),
    mar = c(2.5,3.5,1,1))
#load, workup and plot positive mortality feedbacks.----
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
N <- 300
y.am.am <- rnorm(N, am.am$fit, am.am$se.fit)
y.am.em <- rnorm(N, am.em$fit, am.em$se.fit)
y.em.am <- rnorm(N, em.am$fit, em.am$se.fit)
y.em.em <- rnorm(N, em.em$fit, em.em$se.fit)
y.am <- exp(c(y.am.am, y.em.am)) #undo log link to get to scale of recruitment.
y.em <- exp(c(y.am.em, y.em.em))

#setup x positions.
jit <- 0.05
x.pos <- c(0.5,1,1.5,2)
x1 <- rnorm(N,x.pos[1],jit)
x2 <- rnorm(N,x.pos[2],jit)
x3 <- rnorm(N,x.pos[3],jit)
x4 <- rnorm(N,x.pos[4],jit)
x.am <- c(x1,x2)
x.em <- c(x3,x4)

#plot mortality effects.----
#par(mar = c(4,5,1,1))
#cols <- c('#00acd9','#cfe83c')
cols <- c('green','purple')
col.1 <- c(rep(cols[1],N),rep(cols[2],N))
limx <- c(min(c(x.am,x.em)), max(c(x.am,x.em)))
limy <- c(0, max(c(y.am,y.em)))
#plot.
plot(y.am ~ x.am, pch = 16, xlim = limx, ylim = limy, cex = 0.3, bty = 'l', ylab = NA, xlab = NA, xaxt = 'n', col = col.1)
points(y.em ~ x.em, pch = 16, cex = 0.3, col = col.1)
abline(v = 1.25, lwd=2, lty = 3, col = 'light gray')
#labels.
mtext('Mortality Probability', side = 2, cex = 1.0, line = 2.5)
legend('bottomleft',legend = c('ectomycorrhizal','arbuscular mycorrhizal'), 
       pch = 16, col = rev(cols), bty = 'n', cex = 0.9,
       x.intersp = .75, xpd = T, 
       horiz = F)
mtext('AM Forest',side = 1, line = 1, adj = 0.125, cex = 1)
mtext('EM Forest',side = 1, line = 1, adj = 0.9, cex = 1)
mtext('A', side = 1, line = -1.75, adj = 0.95, cex = 1.0, font = 2)

#load, workup and plot positive recruitment feedbacks.----
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
y.am <- exp(c(y.am.am, y.em.am)) #undo log link to get to scale of recruitment.
y.em <- exp(c(y.am.em, y.em.em))
#setup x positions.
jit <- 0.05
x.pos <- c(0.5,1,1.5,2)
x1 <- rnorm(N,x.pos[1],jit)
x2 <- rnorm(N,x.pos[2],jit)
x3 <- rnorm(N,x.pos[3],jit)
x4 <- rnorm(N,x.pos[4],jit)
x.am <- c(x1,x2)
x.em <- c(x3,x4)

#plot recruitment effects.----
#par(mar = c(4,5,1,1))
#cols <- c('#00acd9','#cfe83c')
cols <- c('green','purple')
col.1 <- c(rep(cols[1],N),rep(cols[2],N))
labx <- c('AM in EM Forest','AM in AM Forest','EM in AM Forest','EM in EM Forest')
limx <- c(min(c(x.am,x.em)), max(c(x.am,x.em)))
limy <- c(0, max(c(y.am,y.em)))
#plot.
plot(y.am ~ x.am, pch = 16, xlim = limx, ylim = limy, cex = 0.3, bty = 'l', ylab = NA, xlab = NA, xaxt = 'n', col = col.1)
points(y.em ~ x.em, pch = 16, cex = 0.3, col = col.1)
abline(v = 1.25, lwd=2, lty = 3, col = 'light gray')
#labels.
lab <- expression(paste('Recruitment - Poisson ',lambda))
mtext(lab, side = 2, cex = 1.0, line = 2.5)
mtext('AM Forest',side = 1, line = 1, adj = 0.125, cex = 1)
mtext('EM Forest',side = 1, line = 1, adj = 0.9, cex = 1)
mtext('B', side = 1, line = -1.75, adj = 0.95, cex = 1.0, font = 2)

#end plots.----
dev.off()
