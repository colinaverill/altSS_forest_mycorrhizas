#visualize species level growth responses.
rm(list=ls())
source("paths.r")
source("project_functions/predict_gam_well.r")
library(mgcv)

#set output path.----
output.path <- 'figures/Extended_Data_Fig._5._growth_responses.pdf'

#load and workup functional group level growth.----
N <- 300 #number of posterior draws.
fit <- readRDS(demographic_fits_gam_separate_plus_re_county.path)
env <- fit$env.cov
fit.am <- fit$y.feedback$G.mod.am
fit.em <- fit$y.feedback$G.mod.em

#Estimate impact of relEM on AM vs. EM mortality.
am.am.dat <- c(0,0,25000,30,22,env)
names(am.am.dat)[1:5] <- c('em','relEM','BASAL.plot','stem.density','PREVDIA.cm')
#add county level random effect that will be ignored.
check1 <- fit.am$model$county.ID
check2 <- fit.em$model$county.ID
check <- check1[check1 %in% check2]
am.am.dat <- data.frame(t(am.am.dat))
am.am.dat$county.ID <- as.factor(check[1])

#make the rest of the data sets.
am.em.dat <- am.am.dat
am.em.dat$relEM <- 1
em.em.dat <- am.em.dat
em.em.dat$em <- 1
em.am.dat <- em.em.dat
em.am.dat$relEM <- 0
am.am <- predict_gam_well(fit.am, newdata = am.am.dat, ranef.lab = "county.ID")
am.em <- predict_gam_well(fit.am, newdata = am.em.dat, ranef.lab = "county.ID")
em.am <- predict_gam_well(fit.em, newdata = em.am.dat, ranef.lab = "county.ID")
em.em <- predict_gam_well(fit.em, newdata = em.em.dat, ranef.lab = "county.ID")
#am.am <- predict(fit.am, newdata = am.am.dat, se.fit = T, exclude = c("s(county.ID)"), newdata.guaranteed = T)
#am.em <- predict(fit.am, newdata = am.em.dat, se.fit = T, exclude = c("s(county.ID)"), newdata.guaranteed = T)
#em.am <- predict(fit.em, newdata = em.am.dat, se.fit = T, exclude = c("s(county.ID)"), newdata.guaranteed = T)
#em.em <- predict(fit.em, newdata = em.em.dat, se.fit = T, exclude = c("s(county.ID)"), newdata.guaranteed = T)
y.am.am <- rnorm(N, am.am$fit, am.am$se.fit)
y.am.em <- rnorm(N, am.em$fit, am.em$se.fit)
y.em.am <- rnorm(N, em.am$fit, em.am$se.fit)
y.em.em <- rnorm(N, em.em$fit, em.em$se.fit)
y.am.grow <- c(y.am.am, y.em.am) - am.am.dat$PREVDIA.cm
y.em.grow <- c(y.am.em, y.em.em) - am.am.dat$PREVDIA.cm


#load and workup species level growth.----
d <- readRDS(demographic_fits_gam_species.path)
spp.key <- d$spp.key
d$spp.key <- NULL
g.mu <- rep(0, length(d))
g.sd <- rep(0, length(d))
r.mu <- rep(0, length(d))
r.sd <- rep(0, length(d))
m.mu <- rep(0, length(d))
m.sd <- rep(0, length(d))
for(i in 1:length(d)){
  g.mu[i] <- d[[i]]$post.contrast$g.mu
  g.sd[i] <- d[[i]]$post.contrast$g.sd
  m.mu[i] <- d[[i]]$post.contrast$m.mu
  m.sd[i] <- d[[i]]$post.contrast$m.sd
  r.mu[i] <- d[[i]]$post.contrast$r.mu
  r.sd[i] <- d[[i]]$post.contrast$r.sd
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

#setup x positions for first panel.-----
jit <- 0.05
x.pos <- c(0.5,1,1.5,2)
x1 <- rnorm(N,x.pos[1],jit)
x2 <- rnorm(N,x.pos[2],jit)
x3 <- rnorm(N,x.pos[3],jit)
x4 <- rnorm(N,x.pos[4],jit)
x.am <- c(x1,x2)
x.em <- c(x3,x4)

#png save line.----
#png(output.path, width = 11, height = 5, units = 'in', res = 300)
pdf(output.path, width = 11, height = 5)

#panels.----
par(mfrow = c(1,2),
    mar = c(1.5,6.5,5,1),
    oma = c(1,1, 1,1))

#functional group level  level survival.-----
cols <- c('green','purple')
col.1 <- c(rep(cols[1],N),rep(cols[2],N))
limx <- c(min(c(x.am,x.em)), max(c(x.am,x.em)))
#limy <- c(min(c(y.am.grow,y.em.grow)), max(c(y.am.grow,y.em.grow)))
limy <- c(0, max(c(y.am.grow,y.em.grow)))
#plot.
plot(y.am.grow ~ x.am, pch = 16, xlim = limx, ylim = limy, cex = 0.3, bty = 'l', ylab = NA, xlab = NA, xaxt = 'n', col = col.1)
points(y.em.grow ~ x.em, pch = 16, cex = 0.3, col = col.1)
abline(v = 1.25, lwd=2, lty = 3, col = 'light gray')
#labels.
mtext('Predicted Diameter Increment (cm)', side = 2, cex = 1.0, line = 2.5)
legend('topleft',legend = c('All EM trees','All AM trees'), 
       pch = 16, col = rev(cols), bty = 'n', cex = 0.9,
       x.intersp = .75, xpd = T, 
       horiz = F)
mtext('AM Forest',side = 1, line = 1, adj = 0.175, cex = 1)
mtext('EM Forest',side = 1, line = 1, adj = 0.825, cex = 1)
mtext('A', side = 1, line = -1.65, adj = 0.98, cex = 0.9, font = 1)

#species level growth.----
spp.key <- spp.key[order(spp.key$g.mu),]
x    <- spp.key$g.mu
x.sd <- spp.key$g.sd
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
mtext('Growth Advantage\namong AM Trees' , side = 3, line = 2, at = limx[1]*0.5, adj = 0.5, cex = 0.7)
mtext('Growth Advantage\namong EM Trees', side = 3, line = 2, at = limx[2]*0.5, adj = 0.5, cex = 0.7)
mtext('B', side = 1, line = -1.65, adj = 0.98, cex = 0.9, font = 1)
#legend
legend(x = 0.5, y = 8, legend = c('AM hardwood','AM conifer', 'EM hardwood','EM confier'), pch = c(1,2,16,17) , ncol=, bg = NA, bty = 'n', cex = 0.8)

#end plot.----
dev.off()
