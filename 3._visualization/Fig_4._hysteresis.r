#AM to EM hysteresis conceptual figure
rm(list=ls())
source('paths.r')

#set output path.----
output.path <- Fig_4.path

#generate conceptual figure data.-----
x  <- seq(0,14,length.out=100)
y1 <- (-1/ (1+ exp(-(x-10.0))))+1
y2 <- (-1/ (1+ exp(-(x- 5.0))))+1
y3 <- (-1/ (1+ exp(-(x- 7.5))))+1

#setup plot save destination, general figure parameters.----
png(filename=output.path, width=10, height=4, units='in', res=300)
#Set number of panels, overall margins.
par(mfrow=c(1,3), oma=c(5.5,5,2,1), mar=c(0,0,0,0))
#set vertical adjustment to add room for titles.
adjust <- 1.07

#Panel 1. null hypothesis.----
limy <- max(y3) * adjust
plot(y3~x,ylab=NA,xlab=NA,cex=0,cex.lab=1.2, xaxt='n', yaxt = 'n', bty = 'l', ylim = c(0, limy))
lines(smooth.spline(y3~x,spar=0.35), lwd=3,lty=1)
arrows(7.5, 0.7, 9.25, 0.3, lwd = 2, length = 0.1)
arrows(7.5, 0.29, 5.75, 0.69, lwd = 2, length = 0.1)
mtext('A',side=1,line=-1.5,adj=0.05, font = 2)
#mtext('Probability EM dominated',side=2,line = 2.2)
mtext('Frequency Ectomycorrhizal Forests',side=2,line = 3)
text(x = -1.5, y = 1, 'high', xpd = NA, cex = 1.2)
text(x = -1.5, y = 0.025, 'low' , xpd = NA, cex = 1.2)
text(x =  1.5, y = -0.1, 'low' , xpd = NA, cex = 1.2)
text(x = 13.5, y = -0.1, 'high', xpd = NA, cex = 1.2)
mtext('Nitrogen Deposition',side = 1, line = 2.5 )
mtext('H0: Environmental Filtering', side = 3, adj = 0.05, line = -1.5)

#Panel 2. hysteresis hypothesis.----
plot(y1~x,ylab=NA,xlab=NA,cex=0,cex.lab=1.2, yaxt='n',xaxt='n', bty = 'l', ylim = c(0, limy))
lines(smooth.spline(y1~x,spar=0.35), lwd=3,lty=2)
lines(smooth.spline(y2~x,spar=0.35), lwd=3,lty=3)
arrows(10.0, 0.7, 11.75, 0.3, lwd = 2, length = 0.1)
arrows(5.0, 0.29, 3.25, 0.69, lwd = 2, length = 0.1)
mtext('B',side=1,line=-1.5,adj=0.05, font = 2)
text(x =  1.5, y = -0.1, 'low' , xpd = NA, cex = 1.2)
text(x = 13.5, y = -0.1, 'high', xpd = NA, cex = 1.2)
mtext('Nitrogen Deposition',side = 1, line = 2.5 )
mtext('H1: Alt. Stable States', side = 3, adj = 0.05, line = -1.5)

#legend.
#legend(x=8.5,y= 1.05,legend=c("EM to AM","AM to EM"),lty=c(2,3),lwd=c(3,3), box.lwd=0, bg="transparent",cex=0.9)

#Lower label.
#mtext('Nitrogen Deposition', side=1, out=T, cex=1.5, line = 3)

#Hysteresis data load and workup.----
#load data.
d <- readRDS(factorial_hysteresis_simulation.path)
#total number of plots simulated.
n.tot <- nrow(d$ramp.up$nul$l1$plot.table)

#Pick your models.
up <- d$ramp.up  $alt.GRM
down <- d$ramp.down$alt.GRM
n.lev <- length(up)

#count number of EM plots across simulations, and get a 95% bootstrap CI.
#Ramp up N dep estimates.
em.up    <- list()
em.up.CI <- list()
n.straps <- 1000
for(i in 1:length(up)){
  #grab plot table relEM vector
  plot.vec <- up[[i]]$plot.table$relEM
  em.up[[i]]  <- length(plot.vec[plot.vec > 0.9])
  #perform bootstrap.
  boot.out <- list()
  for(k in 1:n.straps){
    plot.sample <- sample(plot.vec, size = length(plot.vec), replace = T)
    boot.out[[k]] <- length(plot.sample[plot.sample > 0.9])
  }
  boot.out <- unlist(boot.out)
  #calculate 95% CI.
  em.up.CI[[i]] <- quantile(boot.out, probs = c(0.025, 0.975))
}
em.up <- unlist(em.up)
em.up.CI <- do.call(rbind, em.up.CI)
up <- data.frame(cbind(em.up, em.up.CI))
colnames(up) <- c('n','lo95','hi95')
up$ndep <- c(1:n.lev)

#convert to frequency
up$n    <- up$n    / n.tot
up$lo95 <- up$lo95 / n.tot
up$hi95 <- up$hi95 / n.tot

#Ramp down Ndep estimates.
em.down    <- list()
em.down.CI <- list()
n.straps <- 1000
for(i in 1:length(down)){
  #grab plot table relEM vector
  plot.vec <- down[[i]]$plot.table$relEM
  em.down[[i]]  <- length(plot.vec[plot.vec > 0.9])
  #perform bootstrap.
  boot.out <- list()
  for(k in 1:n.straps){
    plot.sample <- sample(plot.vec, size = length(plot.vec), replace = T)
    boot.out[[k]] <- length(plot.sample[plot.sample > 0.9])
  }
  boot.out <- unlist(boot.out)
  #calculate 95% CI.
  em.down.CI[[i]] <- quantile(boot.out, probs = c(0.025, 0.975))
}
em.down <- unlist(em.down)
em.down.CI <- do.call(rbind, em.down.CI)
down <- data.frame(cbind(em.down, em.down.CI))
colnames(down) <- c('n','lo95','hi95')
down$ndep <- c(n.lev:1)

#convert to frequency
down$n    <- down$n    / n.tot
down$lo95 <- down$lo95 / n.tot
down$hi95 <- down$hi95 / n.tot

#Panel 3.  hysteresis simulation ramp up / ramp down.----
par(mar = c(0,3,0,0))
cols <- c('purple','light pink')
trans <- 0.3
#plot ramp up, transitioning away from EM dominated forests.
color <- cols[1]
limy <- max(up$hi95) * adjust
#plot(up$n ~ up$ndep, bty = 'l', cex = 0, ylab = NA, xlab = NA, ylim = c(0.9*min(c(up$lo95, down$lo95)),max(up$n)*1.05))
plot(up$n ~ up$ndep, bty = 'l', cex = 0, ylab = NA, xlab = NA, ylim = c(0, limy))
lines(smooth.spline(up$n ~ up$ndep, spar = .1), lwd = 3, col = color)
polygon(c(up$ndep, rev(up$ndep)),c(up$hi95, rev(up$lo95)), col=adjustcolor(color, trans), lty=0)

#plot ramp down.
color <- cols[2]
lines(smooth.spline(down$n ~ down$ndep, spar = 0.1), lwd = 2, col = adjustcolor(color, 0.6))
polygon(c(down$ndep, rev(down$ndep)),c(down$hi95, rev(down$lo95)), col=adjustcolor(color, trans), lty=0)

#arrows
arrows(x0 = 5, y0 = .375, x1 = 9.0, y1 = .225, length = 0.1, lwd = 1.5)
arrows(x0 = 7, y0 = 0.08, x1 = 3, y1 = 0.16, length = 0.1, lwd = 1.5)

#outer labels.
#mtext('Number EM Dominated Forests', side =2, cex=1.2, line = 2.5)
mtext('Demographic Simulation', side = 3, adj = 0.05, line = -1.5)
mtext('Nitrogen Deposition', side = 1, line = 2.5)
mtext('C', side = 1, line = -1.5, adj = 0.05, font = 2)
mtext(expression(paste("kg N ha"^"-1"," yr"^"-1",sep="")), side = 1, line = 4.5,  cex = 0.8)


#end plot.----
dev.off()

