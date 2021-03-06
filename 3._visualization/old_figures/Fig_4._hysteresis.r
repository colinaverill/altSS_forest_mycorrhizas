#AM - EM hysteresis conceptual figure and simulation result.
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
#arrows(7.5, 0.7, 9.25, 0.3, lwd = 2, length = 0.1)
#arrows(7.5, 0.29, 5.75, 0.69, lwd = 2, length = 0.1)
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
#arrows(10.0, 0.7, 11.75, 0.3, lwd = 2, length = 0.1)
#arrows(5.0, 0.29, 3.25, 0.69, lwd = 2, length = 0.1)
mtext('B',side=1,line=-1.5,adj=0.05, font = 2)
text(x =  1.5, y = -0.1, 'low' , xpd = NA, cex = 1.2)
text(x = 13.5, y = -0.1, 'high', xpd = NA, cex = 1.2)
mtext('Nitrogen Deposition',side = 1, line = 2.5 )
mtext('H1: Alternative Stable States', side = 3, adj = 0.05, line = -1.5)

#legend.
legend(x= 9.5,y= 0.95,legend=c("initially EM","initially AM"),lty=c(2,3),lwd=c(3,3), 
       box.lwd=0, bg="transparent",cex=1)

#hysteresis data load and workup.----
#load data.
d <- readRDS(initial_condition_hysteresis_simulation.path)
#grab results after 200 years, make these final plot tables to work with other code.
time.step <- 41 #41 time steps = 200 years.
for(i in 1:length(d$all.em$nul)){
  d$all.em$nul    [[i]]$plot.table <- d$all.em$nul    [[i]]$super.table[[time.step]]
  d$all.em$alt.GRM[[i]]$plot.table <- d$all.em$alt.GRM[[i]]$super.table[[time.step]]
  d$all.am$nul    [[i]]$plot.table <- d$all.am$nul    [[i]]$super.table[[time.step]]
  d$all.am$alt.GRM[[i]]$plot.table <- d$all.am$alt.GRM[[i]]$super.table[[time.step]]
}
n.tot <- nrow(d$all.em$alt.GRM$l1$plot.table)
n.lev <- length(d$all.am$alt.GRM)

#get number of EM and AM plots in feedback simulations (>90%) and bootstrap CI.
#starting from EM state.
em.ramp <- list()
for(i in 1:length(d$all.em$alt.GRM)){
  plot.tab <- d$all.em$alt.GRM[[i]]$plot.table
  N.em     <- nrow(plot.tab[plot.tab$relEM > 0.9,])
  N.am     <- nrow(plot.tab[plot.tab$relEM < 0.1,])
  #Get boostrap CI on mean.
  boot.dat <-list()
  for(k in 1:1000){
    dat <- plot.tab[sample(nrow(plot.tab), nrow(plot.tab), replace = T),]
    boot.em <- nrow(dat[dat$relEM > 0.9,])
    boot.am <- nrow(dat[dat$relEM < 0.1,])
    boot.dat[[k]] <- c(boot.em,boot.am)
  }
  boot.dat <- do.call(rbind, boot.dat)
  em.95 <- quantile(boot.dat[,1], probs = c(0.025, 0.975))
  am.95 <- quantile(boot.dat[,2], probs = c(0.025, 0.975))
  #return output.
  em.ramp[[i]] <- c(N.em, N.am, em.95,am.95)
}
em.ramp <- data.frame(do.call(rbind, em.ramp))
colnames(em.ramp) <- c('N.em','N.am','em.lo95','em.hi95','am.lo95','am.hi95')
#starting from AM state.
am.ramp <- list()
for(i in 1:length(d$all.am$alt.GRM)){
  plot.tab <- d$all.am$alt.GRM[[i]]$plot.table
  N.em     <- nrow(plot.tab[plot.tab$relEM > 0.9,])
  N.am     <- nrow(plot.tab[plot.tab$relEM < 0.1,])
  #Get boostrap CI on mean.
  boot.dat <-list()
  for(k in 1:1000){
    dat <- plot.tab[sample(nrow(plot.tab), nrow(plot.tab), replace = T),]
    boot.em <- nrow(dat[dat$relEM > 0.9,])
    boot.am <- nrow(dat[dat$relEM < 0.1,])
    boot.dat[[k]] <- c(boot.em,boot.am)
  }
  boot.dat <- do.call(rbind, boot.dat)
  em.95 <- quantile(boot.dat[,1], probs = c(0.025, 0.975))
  am.95 <- quantile(boot.dat[,2], probs = c(0.025, 0.975))
  #return output.
  am.ramp[[i]] <- c(N.em, N.am, em.95,am.95)
}
am.ramp <- data.frame(do.call(rbind, am.ramp))
colnames(am.ramp) <- c('N.em','N.am','em.lo95','em.hi95','am.lo95','am.hi95')

#convert to frequencies, assign n levels, rename as up and down.
#renaming
up   <- em.ramp[,c('N.em','em.lo95','em.hi95')]
down <- am.ramp[,c('N.em','em.lo95','em.hi95')]
colnames(  up) <- c('n','lo95','hi95')
colnames(down) <- c('n','lo95','hi95')
#convert to frequencies.
  up <-   up / n.tot
down <- down / n.tot
#Assign ndep.
up  $ndep <- c(1:n.lev)
down$ndep <- c(1:n.lev)

#Panel 3.  hysteresis simulation ramp up / ramp down.----
par(mar = c(0,3,0,0))
cols <- c('purple','green')
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
#arrows(x0 = 5, y0 = .375, x1 = 9.0, y1 = .225, length = 0.1, lwd = 1.5)
#arrows(x0 = 7, y0 = 0.08, x1 = 3, y1 = 0.16, length = 0.1, lwd = 1.5)

#legend.
legend(x=8.5,y= 0.8,legend=c("initially EM","initially AM"),lty=1,lwd=2, col = cols, 
         box.lwd=0, bg="transparent",cex= 1)

#outer labels.
#mtext('Number EM Dominated Forests', side =2, cex=1.2, line = 2.5)
mtext('Demographic Simulation', side = 3, adj = 0.05, line = -1.5)
mtext('Nitrogen Deposition', side = 1, line = 2.5)
mtext('C', side = 1, line = -1.5, adj = 0.05, font = 2)
mtext(expression(paste("kg N ha"^"-1"," yr"^"-1",sep="")), side = 1, line = 4.5,  cex = 0.8)


#end plot.----
dev.off()

