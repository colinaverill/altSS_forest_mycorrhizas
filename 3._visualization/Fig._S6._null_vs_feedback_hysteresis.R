#AM - EM hysteresis conceptual figure and simulation result.
#code below flips am vs. em ramps, because i flipped the initial conditions in the actual simulation script.
rm(list=ls())
source('paths.r')

#set output path.----
output.path <- 'figures/Extended_Data_Fig._6._hysteresis_null_vs._feedback.pdf'

#setup plot save destination, general figure parameters.----
#jpeg(filename=output.path, width=8, height=5, units='in', res=300)
pdf(file=output.path, width=8, height=5)

#Set number of panels, overall margins.-----
par(mfrow=c(1,2), oma=c(5.5,3,2,1), mar=c(0,0,0,0))
#set vertical adjustment to add room for titles.
adjust <- 1.07

#hysteresis data load and workup feedback results.----
#load data.
d <- readRDS(initial_condition_hysteresis_simulation.path)
#grab results after 200 years, make these final plot tables to work with other code.
time.step <- 81 #41 time steps = 200 years.
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
  relEM    <- mean(plot.tab$relEM)
  #Get boostrap CI on mean.
  boot.dat <-list()
  for(k in 1:1000){
    dat <- plot.tab[sample(nrow(plot.tab), nrow(plot.tab), replace = T),]
    boot.em    <- nrow(dat[dat$relEM > 0.9,])
    boot.am    <- nrow(dat[dat$relEM < 0.1,])
    boot.relEM <- mean(dat$relEM)
    boot.dat[[k]] <- c(boot.relEM,boot.em,boot.am)
  }
  boot.dat <- do.call(rbind, boot.dat)
  relEM.95 <- quantile(boot.dat[,1], probs = c(0.025, 0.975))
  em.95 <- quantile(boot.dat[,2], probs = c(0.025, 0.975))
  am.95 <- quantile(boot.dat[,3], probs = c(0.025, 0.975))
  
  #return output.
  em.ramp[[i]] <- c(relEM,N.em,N.am,relEM.95,em.95,am.95)
}
em.ramp <- data.frame(do.call(rbind, em.ramp))
colnames(em.ramp) <- c('relEM','N.em','N.am','relEM.lo95','relEM.hi95','em.lo95','em.hi95','am.lo95','am.hi95')

#starting from AM state.
am.ramp <- list()
for(i in 1:length(d$all.am$alt.GRM)){
  plot.tab <- d$all.am$alt.GRM[[i]]$plot.table
  N.em     <- nrow(plot.tab[plot.tab$relEM > 0.9,])
  N.am     <- nrow(plot.tab[plot.tab$relEM < 0.1,])
  relEM    <- mean(plot.tab$relEM)
  #Get boostrap CI on mean.
  boot.dat <-list()
  for(k in 1:1000){
    dat <- plot.tab[sample(nrow(plot.tab), nrow(plot.tab), replace = T),]
    boot.em    <- nrow(dat[dat$relEM > 0.9,])
    boot.am    <- nrow(dat[dat$relEM < 0.1,])
    boot.relEM <- mean(dat$relEM)
    boot.dat[[k]] <- c(boot.relEM,boot.em,boot.am)
  }
  boot.dat <- do.call(rbind, boot.dat)
  relEM.95 <- quantile(boot.dat[,1], probs = c(0.025, 0.975))
  em.95 <- quantile(boot.dat[,2], probs = c(0.025, 0.975))
  am.95 <- quantile(boot.dat[,3], probs = c(0.025, 0.975))
  #return output.
  am.ramp[[i]] <- c(relEM,N.em,N.am,relEM.95,em.95,am.95)
}
am.ramp <- data.frame(do.call(rbind, am.ramp))
colnames(am.ramp) <- c('relEM','N.em','N.am','relEM.lo95','relEM.hi95','em.lo95','em.hi95','am.lo95','am.hi95')

#Get EM dominated forest frequencies.
#renaming
#up   <- em.ramp[,c('N.em','em.lo95','em.hi95')]
#down <- am.ramp[,c('N.em','em.lo95','em.hi95')]
#colnames(  up) <- c('n','lo95','hi95')
#colnames(down) <- c('n','lo95','hi95')
#convert to frequencies.
#  up <-   up / n.tot
#down <- down / n.tot

#Use average relative abundance EM trees instead.
up <- em.ramp[,c('relEM','relEM.lo95','relEM.hi95')]
down <- am.ramp[,c('relEM','relEM.lo95','relEM.hi95')]
colnames(  up) <- c('n','lo95','hi95')
colnames(down) <- c('n','lo95','hi95')

#Assign ndep.
up  $ndep <- c(1:n.lev)
down$ndep <- c(1:n.lev)
  up.feed <- up
down.feed <- down


#workup null results.----
#get number of EM and AM plots in feedback simulations (>90%) and bootstrap CI.
#starting from EM state.
em.ramp <- list()
for(i in 1:length(d$all.em$nul)){
  plot.tab <- d$all.em$nul[[i]]$plot.table
  N.em     <- nrow(plot.tab[plot.tab$relEM > 0.9,])
  N.am     <- nrow(plot.tab[plot.tab$relEM < 0.1,])
  relEM    <- mean(plot.tab$relEM)
  #Get boostrap CI on mean.
  boot.dat <-list()
  for(k in 1:1000){
    dat <- plot.tab[sample(nrow(plot.tab), nrow(plot.tab), replace = T),]
    boot.em    <- nrow(dat[dat$relEM > 0.9,])
    boot.am    <- nrow(dat[dat$relEM < 0.1,])
    boot.relEM <- mean(dat$relEM)
    boot.dat[[k]] <- c(boot.relEM,boot.em,boot.am)
  }
  boot.dat <- do.call(rbind, boot.dat)
  relEM.95 <- quantile(boot.dat[,1], probs = c(0.025, 0.975))
  em.95 <- quantile(boot.dat[,2], probs = c(0.025, 0.975))
  am.95 <- quantile(boot.dat[,3], probs = c(0.025, 0.975))
  
  #return output.
  em.ramp[[i]] <- c(relEM,N.em,N.am,relEM.95,em.95,am.95)
}
em.ramp <- data.frame(do.call(rbind, em.ramp))
colnames(em.ramp) <- c('relEM','N.em','N.am','relEM.lo95','relEM.hi95','em.lo95','em.hi95','am.lo95','am.hi95')

#starting from AM state.
am.ramp <- list()
for(i in 1:length(d$all.am$nul)){
  plot.tab <- d$all.am$nul[[i]]$plot.table
  N.em     <- nrow(plot.tab[plot.tab$relEM > 0.9,])
  N.am     <- nrow(plot.tab[plot.tab$relEM < 0.1,])
  relEM    <- mean(plot.tab$relEM)
  #Get boostrap CI on mean.
  boot.dat <-list()
  for(k in 1:1000){
    dat <- plot.tab[sample(nrow(plot.tab), nrow(plot.tab), replace = T),]
    boot.em    <- nrow(dat[dat$relEM > 0.9,])
    boot.am    <- nrow(dat[dat$relEM < 0.1,])
    boot.relEM <- mean(dat$relEM)
    boot.dat[[k]] <- c(boot.relEM,boot.em,boot.am)
  }
  boot.dat <- do.call(rbind, boot.dat)
  relEM.95 <- quantile(boot.dat[,1], probs = c(0.025, 0.975))
  em.95 <- quantile(boot.dat[,2], probs = c(0.025, 0.975))
  am.95 <- quantile(boot.dat[,3], probs = c(0.025, 0.975))
  #return output.
  am.ramp[[i]] <- c(relEM,N.em,N.am,relEM.95,em.95,am.95)
}
am.ramp <- data.frame(do.call(rbind, am.ramp))
colnames(am.ramp) <- c('relEM','N.em','N.am','relEM.lo95','relEM.hi95','em.lo95','em.hi95','am.lo95','am.hi95')

#Get EM dominated forest frequencies.
#renaming
#up   <- em.ramp[,c('N.em','em.lo95','em.hi95')]
#down <- am.ramp[,c('N.em','em.lo95','em.hi95')]
#colnames(  up) <- c('n','lo95','hi95')
#colnames(down) <- c('n','lo95','hi95')
#convert to frequencies.
#  up <-   up / n.tot
#down <- down / n.tot

#Use average relative abundance EM trees instead.
  up.null <- em.ramp[,c('relEM','relEM.lo95','relEM.hi95')]
down.null <- am.ramp[,c('relEM','relEM.lo95','relEM.hi95')]
colnames(  up.null) <- c('n','lo95','hi95')
colnames(down.null) <- c('n','lo95','hi95')

#Assign ndep.
up.null  $ndep <- c(1:n.lev)
down.null$ndep <- c(1:n.lev)

#Panel 1.  null simulation ramp up / ramp down.----
up <- up.null
down <- down.null
par(mar = c(0,3,0,0))
cols <- c('purple','green')
trans <- 0.3
#plot ramp up, transitioning away from EM dominated forests.
color <- cols[1]
max.y <- max(c(up.null$hi95,up.feed$hi95)) * 1.07
limy <- c(0,max.y)
#plot(up$n ~ up$ndep, bty = 'l', cex = 0, ylab = NA, xlab = NA, ylim = c(0.9*min(c(up$lo95, down$lo95)),max(up$n)*1.05))
plot(up$n ~ up$ndep, bty = 'l', cex = 0, ylab = NA, xlab = NA, ylim = limy)
lines(smooth.spline(up$n ~ up$ndep, spar = .1), lwd = 3, col = color)
polygon(c(up$ndep, rev(up$ndep)),c(up$hi95, rev(up$lo95)), col=adjustcolor(color, trans), lty=0)

#plot ramp down.
color <- cols[2]
lines(smooth.spline(down$n ~ down$ndep, spar = 0.1), lwd = 2, col = adjustcolor(color, 0.6))
polygon(c(down$ndep, rev(down$ndep)),c(down$hi95, rev(down$lo95)), col=adjustcolor(color, trans), lty=0)

#legend.
legend(x=2,y= 0.2,legend=c("initially EM","initially AM"),lty=1,lwd=2, col = cols, 
       box.lwd=0, bg="transparent",cex= 1)

#outer labels.
mtext('Relative Abundance Ectomycorrhizal Trees', side =2, cex=1.2, line = 2.5)
mtext('Null Simulation', side = 3, adj = 0.05, line = -1.5)
mtext('Nitrogen Deposition', side = 1, line = 2.5)
mtext('C', side = 1, line = -1.5, adj = 0.05, font = 2)
mtext(expression(paste("kg N ha"^"-1"," yr"^"-1",sep="")), side = 1, line = 3.5,  cex = 0.7)


#Panel 2. hysteresis feedback simulation ramp up / ramp down.----
  up <-   up.feed
down <- down.feed
par(mar = c(0,3,0,0))
cols <- c('purple','green')
trans <- 0.3
#plot ramp up, transitioning away from EM dominated forests.
color <- cols[1]
plot(up$n ~ up$ndep, bty = 'l', cex = 0, ylab = NA, xlab = NA, ylim = limy)
lines(smooth.spline(up$n ~ up$ndep, spar = .1), lwd = 3, col = color)
polygon(c(up$ndep, rev(up$ndep)),c(up$hi95, rev(up$lo95)), col=adjustcolor(color, trans), lty=0)

#plot ramp down.
color <- cols[2]
lines(smooth.spline(down$n ~ down$ndep, spar = 0.1), lwd = 2, col = adjustcolor(color, 0.6))
polygon(c(down$ndep, rev(down$ndep)),c(down$hi95, rev(down$lo95)), col=adjustcolor(color, trans), lty=0)

#outer labels.
#mtext('Number EM Dominated Forests', side =2, cex=1.2, line = 2.5)
mtext('Feedback Simulation', side = 3, adj = 0.05, line = -1.5)
mtext('Nitrogen Deposition', side = 1, line = 2.5)
mtext('C', side = 1, line = -1.5, adj = 0.05, font = 2)
mtext(expression(paste("kg N ha"^"-1"," yr"^"-1",sep="")), side = 1, line = 3.5,  cex = 0.7)

#end plot.----
dev.off()

