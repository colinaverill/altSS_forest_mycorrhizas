rm(list=ls())
source('paths.r')

#load data.----
d <- readRDS(initial_condition_hysteresis_simulation.path)

#Gather up initially AM and Em time series across Ndep levels.----
all.em <- list()
all.am <- list()
for(k in 1:length(d$all.em$alt.GRM)){
  am <- d$all.am$alt.GRM[[k]]$super.table
  em <- d$all.em$alt.GRM[[k]]$super.table
  em.series <- list()
  am.series <- list()
  for(i in 1:length(am)){
    plot.tab.am <- am[[i]]
    plot.tab.em <- em[[i]]
    am.N.em <- nrow(plot.tab.am[plot.tab.am$relEM > 0.9,])
    am.N.am <- nrow(plot.tab.am[plot.tab.am$relEM < 0.1,])
    em.N.em <- nrow(plot.tab.em[plot.tab.em$relEM > 0.9,])
    em.N.am <- nrow(plot.tab.em[plot.tab.em$relEM < 0.1,])
    am.series[[i]] <- c(am.N.em,am.N.am)
    em.series[[i]] <- c(em.N.em,em.N.am)
  }
  am.series <- data.frame(do.call(rbind, am.series))
  em.series <- data.frame(do.call(rbind, em.series))
  colnames(am.series) <- c('N.em','N.am')
  colnames(em.series) <- c('N.em','N.am')
  all.am[[k]] <- am.series
  all.em[[k]] <- em.series
}
names(all.am) <- names(d$all.am$alt.GRM)
names(all.em) <- names(d$all.em$alt.GRM)

#Plot AM and EM time series across Ndep levels.----
#plotting abundance of EM dominated plots.
par(mfrow = c(2,4), mar = c(2,2,0,0), oma = c(2,2,2,2))
#initially AM plots.
plot(all.am$l1 $N.em, bty = 'l')
plot(all.am$l5 $N.em, bty = 'l')
plot(all.am$l10$N.em, bty = 'l')
plot(all.am$l14$N.em, bty = 'l')
#initially EM plots.
plot(all.em$l1 $N.em, bty = 'l')
plot(all.em$l5 $N.em, bty = 'l')
plot(all.em$l10$N.em, bty = 'l')
plot(all.em$l14$N.em, bty = 'l')

#Plotting abundance of AM dominated plots.
par(mfrow = c(2,4), mar = c(2,2,0,0), oma = c(2,2,2,2))
#initially AM plots.
plot(all.am$l1 $N.am, bty = 'l')
plot(all.am$l5 $N.am, bty = 'l')
plot(all.am$l10$N.am, bty = 'l')
plot(all.am$l14$N.am, bty = 'l')
#initially EM plots.
plot(all.em$l1 $N.am, bty = 'l')
plot(all.em$l5 $N.am, bty = 'l')
plot(all.em$l10$N.am, bty = 'l')
plot(all.em$l14$N.am, bty = 'l')



#get number of EM and AM plots in feedback simulations (>90%) and bootstrap CI.----
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

