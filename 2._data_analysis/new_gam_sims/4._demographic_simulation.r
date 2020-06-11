#Ttesting demographic simulation model. 
rm(list=ls())
source('paths.r')
library(doParallel)
library(randomForest)
library(mgcv)
source('project_functions/gam.int_forest.sim.R')
source('project_functions/tic_toc.r')

#load results.
d <- readRDS(demographic_fits_gam_separate.path)
p2 <- readRDS(Product_2.subset.path)
p1 <- readRDS(Product_1.path)
p1 <- p1[p1$PLT_CN %in% p2$PLT_CN,]

#Just run the function.----
tic()
null <- forest.sim(g.mod.am = d$n.feedback$G.mod.am, g.mod.em = d$n.feedback$G.mod.em,
                   r.mod.am = d$n.feedback$R.mod.am, r.mod.em = d$n.feedback$R.mod.em,
                   m.mod.am = d$n.feedback$M.mod.am, m.mod.em = d$n.feedback$M.mod.em,
                   myco.split = 'between_plot',
                   env.cov = d$all.cov,
                   n.cores = 8,
                   n.plots = 500, n.step = 40)
cat('Null simulation complete.\n')
toc()
tic()
feed <- forest.sim(g.mod.am = d$y.feedback$G.mod.am, g.mod.em = d$y.feedback$G.mod.em,
                   r.mod.am = d$y.feedback$R.mod.am, r.mod.em = d$y.feedback$R.mod.em,
                   m.mod.am = d$y.feedback$M.mod.am, m.mod.em = d$y.feedback$M.mod.em,
                   myco.split = 'between_plot',
                   env.cov = d$all.cov,
                   n.cores = 8,
                   n.plots = 500, n.step = 40)
cat('Feedback simulation complete.\n')
toc()



#Diagnostics.----
par(mfrow = c(3,2), mar = c(2,2,0,0))
#EM-AM distribution.
hist(null$plot.table$relEM, breaks = 10, main = 'Relative Abundance EM Trees - NULL', cex = 0.8, ylim = c(0, 250))
hist(feed$plot.table$relEM, breaks = 10, main = 'Relative Abundance EM Trees - FEED', cex = 0.8, ylim = c(0, 250))
#self-thinning.
plot(log(feed$plot.table$BASAL.plot / feed$plot.table$stem.density) ~ log(feed$plot.table$stem.density), bty = 'l', main = 'self-thinning')
#Basal area and stem density distirbution.
hist(feed$plot.table$BASAL.plot, main = 'Basal Area')
hist(feed$plot.table$stem.density, main = 'stem density')

#Declaring environmental variables within the forest.sim_rf function for testing.-----
g.mod.am = d$n.feedback$G.mod.am
g.mod.em = d$n.feedback$G.mod.em
r.mod.am = d$n.feedback$R.mod.am
r.mod.em = d$n.feedback$R.mod.em
m.mod.am = d$n.feedback$M.mod.am
m.mod.em = d$n.feedback$M.mod.em
env.cov  <- d$all.cov
initial_density = 20
n.plots = 100
n.step = 20
disturb_rate = 0.018
step.switch = NA
switch.lev = NA #if changing N level mid run.
n.cores = detectCores()
silent = F
myco.split = 'between_plot'
