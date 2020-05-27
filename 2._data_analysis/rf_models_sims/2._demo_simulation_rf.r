#Ttesting demographic simulation model. 
rm(list=ls())
source('paths.r')
library(doParallel)
library(randomForest)
source('project_functions/rf_forest.sim.R')
source('project_functions/tic_toc.r')

#load results.
d <- readRDS(rf_demographic_fits.path)

#Just run the function.----
tic()
test <- forest.sim(g.mod.am = d$models$grow.mod.am, g.mod.em = d$models$grow.mod.em,
                   r.mod.am = d$models$recr.mod.am, r.mod.em = d$models$recr.mod.em,
                   m.mod.am = d$models$mort.mod.am, m.mod.em = d$models$mort.mod.em,
                   myco.split = 'between_plot',
                   env.cov = d$env.cov,
                   n.cores = 8,
                   n.plots = 200)
toc()

#With intearctive model.----
d <- readRDS(rf_demographic_fits_interactive.path)
tic()
test <- forest.sim(g.mod.am = d$models$grow.mod, g.mod.em = d$models$grow.mod,
                   r.mod.am = d$models$recr.mod.am, r.mod.em = d$models$recr.mod.em,
                   m.mod.am = d$models$mort.mod, m.mod.em = d$models$mort.mod,
                   myco.split = 'between_plot',
                   env.cov = d$env.cov,
                   n.cores = 8,
                   n.plots = 200)
toc()


#Diagnostics.----
par(mfrow = c(3,2))
#EM-AM distribution.
hist(d$data$recr.dat.all$relEM, breaks = 10, main = 'Relative Abundance EM Trees - OBS',  cex = 0.8)
hist(test$plot.table$relEM, breaks = 10, main = 'Relative Abundance EM Trees - SIM', cex = 0.8)
#self-thinning.
plot(log(test$plot.table$BASAL.plot / test$plot.table$stem.density) ~ log(test$plot.table$stem.density), bty = 'l', main = 'self-thinning')
#Basal area and stem density distirbution.
hist(test$plot.table$BASAL.plot)
hist(test$plot.table$stem.density)

#Declaring environmental variables within the forest.sim_rf function for testing.-----
g.mod.am <- d$models$grow.mod.am
g.mod.em <- d$models$grow.mod.em
r.mod.am <- d$models$recr.mod.am
r.mod.em <- d$models$recr.mod.em
m.mod.am <- d$models$mort.mod.am
m.mod.em <- d$models$mort.mod.em
env.cov  <- d$env.cov
initial_density = 20
n.plots = 100
n.step = 20
disturb_rate = 0.018
step.switch = NA
switch.lev = NA #if changing N level mid run.
n.cores = detectCores()
silent = F
myco.split = 'between_plot'