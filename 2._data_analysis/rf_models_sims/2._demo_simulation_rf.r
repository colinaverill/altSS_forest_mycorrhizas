#Ttesting demographic simulation model. 
rm(list=ls())
source('paths.r')
library(doParallel)
library(randomForest)
source('project_functions/rf_forest.sim.R')
source('project_functions/tic_toc.r')

#set output path.----
output.path <- RF_null_vs_feedback_simulation_output.path
  
#load model fits.----
d <- readRDS(rf_demographic_fits.path)
env.cov <- d$env.cov
env.cov$BASAL.am <- NULL
env.cov$BASAL.em <- NULL

#Running simulations.----
#models w/o feedbacks.
tic()
cat('Fitting null models...\n')
null <- forest.sim(g.mod.am = d$null.models$grow.mod.am, g.mod.em = d$null.models$grow.mod.em,
                   r.mod.am = d$null.models$recr.mod.am, r.mod.em = d$null.models$recr.mod.em,
                   m.mod.am = d$null.models$mort.mod.am, m.mod.em = d$null.models$mort.mod.em,
                   myco.split = 'between_plot',
                   env.cov = env.cov,
                   n.cores = 28,
                   n.plots = 500, n.step = 40)
cat('Null models fit.\n')
toc()

#models w/ feedbacks.
tic()
cat('Fitting feedback models...\n')
feed <- forest.sim(g.mod.am = d$feedback.models$grow.mod.am, g.mod.em = d$feedback.models$grow.mod.em,
                   r.mod.am = d$feedback.models$recr.mod.am, r.mod.em = d$feedback.models$recr.mod.em,
                   m.mod.am = d$feedback.models$mort.mod.am, m.mod.em = d$feedback.models$mort.mod.em,
                   myco.split = 'between_plot',
                   env.cov = d$env.cov,
                   n.cores = 28,
                   n.plots = 500, n.step = 40)
cat('Feedback mdeols fit.\n')
toc()

#wrap and save.----
output <- list(null,feed)
names(output) <- c('null.sim','feed.sim')
saveRDS(output, output.path)

plotting = F
if(plotting == T){
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
}
