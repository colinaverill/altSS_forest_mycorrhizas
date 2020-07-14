#Ttesting demographic simulation model. 
rm(list=ls())
source('paths.r')
library(doParallel)
library(randomForest)
library(mgcv)
source('project_functions/gam.int_forest.sim.r')
source('project_functions/tic_toc.r')

#set output path.----
#output.path <- null_vs_feedback_simulation_output_RE.path

#load gam model results.----
all <- readRDS(demographic_fits_gam_ecoregion.path)

#Run null ad feedback simulations.-----
#setup output list.
region.simulation.out <- list()

#begin simulation loop.
for(i in 1:length(all)){
  #grab data subset and announce.
  d <- all[[i]]$county.re
  lab <- names(all)[i]
  cov <- d$all.cov
  cov <- cov[,!(names(cov) %in% c('BASAL.plot','BASAL.am','BASAL.em','stem.density'))]
  cat(paste0('Simulating '),lab,'...');tic()
  
  #Null simulation.
  #null <- forest.sim(g.mod.am = d$n.feedback$G.mod.am, g.mod.em = d$n.feedback$G.mod.em,
  #                   r.mod.am = d$n.feedback$R.mod.am, r.mod.em = d$n.feedback$R.mod.em,
  #                   m.mod.am = d$n.feedback$M.mod.am, m.mod.em = d$n.feedback$M.mod.em,
  #                   myco.split = 'between_plot',
  #                   env.cov = d$all.cov,
  #                   n.cores = 8,
  #                   n.plots = 1000, n.step = 40)
  #cat('Null simulation complete.\n');toc()
  
  #Feedback simulation.
  feed <- forest.sim(g.mod.am = d$G.mod.am, g.mod.em = d$G.mod.em,
                     r.mod.am = d$R.mod.am, r.mod.em = d$R.mod.em,
                     m.mod.am = d$M.mod.am, m.mod.em = d$M.mod.em,
                     myco.split = 'between_plot',
                     env.cov = cov,
                     n.cores = 28,
                     n.plots = 1000, n.step = 40)
  #cat('Feedback simulation complete.\n')
  
  #wrap output and save in output list.
  #out <- list(null,feed)
  #names(out) <- c('n.feedback','y.feedback')
  region.simulation.out[[i]] <- feed
  
  #report
  cat(paste0(lab,' analyses completed.\n'));toc()
  toc()
  
}

#save output.----
saveRDS(out, output.path)
