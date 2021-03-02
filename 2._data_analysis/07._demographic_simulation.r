#Ttesting demographic simulation model. 
rm(list=ls())
source('paths.r')
library(doParallel)
library(randomForest)
library(mgcv)
source('project_functions/gam.int_forest.sim.r')
source('project_functions/tic_toc.r')

#set output path.----
output.path <- null_vs_feedback_simulation_output_RE.path
  
#load gam model results.----
d <- readRDS(demographic_fits_gam_separate_plus_re_county.path)

#Just run the function.----
tic()
null <- forest.sim(g.mod.am = d$n.feedback$G.mod.am, g.mod.em = d$n.feedback$G.mod.em,
                   r.mod.am = d$n.feedback$R.mod.am, r.mod.em = d$n.feedback$R.mod.em,
                   m.mod.am = d$n.feedback$M.mod.am, m.mod.em = d$n.feedback$M.mod.em,
                   myco.split = 'between_plot',
                   env.cov = d$all.cov,
                   n.cores = 28,
                   n.plots = 1000, n.step = 40)
cat('Null simulation complete.\n')
toc()
tic()
feed <- forest.sim(g.mod.am = d$y.feedback$G.mod.am, g.mod.em = d$y.feedback$G.mod.em,
                   r.mod.am = d$y.feedback$R.mod.am, r.mod.em = d$y.feedback$R.mod.em,
                   m.mod.am = d$y.feedback$M.mod.am, m.mod.em = d$y.feedback$M.mod.em,
                   myco.split = 'between_plot',
                   env.cov = d$all.cov,
                   n.cores = 28,
                   n.plots = 1000, n.step = 40)
cat('Feedback simulation complete.\n')
toc()

#save output.----
out <- list(null,feed)
names(out) <- c('n.feedback','y.feedback')
saveRDS(out, output.path)


#Diagnostics.----
plotting = F
if(plotting == T){
  par(mfrow = c(3,2), mar = c(2,2,0,0))
  #EM-AM distribution.
  hist(null$plot.table$relEM, breaks = 10, main = 'Relative Abundance EM Trees - NULL', cex = 0.8, ylim = c(0, 250))
  hist(feed$plot.table$relEM, breaks = 10, main = 'Relative Abundance EM Trees - FEED', cex = 0.8, ylim = c(0, 250))
  #self-thinning.
  plot(log(feed$plot.table$BASAL.plot / feed$plot.table$stem.density) ~ log(feed$plot.table$stem.density), bty = 'l', main = 'self-thinning')
  #Basal area and stem density distirbution.
  hist(feed$plot.table$BASAL.plot, main = 'Basal Area')
  hist(feed$plot.table$stem.density, main = 'stem density')
}
