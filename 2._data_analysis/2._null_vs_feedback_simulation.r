#running null and feedback simulations drawing sampling plot-level environmental conditions from dataset.
#Takes ~14 minutes to run w/ 2 cores on pecan2.
rm(list=ls())
source('paths.r')
source('project_functions/forest.sim.r')
source('project_functions/tic_toc.r')
source('project_functions/makeitwork.r')
library(mgcv)
library(doParallel)
library(data.table)

#set output path.----
output.path <- null_vs_feedback_simulation_output.path

#RNG seed for reproducibility.----
#7 - close - even AM bars.
#42069 - no.
#42 - just slightly higher on AM extreme. best so far.
#1- next test.
set.seed(42)

#load models and environmental covariates.----
fits <- readRDS(demographic_fits.path) #trying with new density dependence for recruitment.
env.cov <- data.frame(t(fits$env.cov))
env.cov <- fits$all.cov
N.PLOTS <- 1000 #Must be even!

#register parallel environment.----
n.cores <- detectCores() - 1 #minus 1 so your computer keeps running.

#Run drawwing from distribution of plot level environmental conditions.----
#Null model.
tic()
cat('running null simulations...\n')
nul <-     forest.sim(g.mod    = fits$n.feedback$G.mod, 
                      m.mod    = fits$n.feedback$M.mod,
                      r.mod.am = fits$n.feedback$R.mod.am, 
                      r.mod.em = fits$n.feedback$R.mod.em,
                      env.cov = env.cov, 
                      myco.split = 'between_plot',
                      disturb_rate = 0.0476/2,
                      n.plots = N.PLOTS,
                      n.cores = n.cores)
cat('null simulations complete. \n')
toc()
#Feedback model.
tic()
cat('running feedback simulations...\n')
fed <-     forest.sim(g.mod    = fits$y.feedback$G.mod, 
                      m.mod    = fits$y.feedback$M.mod,
                      r.mod.am = fits$y.feedback$R.mod.am, 
                      r.mod.em = fits$y.feedback$R.mod.em,
                      env.cov = env.cov, 
                      myco.split = 'between_plot',
                      disturb_rate = 0.0476/2,
                      n.plots = N.PLOTS,
                      n.cores = n.cores)
cat('feedback simulations complete. \n')
toc()

#save output.----
out <- list(nul,fed)
names(out) <- c('n.feedback','y.feedback')
saveRDS(out, output.path, version = 2)

#visualize.----
visualize = F
if(visualize == T){
  par(mfrow = c(4,2))
  a <- out$n.feedback$plot.table
  b <- out$y.feedback$plot.table
  hist(a$relEM, xlim = c(0,1))
  hist(b$relEM, xlim = c(0,1))
  #basal area.
  hist(a$BASAL.plot)
  hist(b$BASAL.plot)
  #stem density.
  hist(a$stem.density)
  hist(b$stem.density)
  #self thinning. Not sure if I am plotting this correctly.
  plot(log(a$BASAL.plot/a$stem.density) ~ log(a$stem.density), bty = 'l')
  plot(log(b$BASAL.plot/b$stem.density) ~ log(b$stem.density), bty = 'l')
}