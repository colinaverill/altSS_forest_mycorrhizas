#running null and feedback simulations drawing sampling plot-level environmental conditions from dataset.
#~9 minutes using 4 cores on my 2014 macbook pro.
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
set.seed(42069) #looks good!

#load models and environmental covariates.----
fits <- readRDS(demographic_fits.path)
env.cov <- data.frame(t(fits$env.cov))
env.cov <- fits$all.cov
N.PLOTS <- 1000 #Must be even!

#register parallel environment.----
n.cores <- detectCores()

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
                      #disturb_rate = 0.0476/2,
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
visualize = T
if(visualize == T){
  par(mfrow = c(4,2))
  a <- out$n.feedback$plot.table
  b <- out$y.feedback$plot.table
  hist(a$relEM, xlim = c(0,1), main = 'null model')
  hist(b$relEM, xlim = c(0,1), main = 'feedback model')
  #basal area goes too high for deedback, likely driven by the ridiculous stem densities.
  hist(a$BASAL.plot)
  hist(b$BASAL.plot)
  #stem density - stem density goes too high for feedback. 
  hist(a$stem.density)
  hist(b$stem.density)
  #self thinning. Not sure if I am plotting this correctly.
  plot(log(a$BASAL.plot/a$stem.density) ~ log(a$stem.density), bty = 'l')
  plot(log(b$BASAL.plot/b$stem.density) ~ log(b$stem.density), bty = 'l')
}
