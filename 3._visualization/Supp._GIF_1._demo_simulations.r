#Create gifs of null and feedback simulations.
rm(list=ls())
library(magick)
library(purrr)
source('paths.r')
source('project_functions/make_gif.r')

#set output path.----
output.path <- Supp._gif_1.path

#Generate file structure for outputs.----
#set directory paths.
gif.dir <- 'figures/gif/'
nul.gif.dir <- paste0(gif.dir,'nul/')
alt.gif.dir <- paste0(gif.dir,'alt/')
all.gif.dir <- paste0(gif.dir,'all/')
#make directories.
system(paste0('mkdir -p ',nul.gif.dir))
system(paste0('mkdir -p ',alt.gif.dir))
system(paste0('mkdir -p ',all.gif.dir))

#Load null vs. feedback simulation output.----
d <- readRDS(null_vs_feedback_simulation_output.path)

#Null gif figures.----
z <- d$n.feedback$super.table             #grab null simulation supertable.
year <- seq(0, (length(z) - 1)*5, by = 5) #get years as vector.
gif.col <- 'light green'                  #pick histogram bar color.
limx <- c(0,1)                            #xlim always 0,1.
#limy <- c(0, nrow(z[[1]])/2)              #ylim alays 0, number of plots / 2

#plot loop.
for(i in 1:length(year)){
  #Generate filepath to save to.
  if(nchar(i) == 1){num <- paste0('0',i)}
  if(nchar(i)  > 1){num <- i}
  out.path <- paste0(nul.gif.dir,'nul_',num,'.png')
  
  #save line.
  png(out.path, width = 5, height = 5, units = 'in', res = 300)
  
  #set margins.
  par(mar = c(4,4,0,0), oma = c(0,0,0,0))
  
  #grab data for time increment.
  x <- z[[i]]$relEM
  
  #Drop plot.
  hist(x, main = NULL, xlab = NA, ylab = NA, col = gif.col, lty = 'blank', xlim = limx)
  #plot labels.
  mtext('Number of Forests'                       , side = 2, line = 2.4, cex = 1.2)
  mtext('Relative Abundance Ectomycorrhizal Trees', side = 1, line = 2.4, cex = 1.2)
  msg <- paste0('Year ',year[i])
  mtext(msg, side = 3, line = -1.5, adj = 0.05)
  
  #end this particular plot.
  dev.off()
}

#Feedback gif figures.----
z <- d$y.feedback$super.table             #grab null simulation supertable.
year <- seq(0, (length(z) - 1)*5, by = 5) #get years as vector.
gif.col <- 'light green'                  #pick histogram bar color.
limx <- c(0,1)                            #xlim always 0,1.
#limy <- c(0, nrow(z[[1]])/2)              #ylim alays 0, number of plots / 2

#plot loop.
for(i in 1:length(year)){
  #Generate filepath to save to.
  if(nchar(i) == 1){num <- paste0('0',i)}
  if(nchar(i)  > 1){num <- i}
  out.path <- paste0(alt.gif.dir,'alt_',num,'.png')
  
  #save line.
  png(out.path, width = 5, height = 5, units = 'in', res = 300)
  
  #set margins.
  par(mar = c(4,4,0,0), oma = c(0,0,0,0))
  
  #grab data for time increment.
  x <- z[[i]]$relEM
  
  #Drop plot.
  hist(x, main = NULL, xlab = NA, ylab = NA, col = gif.col, lty = 'blank', xlim = limx)
  #plot labels.
  mtext('Number of Forests'                       , side = 2, line = 2.4, cex = 1.2)
  mtext('Relative Abundance Ectomycorrhizal Trees', side = 1, line = 2.4, cex = 1.2)
  msg <- paste0('Year ',year[i])
  mtext(msg, side = 3, line = -1.5, adj = 0.05)
  
  #end this particular plot.
  dev.off()
}

#Null and feedback plots together gif figures.----
z1 <- d$n.feedback$super.table             #grab null simulation supertable.
z2 <- d$y.feedback$super.table             #grab feedback simulation supertable.
year <- seq(0, (length(z1) - 1)*5, by = 5) #get years as vector.
gif.col <- 'light green'                   #pick histogram bar color.
limx <- c(0,1)                             #xlim always 0,1.

#plot loop.
for(i in 1:length(year)){
  #Generate filepath to save to.
  if(nchar(i) == 1){num <- paste0('0',i)}
  if(nchar(i)  > 1){num <- i}
  out.path <- paste0(all.gif.dir,'all_',num,'.png')
  
  #save line.
  png(out.path, width = 10, height = 5, units = 'in', res = 300)
  
  #set margins and number of panels.
  par(mar = c(2,2,0,0), oma = c(2,2,2,0), mfrow = c(1,2))
  
  #grab data for time increment.
  x1 <- z1[[i]]$relEM
  x2 <- z2[[i]]$relEM
  
  #set y-limit.
  n.breaks = 10
  x1.ref   <- cut(x1, n.breaks)
  count.x1 <- table(x1.ref)
  x2.ref   <- cut(x2, n.breaks)
  count.x2 <- table(x2.ref)
  check <- max(c(count.x1, count.x2))
  #This doesn't work very well. Setting increments of ylimit by "hand".
  if(check > 400)               {limy = c(0,500)}
  if(check < 400 & check >= 300){limy = c(0,400)}
  if(check < 300)               {limy = c(0,300)}
  
  #Drop plots and inner labels.
  #null histogram.
  hist(x1, main = NULL, xlab = NA, ylab = NA, col = gif.col, lty = 'blank', 
       xlim = limx, ylim = limy, breaks = n.breaks)
  mtext('Simulation without \ncon-mycorrhizal feedbacks', side = 3, line = -3, adj = 0.05)
  #feedback histogram.
  hist(x2, main = NULL, xlab = NA, ylab = NA, col = gif.col, lty = 'blank', 
       xlim = limx, ylim = limy, breaks = n.breaks)
  mtext('Simulation with \ncon-mycorrhizal feedbacks'   , side = 3, line = -3, adj = 0.05)
  #outer plot labels.
  mtext('Number of Forests'                       , side = 2, line = 0.8, cex = 1.2, outer = T)
  mtext('Relative Abundance Ectomycorrhizal Trees', side = 1, line = 0.8, cex = 1.2, outer = T)
  msg <- paste0('Year ',year[i])
  mtext(msg, side = 3, line = 0.5, adj = 0.05, cex = 1.2, outer = T)
  
  #end this particular plot.
  dev.off()
}

#Make your gif!----
make_gif(all.gif.dir, output.path)

#clean up files.----
cmd <- paste0('rm -r ',gif.dir)
system(cmd)

#end script.----
