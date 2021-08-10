#visualizing null vs. feedback demographic simulations.
#5-panel to show recruitment and mortality positive feedbacks, demo simulation altnerative stable states, and hysteresis.
rm(list=ls())
source('paths.r')
library(mgcv)

#set output path.----
output.path <- Fig_3.path

#load demographic simulations.----
d <- readRDS(null_vs_feedback_simulation_output_RE.path)
y <- d$y.feedback$super.table[[41]] #200 years in.
n <- d$n.feedback$super.table[[41]]

#set colors, calculate y limit.----
#colors.
cols <- c('#00acd9','#cfe83c') #pick colors
cols <- c('light green','light green')
trans <- 0.3 #set transparency

#number of breaks.
n.breaks <- 20

#generate colors.
rbPal <- colorRampPalette(c('green','pink','purple'))
y$col <- rbPal(n.breaks)[as.numeric(cut(y$relEM,breaks = n.breaks))]
d.col <- y[order(y$relEM),]
hist.cols <- unique(d.col$col)


#y-limit.
y$cat <- cut(y$relEM, n.breaks)
count.y <- table(y$cat)
n$cat <- cut(n$relEM, n.breaks)
count.x <- table(n$cat)
limy <- c(0, max(c(count.x, count.y)))


#png save line.----
png(output.path, width = 11, height = 5, units = 'in', res = 300)
par(mfrow = c(1,2))

#plot simulation without feedbacks.----
tab <- n
par( mar = c(4,4,1,2))
#plot line.
hist(tab$relEM, breaks = n.breaks, 
     xlim = c(0,1), ylim = limy, 
     ylab = NA, xlab = NA, main = NA, 
     bty = 'l', col = hist.cols, lty = 'blank')
#label.
mtext(expression(paste("Relative Abundance Ectomycorrhizal Trees")), side = 1, line = 2.75, cex = 1)
mtext('Number of Forests', side = 2, line = 2.5, cex = 1)
msg1 <- expression(paste('Simulation ',bolditalic('without')))
msg2 <- 'con-mycorrhizal feedbacks'
mtext(msg1, side = 3, line = -3, adj = 0.40, cex = 1.2)
mtext(msg2, side = 3, line = -4, adj = 0.40, cex = 1.2)
mtext('A', side = 3, line = -3, adj = 0.05, cex = 1, font = 2)

#plot simulation without feedbacks.----
tab <- y
#plot line.
hist(tab$relEM, breaks = n.breaks, 
     xlim = c(0,1), ylim = limy, 
     ylab = NA, xlab = NA, main = NA, 
     bty = 'l', col = hist.cols, lty = 'blank')
#label.
mtext(expression(paste("Relative Abundance Ectomycorrhizal Trees")), side = 1, line = 2.75, cex = 1)
mtext('Number of Forests', side = 2, line = 2.5, cex = 1)
msg1 <- expression(paste('Simulation ',bolditalic('with')))
msg2 <- 'con-mycorrhizal feedbacks'
mtext(msg1, side = 3, line = -3, adj = 0.40, cex = 1.2)
mtext(msg2, side = 3, line = -4, adj = 0.40, cex = 1.2)
mtext('B', side = 3, line = -3, adj = 0.05, cex = 1, font = 2)

#end plot.----
dev.off()
