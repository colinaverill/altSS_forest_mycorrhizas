rm(list=ls())
source('paths.r')
library(ggplot2)
library(ggpubr)
library(grid)
library(ggalt)

#set output path.----
output.path <- Fig_1.path

#load data.----
d.dat <- readRDS(bimodality_results.path)
d.map <- d.dat$data
d.map$relEM.detrend <- boot::inv.logit(d.map$relEM.detrend)

#code to plot the eastern US and sites.----
#specifiy lat/lon.
lon <- d.map$longitude
lat <- d.map$latitude

#generate colors.
rbPal <- colorRampPalette(c('green','pink','purple'))
d.map$col <- rbPal(30)[as.numeric(cut(d.map$relEM,breaks = 30))]

#map code.
world <- map_data('usa')
map <- ggplot() + geom_cartogram(data = world, map = world,
                                aes(x=long, y = lat, group = group, map_id=region))
map <- map + coord_fixed(1.5, xlim=c(-95, -66.5)) #clip map to eastern US.
#map <- map + coord_proj("+proj=cea +lat_ts=37.5",ylim = c(25,50), xlim = c(-95,-66.5)) #crib map to eastern US.
map <- map + geom_point(aes(x = lon, y = lat), color = d.map$col, size = .2)              #drop points w/ particular colors.
map <- map + theme(axis.line=element_blank(),                                          #turn off a bunch of labels.
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   legend.position="none",
                   panel.border=element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank()
)

#Generate histogram plots.----
#get colors for histogram plot.
d.col <- d.map[order(d.map$relEM),]
hist.cols <- unique(d.col$col)
#drop hist.plot.
hist.plot.raw     <-
  ggplot(d.map, aes(x=relEM)) + geom_histogram(bins = 30, fill = hist.cols) + theme_classic() + 
  xlab('Relative Abundance Ectomycorrhizal Trees') +
  ylab('Number of Forests') + 
  ylim(0, 900) + 
  #add text.
  annotation_custom(grobTree(textGrob("Raw data", x=0.1,  y=0.95, hjust=0,
                          gp=gpar(fontface="italic"))))
#detrended histogram.
hist.plot.detrend <-
  ggplot(d.map, aes(x=relEM.detrend)) + geom_histogram(bins = 30, fill = hist.cols) + theme_classic() + 
  xlab('Relative Abundance Ectomycorrhizal Trees') +
  ylab('Number of Forests') + 
  ylim(0, 900) +
  #add text.
  annotation_custom(grobTree(textGrob("Data with environmental \nsignatures removed", x=0.1,  y=0.90, hjust=0,
                                      gp=gpar(fontface="italic"))))

#setup 3-panel plot.----
#pdf save line.
#png(output.path, width = 12, height = 4.5, units = 'in', res = 300)
#Drop plots.
#ggarrange(map, hist.plot.raw, hist.plot.detrend,
#          labels = c("A", "B","C"),
#          ncol = 3, nrow = 1)
#end plot and save.
#dev.off()

#png(output.path, width=8, height=7, units='in', res=300)
pdf(output.path, width=8, height=7)
ggarrange(map,
          ggarrange(hist.plot.raw, hist.plot.detrend, nrow = 2, labels = c('B','C')),
          ncol = 2,
          widths = c(1, 0.75),
          labels = 'A'
          )
dev.off()
