rm(list=ls())
source('paths.r')
library(usmap)
library(ggplot2)
library(ggpubr)
library(data.table)
source('project_functions/crib_fun.r')

#set output path.----
output.path <- Fig_2.path
  
#load data.----
#plots.
d <- readRDS(Product_2.subset.path)
d <- d[,.(PLT_CN,LAT,LON,STATECD,relEM)]
setkey(d,'PLT_CN')
d <- unique(d)
#states.
states <- read.csv('required_products_utilities/FIA_state_codes_regions.csv')
states.east <- states[states$east_plus == 1,]$state_name
states.east <- as.character(states.east)

#generate histogram data.----
d2 <- readRDS(Product_2.subset.path)
d  <- readRDS(Product_1.path)
d  <- d[d$PLT_CN %in% d2$PLT_CN,]
#Subset to site level.
setkey(d,'PLT_CN')
d <- d[complete.cases(d),]

#fit a model.
#get predictor and y data.
d$relEM <- crib_fun(d$relEM)
mod <- betareg::betareg(relEM ~ mat + map + ndep, data = d)

#Detrend data.
#predicted values
parms <- coef(mod)[1:length(coef(mod)) -1]
obs <- logit(d$relEM)
x <- d[,.(mat,map,ndep)]
pred <- predict(mod, newdata = x)
mean.pred <- colMeans(x, na.rm = T)
mean.pred <- predict(mod, newdata = data.frame(t(mean.pred)))
corrected <- obs - boot::logit(pred) + boot::logit(mean.pred)[1]
d$corrected <- inv.logit(corrected)

#code to plot the eastern US and sites.----
lon <- d$longitude
lat <- d$latitude
#generate colors.
rbPal <- colorRampPalette(c('green','pink','purple'))
d$col <- rbPal(10)[as.numeric(cut(d$corrected,breaks = 10))]
d$col <- rbPal(30)[as.numeric(cut(d$corrected,breaks = 30))]
world <- map_data('usa')
map <- ggplot() + geom_cartogram(data = world, map = world, 
                                 aes(x=long, y = lat, group = group, map_id=region))
map <- map + coord_proj("+proj=cea +lat_ts=37.5",ylim = c(25,50), xlim = c(-95,-66.5))
map <- map + geom_point(aes(x = lon, y = lat), color = d$col, size = .2)
map <- map + theme(axis.line=element_blank(),
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
#Legend? Not working.
#map <- map + scale_colour_gradient2('%EM',low='green', mid = 'pink', high='purple', midpoint = 0.5)

#Generate histogram plot.----
#get colors for histogram plot.
d.col <- d[order(d$corrected),]
hist.cols <- unique(d.col$col)
#drop hist.plot.
hist.plot <-
ggplot(d, aes(x=corrected)) + geom_histogram(bins = 30, fill = hist.cols) + theme_classic() + 
  xlab('Relative Abundance Ectomycorrhizal Trees') +
  ylab('Number of Forests')# + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))

#setup 2-panel plot.----
#png save line.
png(output.path, width = 8, height = 4, units = 'in', res = 300)
#Drop plots.
ggarrange(map, hist.plot,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
#end plot and save.
dev.off()