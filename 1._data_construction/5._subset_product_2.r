#Subset p1 and p2 to get a random stratified sampling for fitting models quicky.
#CURRENTLY SUBSETS TO NORTHEAST US.
rm(list=ls())
source('paths.r')
library(ggplot2)
library(ggalt)
library(data.table)

#set output path.----
output.path <- Product_2.subset.path

#load data.----
p2 <- readRDS(Product_2.path)
states <- read.csv('required_products_utilities/FIA_state_codes_regions.csv')
states.east <- states[states$east_plus == 1,]
states$northeast <- ifelse(states$state_name %in% c('Connecticut','New York','Massachusetts','Rhode Island','New Hampshire','Vermont','Maine'),1,0)

#Some final data cleaning (should be moved?).-----
#Would be nice to move these to full data filtering script (#2).
p2 <- p2[p2$REMPER >=4.9 & p2$REMPER <= 5.1,]
p2 <- p2[p2$n.trees >  5,]
#Complete cases of environmental covariates.
env.drop <- p2[,c('PLT_CN','ndep','mat','map','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')]
env.drop <- env.drop[!complete.cases(env.drop),]
env.drop <- unique(env.drop$PLT_CN)
p2 <- p2[!(p2$PLT_CN %in% env.drop),]

#Grab a plot table with PLT_CN, lat-lon and STATECD.----
d <- p2[,.(PLT_CN,LAT,LON,STATECD,relEM,REMPER)]
setkey(d, 'PLT_CN')
d <- unique(d)

#Subset.----
#subset to eastern US.
d <- d[d$STATECD %in% states.east$STATECD,]
#NORTHEAST FOR TESTING.
d <- d[d$STATECD %in% states[states$northeast == 1,]$STATECD,]
#set.seed(42069)
#n.plots <- 4000 #set number of plots to subsample.
#d <- d[sample(nrow(d), n.plots),]
p2 <- p2[PLT_CN %in% d$PLT_CN,]



#plot subset to make sure its representative.----
plot = T
if(plot == T){
  library(rworldmap)
  lon <- d$LON
  lat <- d$LAT
  newmap <- getMap(resolution = 'low')
  plot(newmap, xlim = c(-96, -69), ylim = c(25, 50))
  points(lon, lat, col = 'purple', cex = 0.2)
}

#Save output.----
saveRDS(p2, output.path)