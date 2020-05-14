#Assign mat, map and n-deposition from spatial products.
rm(list=ls())
library(data.table)
library(raster)
library(rgdal)
library(sp)
source('project_functions/extract_ndep.r')
source('project_functions/worldclim2_grab.r')
source('paths.r')

#Load data.----
p1 <- data.table(readRDS(Product_1.path))
p2 <- data.table(readRDS(Product_2.path))
composite <- data.table(read.csv(data_from_composite.path))

#remove spatial columns if already present.----
#to.drop <- c("n.dep","ndep","wet.dep","dry.dep",
#             "mat","t_seas","t_range","t_wet_q","t_warm_q","map","p_wet_q","p_warm_q",
#             "decomp_wet_q","decomp_warm_q")
#p1 <- p1[,c(to.drop) := NULL]
#p2 <- p2[,c(to.drop) := NULL]

#pull N deposition data. ignore warnings.-----
p1 <- cbind(p1,extract_ndep(p1$longitude, p1$latitude))
setnames(p1, 'n.dep','ndep')
p2 <- merge(p2,p1[,c('PLT_CN','ndep','dry.dep','wet.dep')])

#pull worldclim2 climate data.----
#wc2 <- worldclim2_grab(p1$latitude, p1$longitude)
#p1 <- cbind(p1, wc2)
#col.names <- c('PLT_CN', colnames(wc2))
#p1 <- as.data.frame(p1)
#p2 <- merge(p2, p1[,col.names])

#calculate decomposition in warmest and wettest quarters, following Steidinger et al. 2019, Nature.----
#calculation from Yasso7 model.
#decomp_calc <- function(temp, prec){
#  k <- exp(0.095*temp - 0.00014 * temp^2) * (1 - exp(-1.21 * prec))
#  return(k)
#}
#p1$decomp_wet_q  <- decomp_calc(p1$t_wet_q , p1$p_wet_q )
#p1$decomp_warm_q <- decomp_calc(p1$t_warm_q, p1$p_warm_q)
#p2 <- merge(p2, p1[,c('PLT_CN','decomp_wet_q','decomp_warm_q')])

#save output.----
saveRDS(p1, Product_1.path, version = 2)
saveRDS(p2, Product_2.path, version = 2)
