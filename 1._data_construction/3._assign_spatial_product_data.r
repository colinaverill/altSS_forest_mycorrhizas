#Assign mat, map and n-deposition from spatial products.
rm(list=ls())
library(data.table)
library(raster)
library(rgdal)
library(sp)
#detach("package:runjags", unload=TRUE) #make sure runjags is not loaded.
source('project_functions/extract_ndep.r')
source('project_functions/worldclim2_grab.r')
source('paths.r')

#Load data.----
p1 <- data.table(readRDS(Product_1.path))
p2 <- data.table(readRDS(Product_2.path))

#remove spatial columns if already present.----
to.drop <- c('n.dep','wet.dep','dry.dep','mat','map','mat_CV','map_CV','mat_sd','map_sd','mdr')
p1 <- p1[,c(to.drop) := NULL]
p2 <- p2[,c(to.drop) := NULL]

#pull N deposition data. ignore warnings.-----
p1 <- cbind(p1,extract_ndep(p1$longitude, p1$latitude))
setnames(p1, 'n.dep','ndep')
p2 <- merge(p2,p1[,c('PLT_CN','ndep','dry.dep','wet.dep')])

#pull worldclim2 climate data.----
p1 <- cbind(p1,worldclim2_grab(p1$latitude,p1$longitude))
p2 <- merge(p2, p1[,c('PLT_CN','mat','map')])

#save output.----
saveRDS(p1, Product_1.path, version = 2)
saveRDS(p2, Product_2.path, version = 2)
