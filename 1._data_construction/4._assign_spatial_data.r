#PCA of environmental factors.
#Should re-run only using the plots you will analyze. May increase variance explained by first 10 eigenvectors.
#Will add in Ndep separately.
rm(list = ls())
source('paths.r')
source('project_functions/tic_toc.r')
source('project_functions/extract_ndep.r')
library(factoextra)
library(data.table)
library(raster)

#set output paths.----
#they are the same as the input paths. Colin, you should know better.

#load data.----
p1 <- readRDS(Product_1.path)
p2 <- readRDS(Product_2.path)
composite <- fread(data_from_composite.path)
#Subset to plots actually in data we're working with.
composite$PLT_CN <-as.numeric(composite$PLT_CN)
composite <- composite[composite$PLT_CN %in% p1$PLT_CN,]

#drop columns if you've already assigned and are reassigning for some reason.----
to.drop <- colnames(p1)[grep(c('mat|map|ndep|dry.dep|wet.dep|PC|ecoregion'),colnames(p1))]
p1 <- p1[,!(colnames(p1) %in% to.drop)]
to.drop <- colnames(p2)[grep(c('mat|map|ndep|dry.dep|wet.dep|PC|ecoregion'),colnames(p2))]
to.drop <- to.drop[!(to.drop %in% c('SPCD'))]
p2 <- as.data.frame(p2)
p2 <- p2[,!(colnames(p2) %in% to.drop)]

#Drop a few things from the composite linked to the soil feedbacks, or from maps derived from composite.----
composite <- composite[, grep("CHELSA" , colnames(composite)):=NULL]        #CHELSEA variables redundant w/ worldclim.
composite <- composite[, grep("Biomass", colnames(composite)):=NULL]        #tree and nematode biomass.
composite <- composite[, grep("Hansen", colnames(composite)):=NULL]         #dont care about year 2000 tree cover.
composite <- composite[, grep("LandCoverClass", colnames(composite)):=NULL] #all are forested.
composite <- composite[, grep("Mycorrhizal", colnames(composite)):=NULL]    #mycorrhizal info is ground sourced.
composite <- composite[, grep("Nematodes", colnames(composite)):=NULL]      #no nematodes.
composite <- composite[, grep("Resolve_Biome", colnames(composite)):=NULL]  #all are forests.
composite <- composite[, grep("Pixel_Area", colnames(composite)):=NULL]     #pixel area not important here.
composite <- composite[, grep("Soil_pH", colnames(composite)):=NULL]        #soil pH is a hypothesized feedback.
composite <- composite[, grep("Soil_CNratio", colnames(composite)):=NULL]   #soil C:N is a hypothesized feedback.
composite <- composite[, grep("Soil_N_", colnames(composite)):=NULL]        #soil N content and density are hypothesized feedbacks.
composite <- composite[, grep("Tree_Density", colnames(composite)):=NULL]   #tree density is known.
composite <- composite[, grep("system:index", colnames(composite)):=NULL]   #this is nothing.
composite <- composite[, grep(".geo", colnames(composite)):=NULL]           #this is nothing.
composite <- composite[, grep("Abs_Lat", colnames(composite)):=NULL]        #latitude not going in.
composite <- composite[, grep("Nadir", colnames(composite)):=NULL]          #wtf is nadir reflectance. 
composite <- composite[, grep("Fpar", colnames(composite)):=NULL]           #this is nothing.
composite <- composite[, grep("Npp", colnames(composite)):=NULL]            #NPP not relevant at plot level where we have direct observation.
composite <- composite[, grep("Gpp", colnames(composite)):=NULL]            #GPP not relevant at plot level where we have direct observation.
composite <- composite[, grep("EVI", colnames(composite)):=NULL]            #this is a billion vegetation indices.
composite <- composite[, grep("EnvTexture", colnames(composite)):=NULL]     #this is a few diversity indices.

#Drop plots with missing values.----
composite <- composite[complete.cases(composite),]

#pull out plot code, lat and long.----
plot.ref   <- composite[,.(PLT_CN,Pixel_Lat,Pixel_Long)]
colnames(plot.ref) <- c('PLT_CN','latitude','longitude')
composite  <- composite[,c('PLT_CN','Pixel_Lat','Pixel_Long') := NULL]
#store column names for later.
comp.varnames <- colnames(composite)


#run PCA on full composite. Takes ~20s.-----
tic()
env.pca <- prcomp(composite, scale = T)
var.10 <- sum(summary(env.pca)$importance[2,1:10]) * 100
cat(paste0(round(var.10,2),'% of variation capture in by the first 10 principle components.\n'))
toc()

#visualize. First 10 predictors capture >80%.----
fviz_screeplot(env.pca, addlabels = TRUE, ylim = c(0, 30)) #factoextra package.

#grab first ten PCA values, append to PLT_CN codes w/ lat-lon. Also pull MAT and MAP.----
output <- cbind(plot.ref, 
                composite[,.(WorldClim2_Annual_Mean_Temperature,WorldClim2_Annual_Precipitation)], 
                env.pca$x[,1:10])
#drop in mat-map independently, just because these are useful.
setnames(output,c('WorldClim2_Annual_Mean_Temperature','WorldClim2_Annual_Precipitation'),c('mat','map'))

#assign N deposition.----
output <- cbind(output,extract_ndep(output$longitude,output$latitude))

#assign EPA ecoregions.----
ecoregions <- raster::shapefile(EPA_L2_ecoregions_raster.path)           #load shape file.
ecoregions <- spTransform(ecoregions, CRS('+proj=longlat +datum=WGS84')) #re-project raster to WGS84.
pts <- cbind(output$longitude, output$latitude)                          #bind up your points.
ecoreg.assign <- raster::extract(ecoregions, pts)                         #extract raster.
output$ecoregion <- ecoreg.assign$NA_L2NAME


#append environmental stats to p1 and p2, save.----
p1 <- merge(p1, output, all.x = T, by = 'PLT_CN')
p2 <- merge(p2, output, all.x = T, by = 'PLT_CN')
saveRDS(p1, Product_1.path)
saveRDS(p2, Product_2.path)
write.csv(comp.varnames, composite_variable_names.path)
