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

#load data.----
d <- readRDS(Product_1.path)
composite <- fread(data_from_composite.path)
#Subset to plots actually in data we're working with.
#composite <- output[output$PLT_CN %in% d$PLT_CN,] #broken because the PLT_CN I sent vs. received don't match. Asking Johan.


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
composite <- composite[, grep("Npp", colnames(composite)):=NULL]            #this is NPP not relevant at plot level where we have direct observation.
composite <- composite[, grep("Gpp", colnames(composite)):=NULL]            #this is NPP not relevant at plot level where we have direct observation.
composite <- composite[, grep("_EVI", colnames(composite)):=NULL]           #this is a billion vegetation indices.
composite <- composite[, grep("EnvTexture", colnames(composite)):=NULL]     #this is a few diversity indices.

#Drop plots with missing values.----
composite <- composite[complete.cases(composite),]

#pull out plot code, lat and long.----
plot.ref   <- composite[,.(PLT_CN,Pixel_Lat,Pixel_Long)]
colnames(plot.ref) <- c('PLT_CN','latitude','longitude')
composite  <- composite[,c('PLT_CN','Pixel_Lat','Pixel_Long') := NULL]

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

#save output.----
