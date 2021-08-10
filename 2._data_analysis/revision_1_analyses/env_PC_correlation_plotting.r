#plotting env correlations with principle components.
#compare decomposition indices to existing principle components.
rm(list = ls())
source('paths.r')
source('project_functions/tic_toc.r')
library(corrplot)
library(ggplot2)
library(GGally)
library(data.table)

#load data.----
#plot data.
d <- readRDS(Product_1.path)

#load and clean comp.----
comp <- fread(data_from_composite.path)
comp$PLT_CN <-as.numeric(comp$PLT_CN)
comp <- comp[comp$PLT_CN %in% d$PLT_CN,]
#Drop a few things from the comp linked to the soil feedbacks, or from maps derived from comp.
comp <- comp[, grep("CHELSA" , colnames(comp)):=NULL]        #CHELSEA variables redundant w/ worldclim.
comp <- comp[, grep("Biomass", colnames(comp)):=NULL]        #tree and nematode biomass.
comp <- comp[, grep("Hansen", colnames(comp)):=NULL]         #dont care about year 2000 tree cover.
comp <- comp[, grep("LandCoverClass", colnames(comp)):=NULL] #all are forested.
comp <- comp[, grep("Mycorrhizal", colnames(comp)):=NULL]    #mycorrhizal info is ground sourced.
comp <- comp[, grep("Nematodes", colnames(comp)):=NULL]      #no nematodes.
comp <- comp[, grep("Resolve_Biome", colnames(comp)):=NULL]  #all are forests.
comp <- comp[, grep("Pixel_Area", colnames(comp)):=NULL]     #pixel area not important here.
comp <- comp[, grep("Soil_pH", colnames(comp)):=NULL]        #soil pH is a hypothesized feedback.
comp <- comp[, grep("Soil_CNratio", colnames(comp)):=NULL]   #soil C:N is a hypothesized feedback.
comp <- comp[, grep("Soil_N_", colnames(comp)):=NULL]        #soil N content and density are hypothesized feedbacks.
comp <- comp[, grep("Tree_Density", colnames(comp)):=NULL]   #tree density is known.
comp <- comp[, grep("system:index", colnames(comp)):=NULL]   #this is nothing.
comp <- comp[, grep(".geo", colnames(comp)):=NULL]           #this is nothing.
comp <- comp[, grep("Abs_Lat", colnames(comp)):=NULL]        #latitude not going in.
comp <- comp[, grep("Nadir", colnames(comp)):=NULL]          #wtf is nadir reflectance. 
comp <- comp[, grep("Fpar", colnames(comp)):=NULL]           #this is nothing.
comp <- comp[, grep("Npp", colnames(comp)):=NULL]            #NPP not relevant at plot level where we have direct observation.
comp <- comp[, grep("Gpp", colnames(comp)):=NULL]            #GPP not relevant at plot level where we have direct observation.
comp <- comp[, grep("EVI", colnames(comp)):=NULL]            #this is a billion vegetation indices.
comp <- comp[, grep("EnvTexture", colnames(comp)):=NULL]     #this is a few diversity indices.

#complete case.
comp <- comp[complete.cases(comp),]

#Merge in PC values.
d <- d[d$PLT_CN %in% comp$PLT_CN,]
d <- d[order(match(d$PLT_CN, comp$PLT_CN)),]
sub <- cbind(d[,c(paste0('PC',1:10))], comp)
sub$PLT_CN <- NULL
sub$Pixel_Lat <- NULL
sub$Pixel_Long <- NULL

#make a subsetted corr plot.-----
colnames(sub)[11:ncol(sub)] <- paste0('v',11:ncol(sub))
par(mfrow = c(1,3))
M<-cor(sub[,11:40],sub[,1:10])
corrplot(M)
M<-cor(sub[,41:80],sub[,1:10])
corrplot(M)

