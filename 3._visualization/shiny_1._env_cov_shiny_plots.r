#preparing data and code for shiny application.
rm(list=ls())
source('paths.r')
library(mgcv)
library(data.table)
library(ggplot2)
library(gratia)
library(boot)


#load data.----
d1 <- data.table(readRDS(Product_1.path))
d2 <- data.table(readRDS(Product_2.subset.path))
d1 <- d1[d1$PLT_CN %in% d2$PLT_CN,]
#load fitted model.
model <- readRDS(bimodality_results.path)$model

#grab composite predictors.----
composite <- fread(data_from_composite.path)
#Subset to plots actually in data we're working with.
composite$PLT_CN <-as.numeric(composite$PLT_CN)
composite <- composite[composite$PLT_CN %in% d1$PLT_CN,]
composite <- composite[order(match(composite$PLT_CN, d1$PLT_CN)),]

#Drop covariates from composite that were not used in the analysis.
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
composite$Pixel_Lat  <- NULL
composite$Pixel_Long <- NULL
rownames(composite) <- composite$PLT_CN
composite$PLT_CN <- NULL
composite <- as.data.frame(composite)

#Grab your PC data.----
PC.dat <- as.data.frame(d1[,c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')])

#For every column in composite....
plot.list <- list()
 fit.list <- list()
for(i in 1:ncol(composite)){
  x <- composite[,i] #here is where index enters (instead of 1, i)
  x_by_PC_plots <- list()
  sub_fit_list <- list()
  for(k in 1:ncol(PC.dat)){
    PC.sub <- PC.dat[,k]
    dat.sub <- data.frame(x, PC.sub)
    plot <- ggplot(dat.sub, aes(x = x, y = PC.sub)) + geom_smooth(method = lm) + 
      xlab(colnames(composite)[i]) + 
      ylab(colnames(PC.dat)[k])
    #store in list.
    x_by_PC_plots[[k]] <- plot
     sub_fit_list[[k]] <- lm(PC.sub ~ x)
  }
  plot.list[[i]] <- x_by_PC_plots
   fit.list[[i]] <- sub_fit_list
}
names(plot.list) <- colnames(PC.dat)
names( fit.list) <- colnames(PC.dat)

#Store plots of relationships between PCs and relEM.
relEM_PC_plot.list <- list()
for(i in 1:ncol(PC.dat)){
  var <- paste0('s(',colnames(PC.dat)[i],')')
  plot <- draw(model, select = var)
  relEM_PC_plot.list[[i]] <- plot
}

#Get plots of each env variable and relEM, by mapping through PC regressions and model.
composite_y_plots <- list()
for(i in 1:ncol(composite)){
  x.pred <- seq(min(composite[,i]), max(composite[,i]), by = (max(composite[,i]) - min(composite[,i]))/99)
  #generate PC values for x.
  PC.dat.sub <- list()
  for(k in 1:length(fit.list[[i]])){
    par <- coef(fit.list[[i]][[k]])
    PC.dat.sub[[k]] <- par[1] + x.pred*par[2]
  }
  PC.dat.sub <- data.frame(do.call(cbind, PC.dat.sub))
  colnames(PC.dat.sub) <- colnames(PC.dat)
  PC.dat.sub$ndep <- mean(d1$ndep)
  PC.dat.sub$county.id <- '9_11'
  PC.dat.sub$pred <- inv.logit(predict(model, newdata = PC.dat.sub))
  PC.dat.sub$x    <- x.pred
  #plot
  plot <- ggplot(PC.dat.sub, aes(x = x, y = pred)) + geom_smooth(col = 'black') +
    xlab(colnames(composite)[i]) + 
    ylab('Relative Abundance EM Trees')
  #store in list.
  composite_y_plots[[i]] <- plot
}



