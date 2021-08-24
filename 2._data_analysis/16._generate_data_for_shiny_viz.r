#preparing data and code for shiny application.
rm(list=ls())
source('paths.r')
library(mgcv)
library(data.table)
library(ggplot2)
library(gratia)
library(boot)
library(mgcv)

#set output path.----
output.path <- shiny_viz_data.path

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

#For every column in composite get linear correlations w/ our 10 PCs.----
#plot.list <- list()
fit.list <- list()
for(i in 1:ncol(composite)){
  x <- composite[,i] #here is where index enters (instead of 1, i)
  sub_fit_list <- list()
  for(k in 1:ncol(PC.dat)){
    PC.sub <- PC.dat[,k]
    dat.sub <- data.frame(x, PC.sub)
    #store in list.
    #sub_fit_list[[k]] <- lm(PC.sub ~ x)
    sub_fit_list[[k]] <- gam(PC.sub ~ s(x, k=7))
  }
  fit.list[[i]] <- sub_fit_list
}
names( fit.list) <- colnames(PC.dat)

#Store plots of relationships between PCs and relEM.----
#This may be unnecessary.
#relEM_PC_plot.list <- list()
#for(i in 1:ncol(PC.dat)){
#  var <- paste0('s(',colnames(PC.dat)[i],')')
#  plot <- draw(model, select = var)
#  relEM_PC_plot.list[[i]] <- plot
#}

#Generate predictions of how each env variable is correlated with relEM, by mapping through PC regressions and then final model.----
composite_preds <- list()
for(i in 1:ncol(composite)){
  x.pred <- seq(min(composite[,i]), max(composite[,i]), by = (max(composite[,i]) - min(composite[,i]))/99)
  #generate PC values for x.
  PC.dat.sub <- list()
  for(k in 1:length(fit.list[[i]])){
    df <- data.frame(x.pred)
    colnames(df) <- 'x'
    PC.dat.sub[[k]] <- predict(fit.list[[i]][[k]], newdata=df)
    #par <- coef(fit.list[[i]][[k]])
    #PC.dat.sub[[k]] <- par[1] + x.pred*par[2]
  }
  PC.dat.sub <- data.frame(do.call(cbind, PC.dat.sub))
  colnames(PC.dat.sub) <- colnames(PC.dat)
  PC.dat.sub$ndep <- mean(d1$ndep)
  PC.dat.sub$county.id <- '9_11'
  PC.dat.sub$x    <- x.pred
  #generate predictions.
  pred.EM <- predict(model, newdata = PC.dat.sub, se.fit = T)
  PC.dat.sub$y    <- inv.logit(pred.EM$fit)
  PC.dat.sub$lo95 <- inv.logit(pred.EM$fit - pred.EM$se.fit)
  PC.dat.sub$hi95 <- inv.logit(pred.EM$fit + pred.EM$se.fit)
  #store output.
  composite_preds[[i]] <- PC.dat.sub[,c('x','y','lo95','hi95')]
}
names(composite_preds) <- colnames(composite)

#save output.----
saveRDS(composite_preds, output.path)

#plot it?----
i=143 #mean annual temperature is 143.
PC.dat.sub <- composite_preds[[i]]
#plot
  limy <- c(min(PC.dat.sub$lo95)*0.95, max(PC.dat.sub$hi95)*1.05)
  plot(y ~ x, data = PC.dat.sub, bty='l', cex = 0, ylim=limy,
       ylab = 'Relative Abundance EM Trees',
       xlab = names(composite_preds)[i])
  lines(smooth.spline(PC.dat.sub$y ~ PC.dat.sub$x), lwd = 2)
  polygon(c(PC.dat.sub$x, rev(PC.dat.sub$x)),c(PC.dat.sub$hi95, rev(PC.dat.sub$lo95)), col=adjustcolor('black', 0.3), lty=0)

