#Supplementary table with environmental layers and associated sources.
rm(list=ls())
source('paths.r')

#load data.----
#load raw layer description from Crowther Lab composite.
comp <- read.csv('required_products_utilities/composite_variable_names.csv', stringsAsFactors = F)

#load variable names actually used in the PCA of the composite.
varnames <- read.csv(composite_variable_names.path, stringsAsFactors = F)
varnames <- varnames[,2]

#subset to variables actually used in analysis.-----
tab <- comp[comp$Old.Band.Name %in% varnames,]
lost.vars <- varnames[!(varnames %in% comp$Old.Band.Name)]

#Some are not completely described in composite. i.e. bands with "month_01" - "Month_12" are entered as "MonthXX".
#This is:
add <- c('EarthEnvCloudCover_MODCF_monthlymean','WorldClim2_H2OVaporPressure_Month','WorldClim2_SolarRadiation_Month','WorldClim2_WindSpeed_Month')
add.out <- list()
for(i in 1:length(add)){
  add.out[[i]] <- comp[grep(add[i], comp$Old.Band.Name),]
}
add.out <- data.frame(do.call(rbind, add.out))
#add into main dataframe.
test <- rbind(tab, add.out)

#Format the table to look line.
tab$Old.Band.Name <- NULL
tab$Layer.Group   <- NULL
tab$Original.Spatial.Resolution <- NULL
tab$Units <- NULL
colnames(tab) <- c('Variable','Description','Source')

#Add in Nitrogen deposition.----
add <- c('Nitrogen Deposition','Sum of 15 year wet and dry nitrogen deposition (as ammonium and nitrate) over 2000-2014','http://nadp.slh.wisc.edu/committees/tdep/reports/TDEPreport15_Final.pdf')
ref <- 'National Atmospheric Deposition Program. (2015). NRSP‐3. Champaign, IL: NADP Program Office, Illinois State Water Survey, University of Illinois.'
#Save output.----
write.csv(tab,'figures/Supplementary_Table_X._environmental_covariates.csv')
