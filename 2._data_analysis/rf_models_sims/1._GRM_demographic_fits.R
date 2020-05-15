#parameterizing mycorrhizal demographic models of growth, recruitment and mortality.
rm(list=ls())
library(data.table)
library(randomForest)
library(rfUtilities)
library(doParallel)
source('paths.r')
source('project_functions/tic_toc.r')
source('project_functions/parallel_rf.r')
 
#set output path.----
output.path <- rf_demographic_fits.path

#load growth/mortality/recruitment data.----
d1 <- data.table(readRDS(Product_1.path))
d2 <- data.table(readRDS(Product_2.subset.path))

#Subset and rename some things.-----
d1 <- d1[d1$PLT_CN %in% d2$PLT_CN,]
setnames(d1,'plot.BASAL','BASAL.plot')
setnames(d2,'plot.BASAL','BASAL.plot')

#organize & complete case data within AM-EM divisions.-----
grow.dat.all <- d2[          ,.(DIA.cm,PREVDIA.cm,relEM,BASAL.plot,stem.density,ndep,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)]
grow.dat.am  <- d2[d2$em == 0,.(DIA.cm,PREVDIA.cm,relEM,BASAL.plot,stem.density,ndep,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)]
grow.dat.em  <- d2[d2$em == 1,.(DIA.cm,PREVDIA.cm,relEM,BASAL.plot,stem.density,ndep,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)]
grow.dat.all <- grow.dat.all[complete.cases(grow.dat.all),]
grow.dat.am  <- grow.dat.am [complete.cases(grow.dat.am ),]
grow.dat.em  <- grow.dat.em [complete.cases(grow.dat.em ),]
#mortality data frames.
mort.dat.all <- d2[          ,.(mortality,PREVDIA.cm,relEM,BASAL.plot,stem.density,ndep,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)]
mort.dat.am  <- d2[d2$em == 0,.(mortality,PREVDIA.cm,relEM,BASAL.plot,stem.density,ndep,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)]
mort.dat.em  <- d2[d2$em == 1,.(mortality,PREVDIA.cm,relEM,BASAL.plot,stem.density,ndep,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)]
mort.dat.all <- mort.dat.all[complete.cases(mort.dat.all),]
mort.dat.am  <- mort.dat.am [complete.cases(mort.dat.am ),]
mort.dat.em  <- mort.dat.em [complete.cases(mort.dat.em ),]
#recruitment data frames.
recr.dat.all <- d1[,.(recruit.am,recruit.em,relEM,BASAL.plot,stem.density,ndep,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)]
recr.dat.all <- recr.dat.all[complete.cases(recr.dat.all),]
recr.dat.am  <- recr.dat.all[,!"recruit.em"]
recr.dat.em  <- recr.dat.all[,!"recruit.am"]

#Fit models.-----
#growth models.
tic()
grow.mod.am <- parallel_rf(DIA.cm ~ ., data = grow.dat.am)
grow.mod.em <- parallel_rf(DIA.cm ~ ., data = grow.dat.em)
cat('Growth models fit.\n')
toc()

#mortality models.
tic()
mort.mod.am <- parallel_rf(mortality ~ ., data = mort.dat.am)
mort.mod.em <- parallel_rf(mortality ~ ., data = mort.dat.em)
cat('Mortality models fit.\n')
toc()

#recuritment models.
tic()
recr.mod.am <- parallel_rf(recruit.am ~ ., data = recr.dat.am)
recr.mod.em <- parallel_rf(recruit.em ~ ., data = recr.dat.em)
cat('Recruitment models fit.\n')
toc()

#Wrap output and save.----
#models.
model.list <- list(grow.mod.am, grow.mod.em,
                   mort.mod.am, mort.mod.em,
                   recr.mod.am, recr.mod.em)
names(model.list) <- c('grow.mod.am', 'grow.mod.em',
                       'mort.mod.am', 'mort.mod.em',
                       'recr.mod.am', 'recr.mod.em')
#data.
 data.list <- list(grow.dat.all, grow.dat.am, grow.dat.em,
                   mort.dat.all, mort.dat.am, mort.dat.em,
                   recr.dat.all, recr.dat.am, recr.dat.em)
names(data.list) <- c('grow.dat.all', 'grow.dat.am', 'grow.dat.em',
                      'mort.dat.all', 'mort.dat.am', 'mort.dat.em',
                      'recr.dat.all', 'recr.dat.am', 'recr.dat.em')
#full output.
output <- list(model.list, data.list)
names(output) <- c('models','data')
saveRDS(output, output.path)
cat('Script complete.\n')

#Visualizing relEM effect.----
#RECRUITMENT.
new.dat <- colMeans(recr.dat.all)
new.dat <- new.dat[!(names(new.dat) %in% c('relEM','recruit.em'))]
new.dat <- data.frame(seq(0,1,by=0.01), t(new.dat))
colnames(new.dat)[1] <- 'relEM'
recr.pred.em <- predict(recr.mod.em, newdata = new.dat)
recr.pred.am <- predict(recr.mod.am, newdata = new.dat)
#plot
par(mfrow = c(3,2))
plot(recr.pred.am ~ new.dat$relEM, bty = 'l', ylim = c(0, max(c(recr.pred.am, recr.pred.em))), main = 'AM Recruitment')
plot(recr.pred.em ~ new.dat$relEM, bty = 'l', ylim = c(0, max(c(recr.pred.am, recr.pred.em))), main = 'EM Recruitment')


#MORTALITY
new.dat <- colMeans(mort.dat.all)
new.dat <- new.dat[!(names(new.dat) %in% c('relEM','mortality'))]
new.dat <- data.frame(seq(0,1,by=0.01), t(new.dat))
colnames(new.dat)[1] <- 'relEM'
mort.pred.em <- predict(mort.mod.em, newdata = new.dat)
mort.pred.am <- predict(mort.mod.am, newdata = new.dat)
#plot
#par(mfrow = c(1,2))
plot(mort.pred.am ~ new.dat$relEM, bty = 'l', ylim = c(0, max(c(mort.pred.am, mort.pred.em))), main = 'AM Mortality Probability')
plot(mort.pred.em ~ new.dat$relEM, bty = 'l', ylim = c(0, max(c(mort.pred.am, mort.pred.em))), main = 'EM Mortality Probability')
plot(mort.pred.am ~ new.dat$relEM, bty = 'l', ylim = c(0, 0.1), main = 'AM Mortality Probability') #Works in more reasonable probability range too. Bounded because if a EM tree is present and dies, relEM cant be 0.
plot(mort.pred.em ~ new.dat$relEM, bty = 'l', ylim = c(0, 0.1), main = 'EM Mortality Probability')


#GROWTH
new.dat <- colMeans(grow.dat.all)
new.dat <- new.dat[!(names(new.dat) %in% c('relEM','DIA.cm'))]
new.dat <- data.frame(seq(0,1,by=0.01), t(new.dat))
colnames(new.dat)[1] <- 'relEM'
grow.pred.em <- predict(grow.mod.em, newdata = new.dat)
grow.pred.am <- predict(grow.mod.am, newdata = new.dat)
#plot
plot(grow.pred.am ~ new.dat$relEM, bty = 'l', ylim = c(min(new.dat$PREVDIA.cm), max(c(grow.pred.am, grow.pred.em))), main = 'AM Predicted Diameter')
plot(grow.pred.em ~ new.dat$relEM, bty = 'l', ylim = c(min(new.dat$PREVDIA.cm), max(c(grow.pred.am, grow.pred.em))), main = 'EM Predicted Diameter')


