#parameterizing mycorrhizal demographic models of growth, recruitment and mortality.
rm(list=ls())
library(data.table)
library(mgcv)
source('paths.r')

#set output path.----
output.path <- demographic_fits.path

#load growth/mortality/recruitment data.----
d1 <- data.table(readRDS(Product_1.path))
d2 <- data.table(readRDS(Product_2.subset.path))

#Subset and rename some things.-----
d1 <- d1[d1$PLT_CN %in% d2$PLT_CN,]
setnames(d1,'plot.BASAL','BASAL.plot')
setnames(d2,'plot.BASAL','BASAL.plot')

#Fit growth, recruitment and mortality models.----
#Environmental models without feedbacks.
G.mod    <- gam(DIA.cm     ~ em + em*ndep      + s(mat, k = 3) + s(map, k = 3) + s(PREVDIA.cm       ) + s(BASAL.plot) + s(stem.density), data = d2[DIA.cm > 0,])
M.mod    <- gam(mortality  ~ em + em*ndep      + s(mat, k = 3) + s(map, k = 3) + s(PREVDIA.cm, k = 5) + s(BASAL.plot) + s(stem.density), data = d2, family = 'binomial')
R.mod.am <- gam(recruit.am ~ am.density + ndep + s(mat, k = 3) + s(map, k = 3)                        + s(BASAL.plot) + s(stem.density), data = d1, family = 'poisson')
R.mod.em <- gam(recruit.em ~ em.density + ndep + s(mat, k = 3) + s(map, k = 3)                        + s(BASAL.plot) + s(stem.density), data = d1, family = 'poisson')
n.feedback <- list(G.mod, M.mod, R.mod.am, R.mod.em)
names(n.feedback) <- c('G.mod','M.mod','R.mod.am','R.mod.em')

#Environmental models with feedbacks.
G.mod    <- gam(DIA.cm     ~ em*relEM + em*ndep      + s(mat, k = 3) + s(map, k = 3) + s(PREVDIA.cm       ) + s(BASAL.plot) + s(stem.density), data = d2[DIA.cm > 0,])
M.mod    <- gam(mortality  ~ em*relEM + em*ndep      + s(mat, k = 3) + s(map, k = 3) + s(PREVDIA.cm, k = 5) + s(BASAL.plot) + s(stem.density), data = d2, family = 'binomial')
R.mod.am <- gam(recruit.am ~ am.density*relEM + ndep + s(mat, k = 3) + s(map, k = 3)                        + s(BASAL.plot) + s(stem.density), data = d1, family = 'poisson')
R.mod.em <- gam(recruit.em ~ em.density*relEM + ndep + s(mat, k = 3) + s(map, k = 3)                        + s(BASAL.plot) + s(stem.density), data = d1, family = 'poisson')
y.feedback <- list(G.mod, M.mod, R.mod.am, R.mod.em)
names(y.feedback) <- c('G.mod','M.mod','R.mod.am','R.mod.em')

#Grab plot environmental covariates for reference.----
#mean plot-level environmental covariates.
cov <- c(mean(d1$mat), mean(d1$map), mean(d1$ndep))
names(cov) <- c('mat','map','ndep')
#all plot-level environmental covariates. we will sample this later when simulating forests.
all.cov <- d1[,c('mat','map','ndep')]

#Save models and size categories.----
output <- list(n.feedback, y.feedback, cov, all.cov)
names(output) <- c('n.feedback', 'y.feedback','env.cov','all.cov')
saveRDS(output, output.path, version = 2)
