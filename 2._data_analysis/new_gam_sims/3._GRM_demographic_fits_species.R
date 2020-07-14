#Fitting gams w/ separate AM-EM models and PC scores across ecoregions.
rm(list=ls())
library(data.table)
library(mgcv)
library(doParallel)
source('paths.r')
source('project_functions/tic_toc.r')

#set output path.----
output.path <- demographic_fits_gam_species.path

#load growth/mortality/recruitment data.----
d1 <- data.table(readRDS(Product_1.path))
d2 <- data.table(readRDS(Product_2.subset.path))
myc <- read.csv('required_products_utilities/mycorrhizal_SPCD_data.csv')
gym <- read.csv('required_products_utilities/gymnosperm_family_genera.csv')

#function to draw randomly from multivariate normal distribution.
rmvn <- function(n,mu,sig) { ## MVN random deviates
  L <- mroot(sig);m <- ncol(L);
  t(mu + L%*%matrix(rnorm(m*n),m,n)) 
}


#Subset and rename some things.-----
d1 <- d1[d1$PLT_CN %in% d2$PLT_CN,]
d2$PLT_CN <- as.factor(d2$PLT_CN)
setnames(d1,'plot.BASAL','BASAL.plot')
setnames(d2,'plot.BASAL','BASAL.plot')
d1$PLT_CN <- as.character(d1$PLT_CN)
d2$PLT_CN <- as.character(d2$PLT_CN)
d1$county.ID <- paste0(d1$STATECD,'_',d1$COUNTYCD)
d2$county.ID <- paste0(d2$STATECD,'_',d1$COUNTYCD)
d1$county.ID <- as.factor(d1$county.ID)
d2$county.ID <- as.factor(d2$county.ID)

#Get 10 most abundant AM and EM tree species.----
check <- aggregate(DIA.cm ~ SPCD, FUN = length, data = d2)
colnames(check)[2] <- 'count'
#merge in names, gymno and AM-EM.
check <- merge(check,myc[,c('SPCD','MYCO_ASSO','GENUS','SPECIES')], all.x = T)
check$gymno <- ifelse(check$GENUS %in% gym$genus, 1, 0)
#check <- check[!(check$GENUS == 'Ulmus'),] #drop elms (dutch elm disease)
check <- check[order(check$count, decreasing = T),]
top.em <- check[check$MYCO_ASSO == 'ECM',]
top.am <- check[check$MYCO_ASSO == 'AM' ,]
n.spp <- 11 #number of species to grab of AM and EM functional types. Will be cut if not 1k plots present.
top.em <- top.em[1:n.spp,]
top.am <- top.am[1:n.spp,]
top.20 <- rbind(top.am, top.em)
top.20 <- top.20[top.20$count > 1000,] #cut out species w/o 1k plots.
top.20$name <- paste(top.20$GENUS, top.20$SPECIES, sep = ' ')
spp.name.lab <- top.20$name

#how much basal area of AM and EM trees is this?----
  em.basal <- sum(d2[d2$em == 1,]$BASAL, na.rm = T)
em10.basal <- sum(d2[d2$em == 1 & d2$SPCD %in% top.20$SPCD,]$BASAL, na.rm = T)
  am.basal <- sum(d2[d2$em == 0,]$BASAL, na.rm = T)
am10.basal <- sum(d2[d2$em == 0 & d2$SPCD %in% top.20$SPCD,]$BASAL, na.rm = T)
em.rel <- round((em10.basal/em.basal)*100, 2)
am.rel <- round((am10.basal/am.basal)*100, 2)
cat(paste0('These species represent ',em.rel,'% of all EM basal area and ',am.rel,'% of all AM basal area.\n'))  

#Set global k for gam fits.----
kk <- 5
bs <- 'cr'

#Fit growth, recruitment and mortality models by species.----
output <- list()
tic()
for(i in 1:length(spp.name.lab)){
  #subset data to species of interest.
  d2.sub <- d2[d2$SPCD == top.20$SPCD[i],]
  d1.sub <- d1[d1$PLT_CN %in% d2.sub$PLT_CN,]
  d2.sub$recruit <- ifelse(is.na(d2.sub$PREVDIA.cm), 1, 0)
  #calculate species recruits.
  recruits <- aggregate(recruit ~ PLT_CN, FUN = sum, data = d2.sub)
  d1.sub$recruit <- NULL
  d1.sub <- merge(d1.sub, recruits, by = 'PLT_CN')
  basal.consp <- aggregate(BASAL ~ PLT_CN, FUN = sum, data=d2.sub[d2.sub$recruit == 0,])
  colnames(basal.consp) <- c('PLT_CN','BASAL.consp')
  d1.sub <- merge(d1.sub,basal.consp, all.x = T)
  d1.sub$BASAL.consp <- ifelse(is.na(d1.sub$BASAL.consp), 0, d1.sub$BASAL.consp)
  
  #fit gam models.
  R.mod <- bam(recruit    ~  s(relEM, k=kk, bs=bs) + s(BASAL.consp, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs)+
                 s(county.ID, bs = 're') + 
                 s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d1.sub, family = 'poisson')
  M.mod <- bam(mortality  ~  s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, k=kk, bs=bs)+
                 s(county.ID, bs = 're') + 
                 s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d2.sub, family = 'binomial')
  G.mod <- bam(DIA.cm   ~  s(relEM, k=kk, bs=bs) + s(ndep, k=kk, bs=bs) + s(BASAL.plot, k=kk, bs=bs) + s(stem.density, k=kk, bs=bs) + s(PREVDIA.cm, kk=kk, bs=bs)+
                 s(county.ID, bs = 're') + 
                 s(PC1, k=kk, bs=bs) + s(PC2, k=kk, bs=bs) + s(PC3, k=kk, bs=bs) + s(PC4, k=kk, bs=bs) + s(PC5, k=kk, bs=bs) + s(PC6, k=kk, bs=bs) + s(PC7, k=kk, bs=bs) + s(PC8, k=kk, bs=bs) + s(PC9, k=kk, bs=bs) + s(PC10, k=kk, bs=bs), 
                  data = d2.sub)
  
  #Grab plot environmental covariates for reference.
  #all plot-level environmental covariates. we will sample this later when simulating forests.
  all.cov <- d1.sub[,c('BASAL.plot','BASAL.consp','stem.density','ndep','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')]
  #mean plot-level environmental covariates.
  cov <- colMeans(all.cov)
  cov <- c(mean(d2.sub$PREVDIA.cm, na.rm=T),cov)
  names(cov)[1] <- 'PREVDIA.cm'
  pred.frame <- data.frame(seq(0 ,1 ,by = 0.01), t(cov))
  colnames(pred.frame)[1] <- 'relEM'
  g.pred <- predict(G.mod, newdata = pred.frame, se.fit = T, exclude = c("s(county.ID)"), newdata.guaranteed = T)
  m.pred <- predict(M.mod, newdata = pred.frame, se.fit = T, exclude = c("s(county.ID)"), newdata.guaranteed = T)
  r.pred <- predict(R.mod, newdata = pred.frame, se.fit = T, exclude = c("s(county.ID)"), newdata.guaranteed = T)
  
  #get contrast by drawing from posterior.
  newdat <- pred.frame[c(1, nrow(pred.frame)),]
  #covert X predictors to match knots.
  Xp.g <- predict(G.mod, newdata = newdat, type = 'lpmatrix', exclude = c("s(county.ID)"), newdata.guaranteed = T)
  Xp.m <- predict(M.mod, newdata = newdat, type = 'lpmatrix', exclude = c("s(county.ID)"), newdata.guaranteed = T)
  Xp.r <- predict(R.mod, newdata = newdat, type = 'lpmatrix', exclude = c("s(county.ID)"), newdata.guaranteed = T)
  
  #Draw matrix of correlated predictors from posterior.
  g.par <- rmvn(1000,coef(G.mod),G.mod$Vp)
  m.par <- rmvn(1000,coef(M.mod),M.mod$Vp)
  r.par <- rmvn(1000,coef(R.mod),R.mod$Vp)
  
  #Multiply predictors by each posterior draw of parameters.
  g.pred <- list()
  m.pred <- list()
  r.pred <- list()
  for(j in 1:1000){
    g.pred[[j]] <- t(Xp.g %*% g.par[j,])
    m.pred[[j]] <- t(Xp.m %*% m.par[j,])
    r.pred[[j]] <- t(Xp.r %*% r.par[j,])
  }
  g.pred <- do.call(rbind, g.pred)
  m.pred <- do.call(rbind, m.pred)
  r.pred <- do.call(rbind, r.pred)
  
  #convert mortality to survival.
  m.pred <- boot::logit(1 - boot::inv.logit(m.pred))
  
  #grab summary stats.
  g.mu <-     mean(g.pred[,2] - g.pred[,1])
  g.sd <-       sd(g.pred[,2] - g.pred[,1])
  g.95 <- quantile(g.pred[,2] - g.pred[,1], probs = c(0.025, 0.975))
  m.mu <-     mean(m.pred[,2] - m.pred[,1])
  m.sd <-       sd(m.pred[,2] - m.pred[,1])
  m.95 <- quantile(m.pred[,2] - m.pred[,1], probs = c(0.025, 0.975))
  r.mu <- mean(r.pred[,2] - r.pred[,1])
  r.sd <-   sd(r.pred[,2] - r.pred[,1])
  r.95 <- quantile(r.pred[,2] - r.pred[,1], probs = c(0.025, 0.975))
  post.contrast <- list(g.mu, g.sd, g.95, m.mu, m.sd, m.95, r.mu, r.sd, r.95)
  names(post.contrast) <- c('g.mu','g.sd','g.95','m.mu','m.sd','m.95','r.mu','r.sd','r.95')
  
  #store in list and report.
  spp.out <- list(G.mod, M.mod, R.mod, g.pred, m.pred, r.pred, cov,post.contrast)
  names(spp.out) <- c('G.mod','M.mod','R.mod','G.pred','M.pred','R.pred','cov','post.contrast')
  output[[i]] <- spp.out
  msg <- paste0(spp.name.lab[i],' fit. ',i,' of ',length(spp.name.lab),' species fits complete.\n')
  cat(msg); toc()
}
names(output) <- spp.name.lab
output$spp.key <- top.20

#Save models and size categories.----
saveRDS(output, output.path, version = 2)
cat('Output saved, script complete.\n')

#Plotting all as panels.----
plotting <- F
if(plotting == T){
  #Plot recruitment models.
  par(mfrow = c(4,5))
  N <- nrow(top.20)
  for(i in 1:N){
    new.dat <- data.frame(seq(0 ,1 ,by = 0.01), t(output[[i]]$cov))
    colnames(new.dat)[1] <- 'relEM'
    pred.r <- predict(output[[i]]$R.mod, newdata = new.dat)
    plot(pred.r ~ new.dat$relEM, bty = 'l', ylim = c(min(c(pred.r,pred.r)), max(c(pred.r,pred.r))), main = spp.name.lab[i])
  }
  
  #Plot mortality models.
  par(mfrow = c(4,5))
  for(i in 1:N){
    new.dat <- data.frame(seq(0 ,1 ,by = 0.01), t(output[[i]]$cov))
    colnames(new.dat)[1] <- 'relEM'
    pred.m <- predict(output[[i]]$M.mod, newdata = new.dat)
    plot(pred.m ~ new.dat$relEM, bty = 'l', ylim = c(min(c(pred.m,pred.m)), max(c(pred.m,pred.m))), main = spp.name.lab[i])
  }
  
  #Plot Growth models.
  par(mfrow = c(4,5))
  for(i in 1:N){
    new.dat <- data.frame(seq(0 ,1 ,by = 0.01), t(output[[i]]$cov))
    colnames(new.dat)[1] <- 'relEM'
    pred <- predict(output[[i]]$G.mod, newdata = new.dat)
    plot(pred ~ new.dat$relEM, bty = 'l', ylim = c(min(c(pred,pred)), max(c(pred,pred))), main = spp.name.lab[i])
  }
}

#Get recruitment and mortality delta values and plot.----
r.delta <- list()
m.delta <- list()
for(i in 1:length(spp.name.lab)){
  new.dat <- data.frame(seq(0 ,1 ,by = 0.01), t(output[[i]]$cov))
  colnames(new.dat)[1] <- 'relEM'
  pred.r <- predict(output[[i]]$R.mod, newdata = new.dat)
  pred.m <- predict(output[[i]]$M.mod, newdata = new.dat)
  r.delta[[i]] <- pred.r[length(pred.r)] - pred.r[1]
  m.delta[[i]] <- pred.m[length(pred.m)] - pred.m[1]
}
r.delta <- unlist(r.delta)
m.delta <- unlist(m.delta)
top.20$r.delta <- r.delta
top.20$m.delta <- m.delta

par(mfrow = c(1,2))
N <- nrow(top.20)
#Recruitment advantage in EM relative to AM forests by species.
top.20 <- top.20[order(top.20$r.delta),]
plot(top.20$r.delta, bty = 'l', ylab = 'Recruit Advantage EM forest'  , xlab = NA, pch = ifelse(top.20$MYCO_ASSO == 'AM',1, 16), xaxt = 'n', cex = 1.3)
abline(h = 0, lwd =2, lty = 2)  
axis(1, at=1:N, labels = F)
text(top.20$name,cex = 0.8, x = c(1:N), y = min(top.20$r.delta)*1.15, xpd = T, srt = 90, adj = 1)


#Mortality advantage in EM relative to AM forests by species.
top.20 <- top.20[order(top.20$m.delta),]
plot(top.20$m.delta, bty = 'l', ylab = 'Mortality Advantage EM forest', xlab = NA, pch = ifelse(top.20$MYCO_ASSO == 'AM',1, 16), xaxt = 'n', cex = 1.3)
abline(h = 0, lwd =2, lty = 2)  
axis(1, at=1:N, labels = F)
text(top.20$name,cex = 0.8, x = c(1:N), y = min(top.20$m.delta)*1.15, xpd = T, srt = 90, adj = 1)
