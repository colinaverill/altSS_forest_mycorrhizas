#simulate density dependence + freq dependence.
rm(list=ls())
library(data.table)
source('paths.r')

#set output path.----
output.path <- recruit_power_analysis.path

#load plot data and growth/mortality/recruitment data.----
d1 <- data.table(readRDS(Product_1.path))
d2 <- data.table(readRDS(Product_2.subset.path))

#Subset plot data and rename some things.-----
d1 <- d1[d1$PLT_CN %in% d2$PLT_CN,]
d2$PLT_CN <- as.factor(d2$PLT_CN)
d2$county.ID <- as.factor(paste0(d2$STATECD,'_',d2$COUNTYCD))
d1$county.ID <- as.factor(paste0(d1$STATECD,'_',d1$COUNTYCD))
setnames(d1,'plot.BASAL','BASAL.plot')
setnames(d2,'plot.BASAL','BASAL.plot')
setnames(d1,'BASAL.ECM','BASAL.em')
setnames(d1,'BASAL.AM' ,'BASAL.am')

#Simulate recruitment at a prescibed effect size, different levels of sampling effort..----
effort <- c(100, 500, 1000, 3000, 6000) #sampling effort (number of plots).
par <- c(-0.6, -0.0005, 2)              #prescribed, true paramters (intercept, basal area EM effect, relative abundance EM effect).
N <- 1000                               #number of simulations to run per level of sampling effort.
effort.out <- list()
for(k in 1:length(effort)){                                                #for each elvel of sampling effort...
  effort.result <- list()
  for(i in 1:N){                                                           #run N simulations...
    sub <- d1[sample(nrow(d1), effort[k]),]                                #where you sample the original dataframe (preserving correlation btwn basal area and relatibe abundance)
    lambda <- par[1] + par[2]*sub$BASAL.em + par[3]*sub$relEM              #multiply predictors by parameters. 
    sub$recruits <- rpois(length(lambda), exp(lambda))                     #pass lambda estimates through poisson distribution, accounting for log-link.
    fit <- glm(recruits ~ BASAL.em + relEM, data = sub, family = poisson)  #fit a poisson model to simulated recruitment.
    effort.result[[i]] <- coef(fit)                                        #save the parameter estimates.
  }
  effort.result <- data.frame(do.call(rbind, effort.result))               #for a given level of sample effort, convert output into a dataframe.
  colnames(effort.result) <- c('m0','m1','m2')
  effort.out[[k]] <- effort.result
  msg <- paste0(k,' of ',length(effort),' levels of fit complete.\n')
  cat(msg)
}

#calculate mean and 95% quantiles of each level of effort for relative EM effect.-----
results <- list()
for(i in 1:length(effort.out)){
  dat <- effort.out[[i]]$m2
  mu  <- mean(dat, na.rm = T)
  lo.hi <- quantile(dat, probs = c(0.025, 0.975))
  results[[i]] <- c(effort[i],mu,lo.hi)
}
results <- data.frame(do.call(rbind, results))
colnames(results) <- c('N','mu','lo95','hi95')

#save output.----
saveRDS(results, output.path)


#plot power analysis.----
plot.it <- F
if(plot.it == T){
  par(mfrow = c(1,1))
  x <- c(1:nrow(results))
  plot(results$mu ~ x, pch = 16, bty='l',cex =2,
       ylim = c(min(results$lo95)*1.05, max(results$hi95)*1.05),
       xaxt='n', ylab=NA, xlab=NA)
  abline(h=2, lwd = 2, lty = 2)
  arrows(x, results$lo95, x, results$hi95, length=0.04, angle=90, code=3)
  mtext('parameter estimate',side = 2, line = 2.3)
  mtext('sampling effort', side = 1, line = 3)
  text(cex=1, x= c(1:nrow(results)), y=-0.0, effort, srt = 45, xpd=TRUE, adj=1)
}
