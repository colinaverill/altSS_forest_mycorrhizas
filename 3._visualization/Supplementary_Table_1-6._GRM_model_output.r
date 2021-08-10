#generate statistical output tables for AM nd EM GRM gam fits.
rm(list=ls())
source('paths.r')
library(tidyverse)
library(pixiedust)

#load data.-----
d <- readRDS(demographic_fits_gam_separate_plus_re_county.path)

#AM Growth.----
fit <- d$y.feedback$G.mod.am
out <-
dust(fit) %>%
  sprinkle(col = 2:4, round = 2) %>%
  sprinkle(col = 5  , round = 3) %>%
  sprinkle_colnames(term = 'Term',
                    statistic = 'F-statistic',
                    p.value = 'P-value',
                    edf = 'Estimated D.F.') %>%
  sprinkle(cols = 'term', replace = c('Relative Abundance EM trees',
                                      'N Deposition',
                                      'Plot Basal Area',
                                      'Stem Density',
                                      'Previous Census Diameter',
                                      'County (Random Effect)',
                                      'PC1','PC2','PC3','PC4','PC5',
                                      'PC6','PC7','PC8','PC9','PC10'))

out <-as.data.frame(out)
out$ref.df <- NULL
AM.growth <- out

#EM Growth.----
fit <- d$y.feedback$G.mod.em
out <-
  dust(fit) %>%
  sprinkle(col = 2:4, round = 2) %>%
  sprinkle(col = 5  , round = 3) %>%
  sprinkle_colnames(term = 'Term',
                    statistic = 'F-statistic',
                    p.value = 'P-value',
                    edf = 'Estimated D.F.') %>%
  sprinkle(cols = 'term', replace = c('Relative Abundance EM trees',
                                      'N Deposition',
                                      'Plot Basal Area',
                                      'Stem Density',
                                      'Previous Census Diameter',
                                      'County (Random Effect)',
                                      'PC1','PC2','PC3','PC4','PC5',
                                      'PC6','PC7','PC8','PC9','PC10'))

out <-as.data.frame(out)
out$ref.df <- NULL
EM.growth <- out

#AM Recruitment.----
fit <- d$y.feedback$R.mod.am
out <-
  dust(fit) %>%
  sprinkle(col = 2:4, round = 2) %>%
  sprinkle(col = 5  , round = 3) %>%
  sprinkle_colnames(term = 'Term',
                    statistic = 'F-statistic',
                    p.value = 'P-value',
                    edf = 'Estimated D.F.') %>%
  sprinkle(cols = 'term', replace = c('Relative Abundance EM trees',
                                      'Basal Area AM Trees',
                                      'N Deposition',
                                      'Plot Basal Area',
                                      'Stem Density',
                                      'County (Random Effect)',
                                      'PC1','PC2','PC3','PC4','PC5',
                                      'PC6','PC7','PC8','PC9','PC10'))

out <-as.data.frame(out)
out$ref.df <- NULL
AM.recruitment <- out

#EM Recruitment.----
fit <- d$y.feedback$R.mod.em
out <-
  dust(fit) %>%
  sprinkle(col = 2:4, round = 2) %>%
  sprinkle(col = 5  , round = 3) %>%
  sprinkle_colnames(term = 'Term',
                    statistic = 'F-statistic',
                    p.value = 'P-value',
                    edf = 'Estimated D.F.') %>%
  sprinkle(cols = 'term', replace = c('Relative Abundance EM trees',
                                      'Basal Area EM Trees',
                                      'N Deposition',
                                      'Plot Basal Area',
                                      'Stem Density',
                                      'County (Random Effect)',
                                      'PC1','PC2','PC3','PC4','PC5',
                                      'PC6','PC7','PC8','PC9','PC10'))

out <-as.data.frame(out)
out$ref.df <- NULL
EM.recruitment <- out

#AM Mortality.----
fit <- d$y.feedback$M.mod.am
out <-
  dust(fit) %>%
  sprinkle(col = 2:4, round = 2) %>%
  sprinkle(col = 5  , round = 3) %>%
  sprinkle_colnames(term = 'Term',
                    statistic = 'F-statistic',
                    p.value = 'P-value',
                    edf = 'Estimated D.F.') %>%
  sprinkle(cols = 'term', replace = c('Relative Abundance EM trees',
                                      'N Deposition',
                                      'Plot Basal Area',
                                      'Stem Density',
                                      'Previous Census Diameter',
                                      'County (Random Effect)',
                                      'PC1','PC2','PC3','PC4','PC5',
                                      'PC6','PC7','PC8','PC9','PC10'))

out <-as.data.frame(out)
out$ref.df <- NULL
AM.mortality <- out

#EM Mortality.----
fit <- d$y.feedback$M.mod.em
out <-
  dust(fit) %>%
  sprinkle(col = 2:4, round = 2) %>%
  sprinkle(col = 5  , round = 3) %>%
  sprinkle_colnames(term = 'Term',
                    statistic = 'F-statistic',
                    p.value = 'P-value',
                    edf = 'Estimated D.F.') %>%
  sprinkle(cols = 'term', replace = c('Relative Abundance EM trees',
                                      'N Deposition',
                                      'Plot Basal Area',
                                      'Stem Density',
                                      'Previous Census Diameter',
                                      'County (Random Effect)',
                                      'PC1','PC2','PC3','PC4','PC5',
                                      'PC6','PC7','PC8','PC9','PC10'))

out <-as.data.frame(out)
out$ref.df <- NULL
EM.mortality <- out

#Write .csv file output.----
#make directory if it doesn't exist
dir <- 'figures/supp._table_2-7/'
system(paste0('mkdir -p ',dir))

#write csvs.
write.csv(AM.growth     , paste0(dir,     'AM_growth.csv'), row.names = F)
write.csv(EM.growth     , paste0(dir,     'EM_growth.csv'), row.names = F)
write.csv(AM.recruitment, paste0(dir,'AM_recruitment.csv'), row.names = F)
write.csv(EM.recruitment, paste0(dir,'EM_recruitment.csv'), row.names = F)
write.csv(AM.mortality  , paste0(dir,  'AM_mortality.csv'), row.names = F)
write.csv(EM.mortality  , paste0(dir,  'EM_mortality.csv'), row.names = F)

#end script.
