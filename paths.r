#paths.r
#Main file dump directories.----
host <- system('hostname', intern = T)
storage.dir <- '/projectnb/talbot-lab-data/caverill/altSS_forest_mycorrhizas_data/'
if(host == 'pecan2')                    {storage.dir <- '/fs/data3/caverill/altSS_forest_mycorrhizas_data/'}
if(host == 'Colins-MacBook-Pro-2.local'){storage.dir <- '/Users/colin/data_storage/altSS_forest_mycorrhizas_data/'}
if(host == 'Colins-MBP-2')              {storage.dir <- '/Users/colin/data_storage/altSS_forest_mycorrhizas_data/'}
cmd <- paste0('mkdir -p ',storage.dir)
system(cmd)

#FIA input paths.----
FIA7.dir.path <- '/fs/data3/caverill/FIA7/'
#FIA7.dir.path <- paste0(storage.dir,'raw_data/')
if(host == 'pecan2'){
  cmd <- paste0('mkdir -p ',FIA7.dir.path)
  system(cmd)
  FIAdb.path <- paste0(FIA7.dir.path,'FIA7.sqlite')
}

#FIA filtered output paths.----
fia.dir <- paste0(storage.dir,'FIA_output/')
cmd <- paste0('mkdir -p ',fia.dir)
system(cmd)
#All FIA data broken up by remeasurement.
all.present.path <- paste0(fia.dir,'FIA.all.present.rds')
all.past1.path   <- paste0(fia.dir,'FIA.all.past1.rds')
all.past2.path   <- paste0(fia.dir,'FIA.all.past2.rds')
all.past3.path   <- paste0(fia.dir,'FIA.all.past3.rds')

#FIA formatted analysis products.----
Product_1.path         <- paste0(fia.dir,"Product_1.rds")
Product_2.path         <- paste0(fia.dir,"Product_2.rds")
Product_2.subset.path  <- paste0(fia.dir,"Product_2.subset.rds")
time_series_dat.path   <- paste0(fia.dir,'time_series_dat.rds')

#GAM and simulation model output paths.----
model.dir <- paste0(storage.dir,'model_output/')
cmd <- paste0('mkdir -p ',model.dir)
system(cmd)

#GAM fits.
demographic_fits.path <- paste0(model.dir,'demographic_fits.rds')
#myco_gam_fits2.path <- paste0(model.dir,'myco_gam_fits2.rds')
 
#simulation output paths.
null_vs_feedback_simulation_output.path          <- paste0(model.dir,'null_vs_feedback_simulation_output.rds')
factorial_hysteresis_simulation.path             <- paste0(model.dir,'factorial_hysteresis_simulation.rds')
factorial_hysteresis_simulation_disturbance.path <- paste0(model.dir,'factorial_hysteresis_simulation_disturbance.rds')
factorial_hysteresis_simulation_uniform.path     <- paste0(model.dir,'factorial_hysteresis_simulation_uniform.path')

#Figure output paths.----
fig.dir <- 'figures/'
cmd <- paste0('mkdir -p ',fig.dir)
system(cmd)
Fig_1.path <- paste0(fig.dir,'Fig._1._conceptual_hysteresis.png')
Fig_2.path <- paste0(fig.dir,'Fig._2._disribution_and_bimodality.png')
Fig_3.path <- paste0(fig.dir,'Fig._3._recruitment_mortality_feedbacks.png')


