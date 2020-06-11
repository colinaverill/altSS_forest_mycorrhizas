#paths.r
#Main file dump directories.----
host <- system('hostname', intern = T)
storage.dir <- '/projectnb/talbot-lab-data/caverill/altSS_forest_mycorrhizas_data/'
if(host == 'pecan2')                    {storage.dir <- '/fs/data3/caverill/altSS_forest_mycorrhizas_data/'}
if(host == 'Colins-MacBook-Pro-2.local'){storage.dir <- '/Users/colin/data_storage/altSS_forest_mycorrhizas_data/'}
if(host == 'colins-MBP')                {storage.dir <- '/Users/colinaverill/Documents/data_storage/altSS_forest_mycorrhizas_data/'}
if(grepl('colins',host) == T)           {storage.dir <- '/Users/colinaverill/Documents/data_storage/altSS_forest_mycorrhizas_data/'}
#check if you're on ETH storage.dir
if(grepl('ethz.ch',host) == T)          {storage.dir <- '/Users/colinaverill/Documents/data_storage/altSS_forest_mycorrhizas_data/'}

cmd <- paste0('mkdir -p ',storage.dir)
system(cmd)

#FIA input paths.----
FIA7.dir.path <- '/fs/data3/caverill/FIA7/'
if(host == 'colins-MBP'      ){FIA7.dir.path <- '/Users/colinaverill/Documents/data_storage/FIA7/'}
if(grepl('ethz.ch',host) == T){FIA7.dir.path <- '/Users/colinaverill/Documents/data_storage/FIA7/'}
#cmd <- paste0('mkdir -p ',FIA7.dir.path)
#system(cmd)
FIAdb.path <- paste0(FIA7.dir.path,'FIA7.sqlite')

#Other raster data outside of main storage directory.----
EPA_L2_ecoregions_raster.path <- '/Users/colinaverill/Documents/data_storage/misc_rasters/NA_level2_ecoregions/NA_CEC_Eco_Level2.shp'

#other (small) database paths.----
nodDB.path <- 'required_products_utilities/nodDB_v1.csv'

#FIA filtered output paths.----
fia.dir <- paste0(storage.dir,'FIA_output/')
cmd <- paste0('mkdir -p ',fia.dir)
system(cmd)
#All FIA data broken up by remeasurement.
all.present.path <- paste0(fia.dir,'FIA.all.present.rds')
all.past1.path   <- paste0(fia.dir,'FIA.all.past1.rds')
all.past2.path   <- paste0(fia.dir,'FIA.all.past2.rds')
all.past3.path   <- paste0(fia.dir,'FIA.all.past3.rds')

#paths for grabbing and returning lab composite environmental data.----
 data_for_composite.path <- paste0(storage.dir, 'data_for_composite_CA.csv')
data_from_composite.path <- paste0(storage.dir,'data_from_composite_CA.csv')

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
rf_demographic_fits.path <- paste0(model.dir,'rf_demographic_fits.rds')
rf_demographic_fits_interactive.path <- paste0(model.dir,'rf_demographic_fits_interactive.rds')
demographic_fits_gam_separate.path <- paste0(model.dir,'demographic_fits_gam_separate.rds')
demographic_fits_gam_ecoregion.path <- paste0(model.dir,'demographic_fits_gam_ecoregion.rds')
  demographic_fits_gam_species.path <- paste0(model.dir,'demographic_fits_gam_species.rds')
  demographic_fits_gam_separate_plus_re.path <- paste0(model.dir,'demographic_fits_gam_separate_plus_re.rds')

#Demographic simulation output paths.
        null_vs_feedback_simulation_output.path <- paste0(model.dir,'null_vs_feedback_simulation_output.rds')
           factorial_hysteresis_simulation.path <- paste0(model.dir,'factorial_hysteresis_simulation.rds')
   initial_condition_hysteresis_simulation.path <- paste0(model.dir,'initial_condition_hysteresis_simulation.rds')
     RF_null_vs_feedback_simulation_output.path <- paste0(model.dir,'RF_null_vs_feedback_simulation_output.rds')
   
#Figure output paths.----
fig.dir <- 'figures/'
cmd <- paste0('mkdir -p ',fig.dir)
system(cmd)
      Fig_1.path <- paste0(fig.dir,'Fig._1._distribution_and_bimodality.png')
      Fig_2.path <- paste0(fig.dir,'Fig._2._recruitment_mortality_feedbacks.png')
      Fig_3.path <- paste0(fig.dir,'Fig._3._demographic_simulations.png')
      Fig_4.path <- paste0(fig.dir,'Fig._4._hysteresis.png')
Supp._gif_1.path <- paste0(fig.dir,'Supp._gif_1.gif')

