#' Extract the following worldclim2 bioclim variables at 30s resolution for a given latitude/longtiude.
#' - mean annual temperature (BIO01)
#' - mean annual precipitation (BIO12)
#' - temperature seasonality (BIO04)
#' - temperature range (BIO07)
#' - temperature wettest quarter (BIO08)
#' - temperature warmest quarter (BIO10)
#' - precipitation wettest quarter (BIO16)
#' - precipitation warmest quarter (BIO18)
#' The precipitation uncertainty scales with elevation. If the user does not supply elevation then we default to 500m (mean in worldclim2 calibration data).
#' Would be nice to add a topography product to this when elevation is not supplied.
#' Only precipitation uncertainty scales with elevation.
#' Depends on output from JAGS models Colin fit. Not totally sure where to store these.
#' 
#'
#' @param latitude         vector of latitude in decimal degrees.
#' @param longitude        vector of longitude in decimal degrees.
#' @param worldclim2_foder path to directory of worldclim2 rasters.
#'
#' @return           returns a datafrme with mat, map, mat_CV, map_CV and mdr
#' @export
#'
#' @examples
#' points <- structure(c(-148.274862, -148.274862, -148.2747566, -148.2747566, 
#'                      -148.2746513, -148.2746513, -148.2744406, -148.2744406, -148.2740191, 
#'                      -148.2740191, 64.766184, 64.766184, 64.766184, 64.766184, 64.766184, 
#'                      64.766184, 64.766184, 64.766184, 64.766184, 64.766184), .Dim = c(10L, 
#'                                                                                      2L))
#' test.out <- worldclim2_grab(points[,2], points[,1])
worldclim2_grab <- function(latitude, longitude,
                            worldclim2_folder = '/projectnb/talbot-lab-data/caverill/WorldClim2/'
                            ){
  #Figure out where dat are stored if not on BU SCC.
  host <- system('hostname', intern = T)
  if(host == 'Colins-MacBook-Pro-2.local'){worldclim2_folder <- '/Users/colin/data_storage/WorldClim2/'}
  if(host == 'Colins-MBP-2')              {worldclim2_folder <- '/Users/colin/data_storage/WorldClim2/'}
  
  #make points an object
  points <- cbind(longitude, latitude)
 
  #load mean annual temperature and precipitation rasters from worldclim2
  mat      <- raster::raster(paste0(worldclim2_folder,'wc2.0_bio_30s_12.tif'))
  t_seas   <- raster::raster(paste0(worldclim2_folder,'wc2.0_bio_30s_04.tif'))
  t_range  <- raster::raster(paste0(worldclim2_folder,'wc2.0_bio_30s_07.tif'))
  t_wet_q  <- raster::raster(paste0(worldclim2_folder,'wc2.0_bio_30s_08.tif'))
  t_warm_q <- raster::raster(paste0(worldclim2_folder,'wc2.0_bio_30s_10.tif'))
  map      <- raster::raster(paste0(worldclim2_folder,'wc2.0_bio_30s_01.tif'))
  p_wet_q  <- raster::raster(paste0(worldclim2_folder,'wc2.0_bio_30s_16.tif'))
  p_warm_q <- raster::raster(paste0(worldclim2_folder,'wc2.0_bio_30s_18.tif'))
  
  #extract worldclim2 predicted climate data.
          mat.obs <- raster::extract(mat     , points)
       t_seas.obs <- raster::extract(t_seas  , points)
      t_range.obs <- raster::extract(t_range , points)
      t_wet_q.obs <- raster::extract(t_wet_q , points)
     t_warm_q.obs <- raster::extract(t_warm_q, points)
          map.obs <- raster::extract(map     , points)
      p_wet_q.obs <- raster::extract(p_wet_q , points)
     p_warm_q.obs <- raster::extract(p_warm_q, points)

  #wrap output to return.
  to_return <- data.frame(
                          cbind(mat.obs,
                                t_seas.obs,
                                t_range.obs,
                                t_wet_q.obs,
                                t_warm_q.obs,
                                map.obs,
                                p_wet_q.obs,
                                p_warm_q.obs)
                          )
  colnames(to_return) <- c('mat','t_seas','t_range','t_wet_q','t_warm_q','map','p_wet_q','p_warm_q')
      
  #return output.
  return(to_return)
  
} #end function.
