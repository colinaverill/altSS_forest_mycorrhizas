#' Extract worldclim2 temperature, precipitation and uncertainties at 30s resolution for a given latitude/longtiude.
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
  
  #make points an object
  points <- cbind(longitude, latitude)
 
  #load mean annual temperature and precipitation rasters from worldclim2
  prec    <- raster::raster(paste0(worldclim2_folder,'wc2.0_bio_30s_12.tif'))
  temp    <- raster::raster(paste0(worldclim2_folder,'wc2.0_bio_30s_01.tif'))
  temp_CV <- raster::raster(paste0(worldclim2_folder,'wc2.0_bio_30s_04.tif'))
  prec_CV <- raster::raster(paste0(worldclim2_folder,'wc2.0_bio_30s_15.tif'))
  mdr     <- raster::raster(paste0(worldclim2_folder,'wc2.0_bio_30s_02.tif'))
  
  
  #extract worldclim2 predicted climate data.
     prec.obs <- raster::extract(prec   , points)
     temp.obs <- raster::extract(temp   , points)
  prec_CV.obs <- raster::extract(prec_CV, points)
  temp_CV.obs <- raster::extract(temp_CV, points)
      mdr.obs <- raster::extract(    mdr, points)
  
  #wrap output to return.
  to_return <- data.frame(cbind(prec.obs,temp.obs,prec_CV.obs,temp_CV.obs,mdr.obs))
  colnames(to_return) <- c('map','mat','map_CV','mat_CV','mdr')
      
  #return output.
  return(to_return)
  
} #end function.
