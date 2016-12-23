

#' @title Extract spatial background from raster
#' @description Obtain a spatial points object of the background given by a raster object.
#' 
#' @param raster Raster object (typically corresponding to a variable used for modelling species distributions).
#' @param projection Object of class CRS with the coordinate reference system. Default is 
#' CRS("+proj=longlat +init=epsg:4326") 
#' 
#'  
#' 
#' @return SpatialPoints object.
#' 
#' 
#' @author M. Iturbide 
#' 
#' @examples
#' 
#' \dontrun{
#' data(biostackENSEMBLES)
#' grid <- background(biostackENSEMBLES$baseline$bio2) 
#' plot(grid)
#' }
#' @export
 
background <- function(raster, projection = CRS("+proj=longlat +init=epsg:4326")){
ac <- xyFromCell(raster, 1:ncell(raster))
ex <- extract(raster, ac)

# Convert to a Spatial object and define projection
sp_grid <- SpatialPoints(ac[-which(is.na(ex)), ])
projection(sp_grid) <- projection
return(sp_grid)
}

#end
