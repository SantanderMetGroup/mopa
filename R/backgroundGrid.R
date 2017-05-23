
#' @title Create background coordinates from raster object
#' @description Creates the background coordinates used for defining the study area and for generating pseudo-absences.
#' @param raster Raster object with projection (\code{?raster::projection} and \code{?crs}) 
#' from which to extract the the point grid (used as background coordinates for generating pseudo-absences)
#' @param spatial.subset Object of class extent (see \code{\link[raster]{extent}}) or a two
#' column data.frame (or matrix) of coordinates (xy, each row is a point).
#' 
#' @details If a data.frame, matrix or a list of the previous is passed to \code{spatial.subset} the bounding coordinates
#' are extracted to delimit the background. For example, to bound the study area to the spatial distribution of a species.
#' 
#' @return A list with a SpatialPolygons object and a matrix of the background coordinates 
#' 
#' @author M. Iturbide 
#' 
#' @importFrom abind abind
#' 
#' @examples 
#' data(Oak_phylo2)
#' data(biostack)
#' projection(biostack$baseline) <- CRS("+proj=longlat +init=epsg:4326")
#' r <- biostack$baseline[[1]]
#' # Background around a set of coordinates
#' bg.species <- backgroundGrid(r, Oak_phylo2)
#' # Background of a subdomain of the study area
#' bg.subdomain <- backgroundGrid(r, extent(c(-10, 30, 35, 65)))
#' # Background of the whole study area
#' bg <- backgroundGrid(r)
#' # plot(bg$xy)
#' # plot(bg.subdomain$xy)
#' # plot(bg.species$xy$H11)
#' plot(bg.species$xy$H01)


#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' @export


backgroundGrid <- function(raster, spatial.subset = NULL){
  #if(is.na(projection(raster))) stop("raster has no projection arguments. Type ?raster::projection")
  if(is.null(spatial.subset)) spatial.subset <- extent(raster)
  if(class(spatial.subset) == "Extent"){
      raster <- crop(raster, spatial.subset)
  }else{
    a <- abind(boundingCoords(spatial.subset), along = 2)
    b <- c(min(a[1,]), max(a[2,]), min(a[3,]), max(a[4,]))
    raster <- crop(raster, extent(b))
  }
  ac <- xyFromCell(raster, 1:ncell(raster))
  ex <- extract(raster, ac)
  # Convert to a Spatial object and define projection
  sp_grid <- SpatialPoints(ac[-which(is.na(ex)), ])
  projection(sp_grid) <- projection(raster)
  if(class(spatial.subset) != "Extent"){
      bcoord <- boundingCoords(spatial.subset)
  }else{
      bcoord <- boundingCoords(coordinates(sp_grid))
  }
  del <-delimit(bounding.coords = bcoord, grid = sp_grid)
  return(del)
}
  
#end





#' @title Bounding box coordinates of xy records
#' @description Creates a vector or a list of vectors of bounding coordinates around point records (xy records).
#' This is an internal function used by function
#' \code{\link{backgroundGrid}}. 
#' @param xy Data frame or list of data frames with coordinates (each row is a point)
#' 
#' @return A vector or a list of vectors with bounding coordinates 
#' in the following form: c(x1, x2, y1, y2). 
#' 
#' @author M. Iturbide 
#' 
#' 
#' @keywords internal
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.



boundingCoords<-function(xy){
  if (class(xy) == "matrix") xy <- as.data.frame(xy)
  if (class(xy) == "data.frame"){
    pres.list<-list(xy)
  }else{pres.list<-xy}
  box <- numeric()
  l.box<-list()
  for (i in 1:length(pres.list)){
    length(l.box)<-i
    box[1]<-min(pres.list[[i]][,1])
    box[2]<-max(pres.list[[i]][,1])
    box[3]<-min(pres.list[[i]][,2])
    box[4]<-max(pres.list[[i]][,2])
    l.box[[i]]<-box
  }
  names(l.box)<-names(pres.list)
  return(l.box)
}

#end




#' @title Delimit study area and background coordinates
#' @description Creation of polygon shapes from bounding coordinates and delimitation of 
#' SpatialPoints data to the defined boundaries. This is an internal function used by function
#' \code{\link{backgroundGrid}}. 
#' 
#' @param bounding.coords A vector or a list of vectors with bounding coordinates in the following form: c(x1, x2, y1, y2). 
#' Also object returned by function \code{\link{boundingCoords}}. 
#' @param grid SpatialPoints object
#' @param names Character. Names or IDs to be given to each shape. If not specified names of the object passed to 
#' bounding.coords will be given. In case this object has no names the background numbers will be used as names. 
#' 
#' @return A list with two components: 
#' \itemize{
#'  \item{polygons}{SpatialPolygons of the bounding boxes} 
#'  \item{xy}{list(s) of matrix(ces) of the background xy coordinates in columns}.
#' }
#' 
#' @details This function is aimed at restricting the study area inside the bounding boxes 
#' of a group of points (presences). The study area can be represented as a SpatialPoints 
#' object of a regular point grid, this grid also represents the background (excluding 
#' presences) from which pseudo-absences can be sampled. 
#' Ideally, this grid is extracted from the raster objects of the variables (centroids)
#' to be used for modelling, see \code{\link[raster]{xyFromCell}}, 
#' \code{\link[raster]{extract}} and example below. bbs.grid is the background inside 
#' bbs shapes.
#' 
#' 
#' 
#' @author M. Iturbide 
#' 
#'
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @keywords internal
#' @import sp
#' @importFrom splancs bboxx

delimit <- function(bounding.coords, grid, names = NULL){
  if (class(bounding.coords) != "list") bounding.coords <- list(bounding.coords)
  boundsmat <- lapply(1:length(bounding.coords), function(x){
    matrix(bounding.coords[[x]], ncol = 2, byrow = T)
  })
  bb0 <- boundsmat
  bb1 <- lapply(bb0, bboxx)
  bb1[[1]] 
  bb2 <- lapply(bb1, function(x) rbind(x, x[1,]))
  bb2[[1]]
  if(!is.null(names)){
    rn <- names
  }else if(is.null(names(bounding.coords))){
    rn <- as.character(paste("background_", 1:length(bounding.coords), sep = ""))
  }else{
    rn <- names(bounding.coords)
  }
  bb3 <- vector(mode="list", length=length(bb2))
  for (i in seq(along=bb3)){
    bb3[[i]] <- Polygons(list(Polygon(bb2[[i]])), ID=rn[i])
  } 
  bb3[[1]]
  bbs <- SpatialPolygons(bb3, proj4string=CRS(projection(grid)))
  bbs.grid<-list()
  projection(bbs)<-projection(grid)
  co<-coordinates(grid)
  
  for (i in 1:length(bbs)){
    length(bbs.grid)<-i
    lo<-over(grid,bbs[i]) 
    oo<-co[which(lo==1),] 
    bbs.grid[[i]]<-oo
  }
  
  if (length(bbs.grid) == 1){bbs.grid <- bbs.grid[[1]]
  }else{
    names(bbs.grid) <- rn
  } 
  return(list("polygons" = bbs, "xy" = bbs.grid))
}


#end

 