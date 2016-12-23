
#' @title Delimit study area and background coordinates
#' @description Creation of polygon shapes from bounding coordinates and delimitation of 
#' SpatialPoints data to the defined boundaries
#' 
#' @param bounding.coords A vector or a list of vectors with bounding coordinates in the following form: c(x1, x2, y1, y2). 
#' Also object returned by function \code{\link{boundingCoords}}. 
#' @param grid Projected SpatialPoints object
#' @param names Character. Names or IDs to be given to each shape. If not specified names of the object passed to 
#' bounding.coords will be given. In case this object has no names the background numbers will be used as names. 
#' 
#' @return A list with two components: 
#' \itemize{
#'  \item{bbs}{SpatialPolygons of the bounding boxes} 
#'  \item{bbs.grid}{list(s) of matrix(ces) of the background xy coordinates in columns}.
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
#' @examples
#' \dontrun{
#' data(Oak_phylo2)
#' data(biostackENSEMBLES)
#' presences <- Oak_phylo2
#' 
#' ##creation of point grid from raster object
#' sp_grid <- background(biostackENSEMBLES$baseline$bio2)
#' 
#' ##delimit study area to the whole study domain for both species
#' bc <- rep(list(boundingCoords(coordinates(sp_grid))), length(presences))
#' del <- delimit(bounding.coords = bc, grid = sp_grid)
#' ## Plot presences and bounding boxes
#' plot(del$bbs, asp = 1)
#' for (i in 1:length(presences)){
#'   points(presences[[i]], col = colors()[i*50])
#' }
#' 
#' ##delimit study area to the bounding coordinates of each species
#' bc <- boundingCoords(presences)
#' del <- delimit(bounding.coords = bc, grid = sp_grid)
#' 
#' ## Plot presences and bounding boxes
#' plot(del$bbs, asp = 1)
#' for (i in 1:length(presences)){
#'   points(presences[[i]], col = colors()[i*50])
#'   }
#'}
#'
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @export
#' 
#' @import sp
#' @importFrom splancs bboxx

delimit<-function(bounding.coords, grid, names = NULL){
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
  return(list("bbs"=bbs, "bbs.grid"=bbs.grid))
}


