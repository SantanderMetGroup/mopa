
#' @title Delimit study area and background 
#' @description Creation of polygon shapes from bounding coordinates and delimitation of 
#' SpatialPoints data to the defined boundaries
#' 
#' @param bounding.coords Object returned by function \code{\link{boundingCoords}}. 
#' Matrix or list of matrixes of bounding coordinates indicating the maximum and minimum 
#' values in columns and xy coordinates in rows
#' @param grid Projected SpatialPoints object
#' @param names Names or IDs to be given to each shape. Default bounding.coords names.
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
#' @author M. Iturbide \email{maibide@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' data(Oak_phylo2)
#' data(biostack)
#' ##creation of point grid from raster object
#' library(raster)
#' ac<-xyFromCell(biostack[[1]],  1:ncell(biostack[[1]]))
#' ex<-extract(biostack[[1]], ac)
#' ##exclude sea and project
#' sp_grid<-SpatialPoints(ac[-which(is.na(ex)),])
#' projection(sp_grid)<-CRS("+proj=longlat +init=epsg:4326")
#' plot(sp_grid)
#' ##delimit study area
#' oak.extension<-boundingCoords(Oak_phylo2)
#' box.grid <- delimit(bounding.coords = oak.extension, grid = sp_grid, names = names(Oak_phylo2))
#'
#' ## Plot presences and bounding boxes
#' plot(box.grid$bbs, asp = 1)
#' for (i in 1:length(Oak_phylo2)){
#'  points(Oak_phylo2[[i]], col = colors()[i*50])
#'  }
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
  
  if (class(bounding.coords) != "list"){
    bounds<-list(bounding.coords)
    names(bounds) <- "background"
  }else{bounds<-bounding.coords}
  
  bb0<-bounds
  bb1 <- lapply(bb0, bboxx)
  bb1[[1]] 
  bb2 <- lapply(bb1, function(x) rbind(x, x[1,]))
  bb2[[1]]
  
  if(!is.null(names)){
  rn <- names
  }else{
    rn <- names(bounds)
  }

  bb3 <- vector(mode="list", length=length(bb2))

    for (i in seq(along=bb3)) bb3[[i]] <- Polygons(list(Polygon(bb2[[i]])),
                                                 ID=rn[i])
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
  
  if (length(bbs.grid) == 1){bbs.grid<-bbs.grid[[1]]
  }else{names(bbs.grid)<-rn} 
  
  return(list("bbs"=bbs, "bbs.grid"=bbs.grid))
}


