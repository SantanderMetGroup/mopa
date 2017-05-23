#' @title Background extent restriction for a sequence of distances 
#' @description Creation of point-grid backgrounds through the establishment of extent 
#' limitations for a sequence of distances, from near presence locations to the length 
#' of the half diagonal of the bounding that encloses the background (study area).
#' 
#' @param xy Data frame or list of data frames with coordinates (each row is a point) 
#' ---typically species presence data--- to be considered as starting points from which different 
#' background extents are created.
#' @param background Matrix or list of matrices of background coordinates.
#' Object derived from function \code{\link[mopa]{OCSVMprofiling}} (component \code{$absence}). 
#' Alternatively, object derived from function \code{\link[mopa]{backgroundGrid}} (component \code{$xy}) 
#' if the environmental profiling step is going to be avoided in the pseudo-absence 
#' generation proccess). 
#' @param start Value for the minimum distance to consider for extent limitations. 
#' Default is 0.166
#' @param by Value of the distance to consider from one extent to the following. 
#' Default is 0.083
#' @param unit Character indicating the coordinate system of the objects. 
#' Default is \code{"decimal degrees"}, alternatively \code{"utm"} can be used
#' 
#' @return List/s of matrixes with xy coordinates, each matrix correspond to a different 
#' background extent.
#' 
#' 
#' @details Argument unit is only used to set extent distances in km. This function is aimed at 
#' creating backgrounds of different extent for pseudo-absence samplig based on an initial 
#' point grid (derived from function \code{\link[mopa]{OCSVMprofiling}} or function 
#' \code{\link[mopa]{backgroundGrid}}). If this function is used for a subsequent application of
#' functions \code{\link[mopa]{pseudoAbsences}} and \code{\link[mopa]{mopaTrain}}, the last will perform species distribution
#' modeling for each of the extents here established, and will return the fitted model that belongs to 
#' the optimum background extent (see references). 
#' 
#' @seealso \code{\link[mopa]{mopaTrain}}, \code{\link[mopa]{pseudoAbsences}}, \code{\link[mopa]{backgroundGrid}}, 
#' \code{\link[mopa]{OCSVMprofiling}}
#' 
#' @author M. Iturbide 
#' 
#' @examples
#' data(Oak_phylo2)
#' data(biostack)
#' projection(biostack$baseline) <- CRS("+proj=longlat +init=epsg:4326")
#' r <- biostack$baseline[[1]]
#' # Background of the whole study area
#' bg <- backgroundGrid(r)
#' bg.extents <- backgroundRadius(xy = Oak_phylo2, background = bg$xy, 
#' start = 0.166, by = 0.083*10, unit = "decimal degrees")
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @export
#' 
#' @importFrom spam nearest.dist
#' 



backgroundRadius <- function(xy, background, start= 0.166, by= 0.083, 
                    unit = c("decimal degrees", "utm")){
  unit <- match.arg(unit, choices = c("decimal degrees", "utm"))
  if(class(xy) == "matrix") xy <- as.data.frame(xy)
  if(class(xy) == "data.frame") xy <- list(xy)
  if(class(background) != "list") background <- rep(list(background),length(xy))
  bounding.coords <- boundingCoords(background)
  bg.rad <- function(xy, background, r){
    background[unique(nearest.dist(x = xy, 
                                   y = background, 
                                   method = "maximum", 
                                   delta = r)@colindices),]
  }
  bg1a <- list()
  for (i in 1:length(xy)){
    print(paste("creating background point-grids for species", i, "out of", length(xy)))
    box <- matrix(bounding.coords[[i]], ncol = 2)
    pr <- xy[[i]]
    gbox <- background[[i]]
    
    a <- (max(box[1,])-min(box[1,]))/2
    b <- (max(box[2,])-min(box[2,]))/2
    c <- sqrt(a*a+b*b)
    radios <- seq(start, c, by)
    r <- radios
    bg1a[[i]] <- sapply(r, FUN=bg.rad, xy = pr, background = gbox)
    
    if(unit == "decimal degrees"){
      names(bg1a[[i]]) <- paste("km", as.character(r/0.0083), sep="")
    } else if (unit == "utm"){
      names(bg1a[[i]]) <- paste("km", as.character(r/1000), sep="")
    } else {
      names(bg1a[[i]]) <- NULL
    }
  }
  # if (length(bg1a) == 1){
  #   bg1a <- bg1a[[1]]
  # }else{
    names(bg1a) <- names(xy)
  # }
  return("bg"= bg1a)
}

#end
