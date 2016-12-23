#' @title Background extent restriction for a sequence of distances 
#' @description Creation of point-grid backgrounds through the establishment of extent 
#' limitations for a sequence of distances, with the length of the half diagonal of the 
#' bounding box around xy records as the maximum distance
#' 
#' @param xy Data frame or list of data frames with coordinates (each row is a point)
#' @param background Matrix or list of matrixes of background coordinates.
#' Object derived from function \code{\link[mopa]{OCSVMprofiling}} ($absence). 
#' Alternatively, object derived from function \code{\link[mopa]{delimit}} ($bbs.grid) 
#' if the environmental profiling step is going to be avoided in the pseudo-absence 
#' generation proccess). 
#' @param start Value for the minimum distance to consider for extent limitations. 
#' Default is 0.166
#' @param by Value of the distance to consider from one extent to the following. 
#' Default is 0.083
#' @param unit Character indicating the coordinate system of the objects. 
#' Default is "Decimal degrees", alternatively "utm" can be used
#' 
#' @return List/s of matrixes with xy coordinates, each matrix correspond to a different 
#' background extent.
#' 
#' 
#' @details Argument unit is only used to set extent names. This function forms part of the 
#' second step in a three-step proccess to generate pseudo-absences, and is aimed at 
#' creating backgrounds of different extent for pseudo-absence samplig based on an initial 
#' point grid (derived from function \code{\link[mopa]{OCSVMprofiling}} or function 
#' \code{\link[mopa]{delimit}}). 
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
#' del <- delimit(bounding.coords = bc, grid = sp_grid, names = names(presences))
#' ## environmental profiling
#' unsuitable.bg <- OCSVMprofiling(xy = presences, varstack = biostackENSEMBLES$baseline, 
#' bbs.grid = del$bbs.grid)
#' ## sequence of 10 km between distances, from 20 km to the length of the 
#' ##half diagonal of the bounding box.
#' ext <- bgRadio(xy = presences,
#' background = unsuitable.bg$absence, start = 0.166, 
#' by = 0.083, unit = "decimal degrees")
#' 
#' ext <- bgRadio(xy = presences,
#' background = del$bbs.grid, start = 0.166, 
#' by = 0.083, unit = "decimal degrees")
#' 
#' # Plot presences for group H11 and background extents of 20, 120 and 520 km
#' plot(ext$H11$km520, col = "green4", pch = "*", asp = 1)
#' points(ext$H11$km120, pch = "*")
#' points(ext$H11$km20, pch = "*", col = "blue")
#' points(Oak_phylo2$H11, col = "red", pch = ".", cex = 1.5)
#' }
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @export
#' 
#' @importFrom spam nearest.dist
#' 



bgRadio <- function(xy, background, start= 0.166, by= 0.083, 
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
  if (length(bg1a) == 1){
    bg1a <- bg1a[[1]]
  }else{
    names(bg1a) <- names(xy)
  }
  return("bg"= bg1a)
}

#end
