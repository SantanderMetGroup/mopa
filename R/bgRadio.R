
#' @title Background extent restriction for a sequence of distances 
#' @description Creation of point-grid backgrounds through the establishment of extent 
#' limitations for a sequence of distances, with the length of the half diagonal of the 
#' bounding box around xy records as the maximum distance
#' 
#' @param xy Data frame or list of data frames with coordinates (each row is a point)
#' @param bounding.coords Object returned by function \code{\link[mopa]{boundingCoords}}.
#' Matrix or list of matrixes of bounding coordinates indicating the maximum and minimum 
#' values in columns and xy coordinates in rows
#' @param bg.absence Object derived from function \code{\link[mopa]{OCSVMprofiling}} ($absence). 
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
#' @author M. Iturbide \email{maibide@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' ##delimit study area
#' data(Oak_phylo2)
#' data(sp_grid)
#' oak.extension<-boundingCoords(Oak_phylo2)
#' box.grid<-delimit(oak.extension, sp_grid, names(Oak_phylo2))
#' ## environmental profiling
#' data(biostack)
#' unsuitable.bg <-OCSVMprofiling(xy = Oak_phylo2, varstack = biostack, 
#' bbs.grid = box.grid$bbs.grid)
#' ## sequence of 10 km between distances, from 20 km to the length of the 
#' ##half diagonal of the bounding box.
#' ext <- bgRadio(xy = Oak_phylo2, bounding.coords = oak.extension, 
#' bg.absence = unsuitable.bg$absence, start = 0.166, 
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



bgRadio<-function(xy, bounding.coords, bg.absence, start= 0.166, by= 0.083, 
                  unit= "decimal degrees"){
  
  bg1a<-list()
  if(class(xy) == "matrix") xy <- as.data.frame(xy)
  if (class(xy) == "data.frame"){
    pres<- list(xy)
  }else{
    pres<-xy
  }
  
  if (class(bounding.coords) != "list"){
    l.box<-rep(list(bounding.coords),length(xy))
  }else{l.box<-bounding.coords}
  
  if (class(bg.absence) != "list"){
    bg.absence<-list(bg.absence)
  }else{bg.absence<-bg.absence}
  
  bg.rad<-function(pres, bg.absence, r){
    bg.absence[unique(nearest.dist(x=pres, y=bg.absence, 
                                   method="maximum", delta=r)@colindices),]
  }
  
    for (i in 1:length(pres)){
      print(paste("creating background point-grids for species", i, "out of", length(pres)))
      box<-l.box[[i]]
      pr<-pres[[i]]
      gbox<-bg.absence[[i]]
    
      a<-(max(l.box[[i]][1,])-min(l.box[[i]][1,]))/2
      b<-(max(l.box[[i]][2,])-min(l.box[[i]][2,]))/2
      c<-sqrt(a*a+b*b)
      radios<-seq(start, c, by)
      r<-radios
      bg1a[[i]]<-sapply(r, FUN=bg.rad, pres=pr, bg.absence=gbox)
      
      if (unit=="decimal degrees"){
        names(bg1a[[i]])<-paste("km", as.character(r/0.0083), sep="")
      } else if (unit=="utm"){
        names(bg1a[[i]])<-paste("km", as.character(r/1000), sep="")
      } else {names(bg1a[[i]])<-NULL}
    }
  
  if (length(bg1a) == 1){bg1a <- bg1a[[1]]
  }else{names(bg1a)<-names(pres)}
  
  return("bg"=bg1a)
}
