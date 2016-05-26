
#' @title Bounding box coordinates of xy records
#' @description Creates the matrix of bounding coordinates around point records (xy records)
#' @param xy Data frame or list of data frames with coordinates (each row is a point)
#' 
#' @return Matrix or list of matrixes with the bounding coordinates.
#' 
#' @author M. Iturbide \email{maibide@@gmail.com}
#' 
#' @examples 
#' \dontrun{
#' data(Oak_phylo2)
#' oak.extension <- boundingCoords(xy = Oak_phylo2)
#' }
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' @export


boundingCoords<-function(xy){
  
  if (class(xy) == "matrix") xy <- as.data.frame(xy)
  if (class(xy) == "data.frame"){
    pres.list<-list(xy)
  }else{pres.list<-xy}
  
  box<-matrix(NA,2,2,dimnames=list(c("x","y"),c("min","max")))
  l.box<-list()
  
    for (i in 1:length(pres.list)){
      length(l.box)<-i
      box[1,1]<-min(pres.list[[i]][,1])
      box[2,1]<-min(pres.list[[i]][,2])
      box[1,2]<-max(pres.list[[i]][,1])
      box[2,2]<-max(pres.list[[i]][,2])
      l.box[[i]]<-box
    }
  
  names(l.box)<-names(pres.list)
  
  if (class(xy)=="data.frame"){
    box<-l.box[[1]]
  }else{
    box<-l.box
  }
  
  return(box)
}
