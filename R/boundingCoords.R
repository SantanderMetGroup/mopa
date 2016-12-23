
#' @title Bounding box coordinates of xy records
#' @description Creates a vector or a list of vectors of bounding coordinates around point records (xy records)
#' @param xy Data frame or list of data frames with coordinates (each row is a point)
#' 
#' @return A vector or a list of vectors with bounding coordinates in the following form: c(x1, x2, y1, y2). 
#' 
#' @author M. Iturbide 
#' 
#' @examples 
#' \dontrun{
#' data(Oak_phylo2)
#' oak.bounds <- boundingCoords(xy = Oak_phylo2)
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
  if (length(l.box) == 1){
    box<-l.box[[1]]
  }else{
    box<-l.box
  }
  return(box)
}