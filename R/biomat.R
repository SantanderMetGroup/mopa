#' @title Matrix with variables for modelling
#' @description Prepares matrix with variables for modeling
#' 
#' @param data Data frame with coordinates in the first two columns and presence/absence 
#' (1=presence, 0=absence) in the third column.
#' @param varstack RasterStack of the variables from which values are extracted for each point in
#' \code{data}. 
#' 
#' @return  2D matrix with the dependent variable (presence/absence) in the first column
#' and the independent variables (extracted from varstack) in the rest.
#' 
#' 
#' @author M. Iturbide 
#' 
#' @examples
#' data(biostack)
#' data(Oak_phylo2)
#' dfp <-cbind(Oak_phylo2[[1]], "pa"= rep(1,nrow(Oak_phylo2[[1]])))
#' dfa <-cbind(Oak_phylo2[[2]], "pa"= rep(0,nrow(Oak_phylo2[[2]])))
#' df3 <-rbind(dfp, dfa)
#' mat <-biomat(df3, biostack$baseline)
#' str(mat)
#' 
#' @keywords internal
#' @import sp
#' @import raster
#' @export

biomat <- function(data, varstack){ 
  bio <- varstack
  coord.esp <- data[,1:2]
  sp.coord.esp <- SpatialPoints(coord.esp)
  # proj4string(sp.coord.esp)<-projection
  # proj4string(bio)<-projection
  z <- extract(bio,sp.coord.esp)
  if(ncol(data) == 3){
    bio.mat <- cbind(data[,3],z)
  }else{
    bio.mat <- z
  }
  bio.df <- as.data.frame(bio.mat)
  return(bio.df)
}

#end
