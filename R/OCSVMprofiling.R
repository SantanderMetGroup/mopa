
#' @title Environmental profiling with One-Classification Support Vector Machine
#' @description Presence-only modeling and classification of coordinates predicted as 
#' presence and absence
#' 
#' @param xy Data frame or list of data frames with coordinates (each row is a point)
#' @param varstack RasterStack of variables for modeling
#' @param background Object derived from function \code{\link[mopa]{backgroundGrid}}. If NULL (default),
#' the background is extracted from varstack.
#' Matrix or list of matrixes of the background xy coordinates in columns. 
#' @param nu Parameter needed for one-classification \code{\link[e1071]{svm}}. 
#' Default is 0.5
#'
#' @seealso \code{\link[e1071]{svm}}, \code{\link[mopa]{backgroundGrid}}
#' 
#' @return  A list with two components:
#'   \item{absence }{Matrix or list of matrixes with xy coordinates predicted 
#' as absence}
#' \item{presence }{Matrix or list of matrixes with xy coordinates predicted as presence}
#' 
#' @details This function constitutes the 
#' first step from a three-step proccess to generate pseudo-absences, and is aimed at 
#' excluding the suitable areas for the species (xy records) from the background 
#' from which pseudo-absences are sampled. 
#' 
#' 
#' 
#' @author M. Iturbide 
#' 
#' @examples
#' data(Oak_phylo2)
#' data(biostack)
#' projection(biostack$baseline) <- CRS("+proj=longlat +init=epsg:4326")
#' r <- biostack$baseline[[1]]
# Background of the whole study area
#' bg <- backgroundGrid(r)
## environmental profiling
#' bg.profiled <- OCSVMprofiling(xy = Oak_phylo2, varstack = biostack$baseline, 
#'                               background = bg$xy)
#' ##plot
#' plot(bg.profiled$absence$H11, pch="*")
#' points(bg.profiled$presence$H11, pch="*", col= "pink")
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @export
#' 
#' @importFrom e1071 svm
 

OCSVMprofiling<-function(xy, varstack, background = NULL, nu=.5){
  if (class(xy) != "list")  xy <- list(xy)
  if(is.null(background)) background <- backgroundGrid(varstack[[1]])$xy
  if (class(background) != "list")  background <- rep(list(background), length(xy))
  
  bioclim <-varstack
  absence <- list()    
  presence <- list()
  for(i in 1:length(xy)){
    length(absence) <- i
    coo <- background[[i]]
    mat <- cbind(xy[[i]], rep(1, nrow(xy[[i]])))
    mat <- biomat(mat, bioclim)
    mod <- svm(mat[,-1], y=NULL, type='one-classification', nu=nu)
    proj <- biomat(cbind(coo,rep(1,nrow(coo))), bioclim)
    pre <- predict(mod, proj[,-1])
    absence[[i]] <- coo[(which(pre==0)),]
    presence[[i]] <- coo[(which(pre!=0)),]
  }
  if(length(absence) == 1){
    absence <- absence[[1]]
    presence <- presence[[1]]
  }else{
    names(absence) <- names(xy)
    names(presence) <- names(xy)
  }
  return(list("absence"=absence, "presence"=presence))
}
