
#' @title Environmental profiling with one-classification support vector machine
#' @description Presence-only modelling and classification of coordinates predicted as 
#' presence and absence
#' 
#' @param xy Data frame or list of data frames with coordinates (each row is a point)
#' @param varstack RasterStack of variables for modelling
#' @param bbs.grid Object derived from function \code{\link[mopa]{delimit}}. 
#' Matrix or list of matrixes of the background xy coordinates in columns. 
#' @param nu parameter needed for one-classification \code{\link[e1071]{svm}}. 
#' Default is 0.5
#'
#' @seealso \code{\link[e1071]{svm}} 
#' 
#' @return  A list with two components:
#'   \item{absence }{Matrix or list of matrixes with xy coordinates predicted 
#' as absence, each matrix correspond to a different background extent}
#' \item{presence }{Matrix or list of matrixes with xy coordinates predicted as presence, 
#' each matrix correspond to a different background extent}
#' 
#' @details This function constitutes the 
#' first step from a three-step proccess to generate pseudo-absences, and is aimed at 
#' separating from background (see \code{\link[mopa]{delimit}}) those areas (xy records) 
#' suitable for the species. 
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
#' ##plot
#' plot(unsuitable.bg$absence$H11, pch="*")
#' points(unsuitable.bg$presence$H11, pch="*", col= "pink4")
#' }
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @export
#' 
#' @importFrom e1071 svm
 

OCSVMprofiling<-function(xy, varstack, bbs.grid, nu=.5){
  if (class(xy) != "list")  xy <- list(xy)
  if (class(bbs.grid) != "list")  bbs.grid <- rep(list(bbs.grid), length(xy))
  bioclim <-varstack
  absence <- list()    
  presence <- list()
  for(i in 1:length(xy)){
    length(absence) <- i
    coo <- bbs.grid[[i]]
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
