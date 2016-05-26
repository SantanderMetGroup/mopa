
#' @title Easy species distribution modelling and cross validation 
#' into backgrounds of different extent
#' @description Species distribution modelling and k-fold cross validation 
#' for a set of presence/absence data per species corresponding to a different background 
#' extent. Algorithms supported are "glm", "svm", "maxent", "mars", "randomForest", "cart.rpart" 
#' and "cart.tree" 
#' 
#' @param data Object returned by function bindPresAbs or list/s of data frames with coordinates in the first two columns and presence/absence 
#' (1=presence, 0=absence) in the third column. 
#' @param varstack RasterStack of variables for modelling
#' 
#' @param k Integer. Number of folds for cross validation. Default is 10
#' @param algorithm Any character of the following: "glm", "svm", "maxent", "mars", "randomForest", "cart.rpart" 
#' or "cart.tree"
#' @param destdir Character of the output path
#' @param projection Object of class CRS with the coordinate reference system. Default is 
#' CRS("+proj=longlat +init=epsg:4326") 
#' 
#'  
#' 
#' @return Named Rdata objects are stored in the 
#' specified path. Each Object is given a name indicating the algorithm, background 
#' extent, and species in this order (if a single species is provided no name is given 
#' for the species). Character object with listed files is returned.
#' Each Rdata consists of a list with six components:
#' 
#'  \item{allmod }{fitted model using all data for training}
#'  \item{auc }{AUC statistic in the cross validation}
#'  \item{kappa }{kappa statistic in the cross validation}
#'  \item{tss }{true skill statistic in the cross validation }
#'  \item{mod }{fitted model with partitioned data}
#'  \item{p }{cross model prediction}
#'  
#' 
#' @details This function calculates the AUC with the function \code{\link[PresenceAbsence]{auc}} from package 
#' \pkg{PresenceAbsence}. \strong{Note:} Package \pkg{SDMTools} must be detached.
#' 
#' 
#' 
#' @author M. Iturbide \email{maibide@@gmail.com}
#' 
#' @examples
#' 
#' \dontrun{
#' data(presaus)
#' data(biostack)
#' ##modelling
#' modirs <- allModeling(data = presaus, varstack = biostack, k = 10, algorithm = "mars")
#' }
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @export
#' 
#' @import raster
#' @import sp
#' @import PresenceAbsence
#' @importFrom sampling strata  
#' @importFrom e1071 best.svm 
#' @importFrom dismo maxent
#' @importFrom earth earth
#' @importFrom rpart rpart
#' @importFrom tree tree
#' @importFrom randomForest randomForest
#' 


allModeling <- function(data, varstack, k = 10, algorithm = c("glm", "svm", "maxent", "mars", "randomForest", "cart.rpart", "cart.tree"), 
                        destdir =getwd(), projection = CRS("+proj=longlat +init=epsg:4326")){
  algorithm <- match.arg(algorithm, choices = c("glm", "svm", "maxent", "mars", "randomForest", "cart.rpart", "cart.tree"))
  biostack <- varstack
  if (class(data[[1]]) != "list"){
    data<-list(data)
    names(data) <- "species"
  }else{
    data <- data
  }
  
  for (i in 1:length(data)){
    sp_01 <- data[[i]]
    
    for(j in 1:length(sp_01)){
      #print(paste("running model for species", i, "considering pseudo-absences inside the extent of", names (sp_01)[j]))
      destfile <- names(sp_01)[j]
      
      sp.bio <- biomat(sp_01[[j]], biostack)
      x <- kfold(k, df = sp.bio)
      xx <- leaveOneOut(x)
      mod <- tryCatch({modelo(kdata = xx, data=sp.bio, algorithm)},
                      error = function(err){xxx = list(rep(NA, k), NA, NA)})
      if (length(data)==1){
       save(list=c("mod"), file=paste(destdir, "/", algorithm,"_bg", destfile, ".rda",sep=""))
      }else{
        save(list=c("mod"), file=paste(destdir, "/", algorithm,"_bg", destfile, "_hg",names(data)[i], ".rda",sep=""))
      }
      rm(mod, xx, x, sp.bio)
    }}
  #collect information
  extents <- list()
  for(j in 1:length(data)){
    extents[[j]] <- names(data[[j]])
  }
  names(extents) <- names(data)
  dirs <- list.files(destdir, full.names = F, pattern = algorithm)
  attr(dirs, "algorithm") <- algorithm
  attr(dirs, "species") <- names(data)
  attr(dirs, "extents") <- extents
  
  return(dirs)
}
