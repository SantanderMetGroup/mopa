
#' @title Load model validation statistics from Rdata generated with function allModelling
#' @description Load model validation statistics from Rdata generated with function allModelling
#' 
#' @param data Object with the same structure as the object returned by function bindPresAbs. 
#' @param test Any character of the following: "auc", "kappa", "tss".
#' @param algorithm Any character of the following: "glm", "svm", "maxent", "mars", "randomForest", 
#' "cart.rpart" or "cart.tree"
#' @param sourcedir Character of the path where Rdata objects are stored
#' 
#'  
#' 
#' @return Matrix with values for each species and background extents.
#' 
#' #@details detail.
#' 
#' 
#' 
#' @author M. Iturbide \email{maibide@@gmail.com}
#' 
#' @examples
#' \dontrun{
#' data(presaus)
#' data(biostack)
#' ##modeling
#' modirs <-allModeling(data = presaus, varstack = biostack, k = 10, "mars") 
#' ##loading
#' auc_mars <-loadTestValues(data = presaus, test = "auc", algorithm = "mars")
#'  
#' library(lattice)
#' levelplot(auc_mars, aspect = 5, 
#'           scales = list(y = list(cex = 0.8, 
#'                                  at = c(1, 49, ncol(auc_mars)), 
#'                                  labels = c(colnames(auc_mars)[1],
#'                                             colnames(auc_mars)[49],
#'                                             colnames(auc_mars)[ncol(auc_mars)]))), 
#'           at = seq(0.6, 1, 0.01), col.regions = bpy.colors,
#'           xlab = "Haplogroups", ylab = "Background extent (km)", main = "AUC")
#' }
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @export



loadTestValues <-function(models, sourcedir = getwd(), test = c("auc", "kappa", "tss")){
  test <- match.arg(test, choices = c("auc", "kappa", "tss"))
  extents <-attr(models, "extents")
  extent.lengths <- lengths(extents)
  max.extents <- extents[[which(extent.lengths == max(extent.lengths))]]
  algorithm <- attr(models, "algorithm")  
  species <- attr(models, "species") 
  
  val_mat <- matrix(NA, nrow= length(species), ncol=length(max.extents), dimnames=list(species, max.extents))
  
  for (i in 1:length(species)){
    g <- species[i]
    r <- extents[[i]]
    print(paste("loading values for species", g))
    for (k in 1:length(r)){
      if (length(species) < 2){
        load(paste(sourcedir,"/", algorithm, "_bg", r[k],".rda",sep=""))
      }else {
        load(paste(sourcedir,"/", algorithm, "_bg", r[k],"_hg",g,".rda",sep=""))
      }
      if (test == "auc"){
        val_mat[i,k]<-mod$auc
      }else if (test == "kappa"){
        val_mat[i,k]<-mod$kappa
      }else if (test == "tss"){
        val_mat[i,k]<-mod$tss
      }
    } 
  }
  attr(val_mat, "algorithm") <- algorithm
  attr(val_mat, "species") <- species
  attr(val_mat, "extents") <- extents
  attr(val_mat, "source directory") <- sourcedir
  return (val_mat)
}
