
#' @title Load model validation statistics from Rdata generated with function allModelling
#' @description Load model validation statistics from Rdata generated with function allModelling
#' 
#' @param models Object with the same structure as the object returned by function allModelling.
#' @param sourcedir Directory of the objects generated and saved by allModelling
#' @param test Any character of the following: "auc", "kappa", "tss".
#' "cart.rpart" or "cart.tree"
#' 
#'  
#' 
#' @return Matrix with values for each species and background extents.
#' 
#' #@details detail.
#' 
#' 
#' 
#' @author M. Iturbide 
#' 
#' @examples
#' \dontrun{
#' data(presaus)
#' data(biostackENSEMBLES)
#' ##modeling
#' modirs <-allModeling(data = presaus, varstack = biostackENSEMBLES$baseline, k = 10, "mars") 
#' ##loading
#' auc_mars <- loadTestValues(models = modirs, test = "auc") 
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



loadTestValues <-function(models, test = c("auc", "kappa", "tss")){
  test <- match.arg(test, choices = c("auc", "kappa", "tss"))
  extents <- models$extents
  extent.lengths <- lengths(extents)
  max.extents <- extents[[which(extent.lengths == max(extent.lengths))[1]]]
  algorithm <- models$algorithm
  species <- models$species
  val_mat <- matrix(NA, nrow= length(species), ncol = length(max.extents), dimnames=list(species, max.extents))
  for (i in 1:length(species)){
    message(paste("loading values for species", species[i]))
    for (k in 1:length(extents[[i]])){
      mod <- get(load(models$dirs[[i]][k]))
      if (test == "auc"){
        val_mat[i,k] <- mod$auc
      }else if (test == "kappa"){
        val_mat[i,k]<-mod$kappa
      }else if (test == "tss"){
        val_mat[i,k]<-mod$tss
      }
    } 
  }
  message("Done.")
  attr(val_mat, "dirs") <- models$dirs
  attr(val_mat, "algorithm") <- algorithm
  attr(val_mat, "species") <- species
  attr(val_mat, "extents") <- extents
  return (val_mat)
}
