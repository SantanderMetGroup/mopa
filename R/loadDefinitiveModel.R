


#' @title Load Rdata files storing definitive fitted model 
#' @description Load information from Rdata generated with function allModelling based 
#' in an index to select the definitive fitted model/s, ideally, the index returned by 
#' function indextent should be used.
#' 
#' @param models Character of the rdata files to be loaded. Including, at least,
#'  the attribute "species" (a character of names corresponding to the modelled species/populations).
#' @param slot Any character of the following: "allmod", auc", "kappa", "tss", "mod", "p"
#' 
#'  
#' 
#' @return Depending on the specified slot:
#'  \item{allmod }{fitted model using all data for training}
#'  \item{auc }{AUC statistic in the cross validation}
#'  \item{kappa }{kappa statistic in the cross validation}
#'  \item{tss }{true skill statistic in the cross validation }
#'  \item{mod }{fitted model with partitioned data}
#'  \item{p }{cross model prediction}
#'  
#' 
#' @details detail.
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
#' auc_mars <-loadTestValues(models = modirs, test = "auc") 
#' ind <- indextent(testmat = auc_mars, diagrams = F)
#' def <- loadDefinitiveModel(models = ind, slot = "allmod")
#' }
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @export


loadDefinitiveModel<-function(models, slot = c("allmod", "auc", "kappa", "tss", "mod", "p")){
  slot <- match.arg(slot, choices = c("allmod", "auc", "kappa", "tss", "mod", "p"))
  species <- attr(models, "species") 
  modelslot<-list()
  for (i in 1:length(species)){
    g <- species[i]
    mod <- get(load(models[i]))
    if (slot == "allmod"){modelslot[[i]]<-mod$allmod
    } else if (slot == "auc"){modelslot[[i]]<-mod$auc
    } else if (slot == "kappa"){modelslot[[i]]<-mod$kappa
    } else if (slot == "tss"){modelslot[[i]]<-mod$tss
    } else if (slot == "mod"){modelslot[[i]]<-mod$mod
    } else if (slot == "p"){modelslot[[i]]<-mod$p}
  } 
  names(modelslot) <- species
  suppressWarnings(
    if(length(species) == 1) modelslot <- modelslot[[1]]
  )
  return (modelslot)
}



