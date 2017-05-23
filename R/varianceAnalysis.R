#' @title Variance analysis of RasterStack objects
#' @description Extract componets of lists of objects 
#' (as returned by function \code{\link[mopa]{mopaPredict}}) 
#' and perform variance analysis to obtain raster objects of the 
#' contribution of each component to the observed varaibility.
#' 
#' 
#' @param predictions listed lists of raster objects as returned by \code{\link[mopa]{mopaPredict}}
#' @param component1 Character. Options are "SP", "PA", "SDM", "baseClim" and "newClim" (see Details). If exist, "foldModel" is 
#' another option. Selected option corresponds to the first component in 
#' the variance analysis.
#' @param component2 Character. Options are "SP", "PA", "SDM", "baseClim" and "newClim" (see Details). If exist, "foldModel" is 
#' another option. Selected option corresponds to the second component in 
#' the variance analysis.
#' @param fixed Optional. Character of the component names corresponding to the components that are not being 
#' analyzed (component3, component4...). One name for each component is provided, components that only have one choice 
#' (e.g. a single species, a single baseline climate, etc.) are internally fixed.
#'  If \code{fixed = NULL} the first element of each component is selected.
#' If \code{fixed} is specified, the selected name must be provided for each of the components that have multi-choices.
#' 
#' @details Rasters are extracted using function \code{\link[base]{grep}}, by matching
#' names in the lists and characters in \code{componen1} and \code{componen2}. 
#' The contribution of componen1 in front component2 to the spread (uncertainty) of the projected probabilities 
#' in \code{predictions} is here assessed using a simple analysis of variance approach, where the total 
#' variance (V) can be decomposed as the summation of the variance explained by component1 (Vcomp1), 
#' component2 (Vcomp2) and the combination of the previous two (Vcomp12):
#' 
#' \eqn{V = Vcomp1 + Vcomp2 + Vcomp2}.
#' 
#' Description of the components:
#' \itemize{
#' \item{SP:} presence data sets 
#' \item{PA:} pseudo-absence realizations 
#' \item{SDM:} modeling algorithms 
#' \item{baseClim:} bseline climate, i.e. sets of vaiables used for model calibration in function
#'  \code{\link[mopa]{mopaTrain}}, 
#' \item{newClim:} new climate, i.e. sets of vaiables used to project models (e.g. future climate projections) in
#' function \code{\link[mopa]{mopaPredict}} .
#' }
#' 
#' 
#' @return A list of two RasterStack objects, the first containing the global mean and standard deviation and the 
#' second containing the percentage of variance correponding to each component in the analysis 
#' (component1, component2 and components 1 and 2). 
#' 
#' 
#' @author M. Iturbide 
#' 
#' @references San Martin, D., Manzanas, R., Brands, S., Herrera, S., & Gutierrez, J.M. (2016) Reassessing 
#' Model Uncertainty for Regional Projections of Precipitation with an Ensemble of Statistical Downscaling Methods. 
#' Journal of Climate 30, 203-223.
#' 
#' 
#' @examples
#' data(Oak_phylo2)
#' presences <- Oak_phylo2$H11
#' data(biostack)
#' projection(biostack$baseline) <- CRS("+proj=longlat +init=epsg:4326")
#' r <- biostack$baseline[[1]]
#' 
#' ## Background of the whole study area
#' bg <- backgroundGrid(r)
#' 
#' ## Considering an unique background extent
#' RS_random <-pseudoAbsences(xy = presences, background = bg$xy, 
#'                            realizations = 10,
#'                            exclusion.buffer = 0.083*5, prevalence = -0.5, kmeans = FALSE)
#' 
#' fittingRS <- mopaTrain(y = RS_random, x = biostack$baseline, k = 10, 
#'                        algorithm = "glm", weighting = TRUE)
#' 
#' modsRS <- extractFromModel(models = fittingRS, value = "model")
#' 
#' # MODEL PREDICTION AND ANALYSIS OF THE VARIABILITY IN PROJECTIONS
#' prdRS.fut <- mopaPredict(models = modsRS, newClim = biostack$future)
#' result <- varianceAnalysis(prdRS.fut, "PA", "newClim")
#' spplot(result$variance, col.regions = rev(get_col_regions()))
#' 
#' 
#' 
#' @export
#' @importFrom stats sd

varianceAnalysis <- function(predictions, component1, component2, fixed = NULL) {
  namescomps <-  c(component1, component2, paste(component1, "and", component2))
  d <- depth(predictions) -1
  choices <- c("SP", "PA", "SDM", "baseClim", "newClim", "foldModel")[1:d]
  component1 <- match.arg(component1, choices = choices)
  component2 <- match.arg(component2, choices = choices)
  if(component1 == component2) stop("component1 and 2 are equal, please select another component")
  wl1 <- which(choices == component1)
  wl2 <- which(choices == component2)
  dl <- depthLength(predictions)
  dl1 <- dl[wl1]
  dl2 <- dl[wl2]
  if(dl1 == 1) stop("there is only one choice for ", component1, " please chose another component")
  if(dl2 == 1) stop("there is only one choice for ", component2, " please chose another component")
  stick <- fixed
  if(is.null(stick)){
    stick <- depthnames1(predictions)[-c(wl1, wl2)]
  }else{
    stick0 <- depthnames1(predictions)[which(dl ==1)]
    stick <- unique(c(stick, stick0))
  }
  component1 <- depthnames(predictions, wl1)
  component2 <- depthnames(predictions, wl2)
  
  comp2 <- list()
  for(i in 1:length(component2)){
    comp2[[i]] <- extractFromPrediction(predictions, component2[i])
  }
  names(comp2) <- component2
  
  comp1 <- list()
  for(i in 1:length(component1)){
    comp1[[i]] <- extractFromPrediction(comp2, component1[i])
  }
  
    a <- length(stick) 
    i <- 0
    if(d - a != 2){
      stop("There are components in prediction with multiple choices, please 
           set argument 'fixed' correctly to select the component member/s that is/are kept 
           constant in the analysis")
    }
    while(a!= 0){
      a <- a-1
      i <- i+1
      comp10 <- stack(unlist(comp1))
      comp1 <- extractFromPrediction(comp10, stick[i])
    }
    
  
  
  bothcomp <- stack(unlist(comp1))
  names(bothcomp)
   if(length(component1)*length(component2) != nlayers(bothcomp)) stop("The specified component names do not match the names in the predictions", call. = FALSE)
  datos <- array(NA, dim = c(length(component1)*length(component2), ncell(bothcomp)), dimnames = list(names(bothcomp)))
  for(i in 1:nlayers(bothcomp)){
    datos[i, ] <- bothcomp[[i]]@data@values
  }
  
  mediaGlobal <- apply(datos, FUN = "mean", MARGIN=2, na.rm = TRUE)
  varGlobal <- apply(datos, FUN = "sd", MARGIN=2, na.rm = TRUE)*sqrt((nrow(datos)-1)/nrow(datos))
  
  
  mediacomp1 <- matrix(data = NA, nrow = length(component1), ncol = ncol(datos))
  for(i in 1:length(component1)){
    indcomp1 <- ((i-1)*length(component2)+1):(i*length(component2))
    mediacomp1[i,] <- apply(datos[indcomp1,], FUN = "mean", MARGIN=2, na.rm = TRUE)
  }
  rownames(mediacomp1) <-component1
  
  mediacomp2 <- matrix(data = NA, nrow = length(component2), ncol = ncol(datos))
  for(i in 1:length(component2)){
    mediacomp2[i,] <- apply(datos[seq(i,nrow(datos),length(component2)),], FUN = mean, MARGIN=2, na.rm = TRUE)
  }
  dos <- apply(mediacomp2, FUN = "sd", MARGIN=2, na.rm = TRUE)*sqrt((length(component2)-1)/length(component2))
  uno <- apply(mediacomp1, FUN = "sd", MARGIN=2, na.rm = TRUE)*sqrt((length(component1)-1)/length(component1))
  
  uno.dos <- matrix(data = NA, nrow = nrow(datos), ncol = ncol(datos))
  for(i in 1:length(component1)){
    for(j in 1:length(component2)){
      uno.dos[(i-1)*length(component2)+j,] <- (datos[(i-1)*length(component2)+j,] - mediacomp2[j,] - mediacomp1[i,] + mediaGlobal)^2
    }
  }
  
  uno.dos <- apply(uno.dos, FUN = "mean", MARGIN=2, na.rm = TRUE)
  
  # plot(dos^2+uno^2+uno.dos, typ = "l")
  # lines(varGlobal^2, col = "red")
  # sd(dos^2+uno^2+uno.dos - varGlobal, na.rm = T)
  # mean(dos^2+uno^2+uno.dos - varGlobal, na.rm = T)
  
  uno100 <- uno^2 *100 / (uno^2+dos^2+uno.dos)
  dos100 <- dos^2 *100 / (uno^2+dos^2+uno.dos)
  uno.dos100 <- uno.dos *100 / (uno^2+dos^2+uno.dos)
  nan.ind <- which(uno100 == "NaN")
  dos100[nan.ind] <- 0
  uno100[nan.ind] <- 0
  uno.dos100[nan.ind] <- 0
  
  
  bothcomp[[1]]@data@values <- uno100
  bothcomp[[2]]@data@values <- dos100
  bothcomp[[3]]@data@values <- uno.dos100
  bothcomp[[4]]@data@values <- varGlobal
  bothcomp[[5]]@data@values <- mediaGlobal
  l1 <- stack(bothcomp[[5]], bothcomp[[4]])
  l2 <- stack(bothcomp[[1]], bothcomp[[2]], bothcomp[[3]])
  
  names(l1) <- c("mean", "sd")
  names(l2) <- namescomps
  
  return(list("mean" = l1, "variance" = l2))
}


#end


#' @title Depth length in a list
#' @description Depth length in a list
#' @param this List
#' @keywords internal
#' 
depthLength <- function(this){
  that <- this
  i <- 0
  len <- numeric()
  while(is.list(that)){
    i <- i + 1
    len[i] <- length(that)
    that <- that[[1]]
  }
  return(len)
}
#end

#' @title Depth names in a list 1
#' @description Depth length in a list 1
#' @param this List
#' @keywords internal
depthnames1 <- function(this){
  that <- this
  i <- 0
  len <- character()
  while(is.list(that)){
    i <- i + 1
    len[i] <- names(that)[1]
    that <- that[[1]]
  }
  return(len)
}
#end

#' @title Depth length in a list
#' @description Depth length in a list
#' @param this List
#' @keywords internal
depthnames <- function(this, ind){
  that <- this
  i <- 0
  while(is.list(that)){
    i <- i + 1
    if(i == ind){
      len <- names(that)
    }
    that <- that[[1]]
  }
  return(len)
}
#end