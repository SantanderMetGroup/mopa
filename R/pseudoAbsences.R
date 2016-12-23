
#' @title Pseudo-absences 
#' @description Pseudo-absence data generation at random or by k-means clustering inside backgrounds of 
#' different extent
#' 
#' @param xy Data frame or list of data frames with coordinates (each row is a point)
#' @param background Matrix or list/s of matrixes with background coordinates in columns.
#' Object derived from function delimit, OCSVMprofiling or bgRadio. 
#' @param exclusion.buffer value of the minimum distance to be kept between presence data and 
#' pseudo-absence data. Default is 0.0166
#' @param prevalence Proportion of presences against absences. Default is 0.5 (equal number of 
#' presences and absences)
#' @param kmeans Logical. If FALSE (default) pseudo-absences are generated at random. If TRUE
#' k-means clustering of the background is done and centroids are extracted as pseudo-absences.
#' @param varstack RasterStack of variables for modelling
#' @param projection Object of class CRS with the coordinate reference system. Default is 
#' CRS("+proj=longlat +init=epsg:4326") 
#' 
#' @return  List/s of data frames 
#' 
#' 
#' @details Details.The application of this function is the third step in a three-step proccess 
#' to generate pseudo-absences. 
#' 
#' 
#' 
#' 
#' @author M. Iturbide 
#' 
#' @examples
#' \dontrun{
#' ##delimit study area
#' data(Oak_phylo2)
#' data(sp_grid)
#' oak.extension<-boundingCoords(Oak_phylo2)
#' box.grid<-delimit(oak.extension, sp_grid, names(Oak_phylo2))
#' ## environmental profiling
#' data(biostack)
#' unsuitable.bg <-OCSVMprofiling(xy = Oak_phylo2, varstack = biostack, 
#' bbs.grid = box.grid$bbs.grid)
#' ## sequence of 10 km between distances, from 20 km to the length of the 
#' ##half diagonal of the bounding box.
#' ext <-bgRadio(xy = Oak_phylo2, bounding.coords = oak.extension, 
#' bg.absence = unsuitable.bg$absence, start = 0.166, by = 0.083, unit = "decimal degrees")
#' ## pseudo-absence generation at random
#' pa_random <-PseudoAbsences(xy = Oak_phylo2, background = ext, 
#' exclusion.buffer = 0.083, prevalence = 0.5, kmeans = FALSE)
#' ##plot
#' plot(ext$H11[[5]], pch="*", col= "grey", cex=.5)
#' points(pa_random$H11[[5]], col="red", pch=".", cex=4)
#' points(Oak_phylo2$H11, col="blue", pch=".", cex=3)
#' }
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @export
#' 
#' @import sp
#' @importFrom spatstat disc
#' @importFrom stats kmeans
#' 
#' 



pseudoAbsences <- function (xy, background, exclusion.buffer = 0.0166, prevalence = 0.5, 
                          kmeans = FALSE, varstack = NULL, projection = CRS("+proj=longlat +init=epsg:4326")){
  polybuffs <- list()
  r <- exclusion.buffer
  prev <- (1 - prevalence) * 2
  if (any(c("data.frame", "matrix") == class(xy))) xy <- list(xy)
  if (any(c("data.frame", "matrix") == class(background))){
    background <- rep(list(background), length(xy))
    if(length(xy) > 1) message("The same background will be used for all presence datasets in xy")
  } 
  if(length(xy) != length(background)) stop("xy and background do not have the same length")
  if(any(c("matrix", "data.frame") == class(background[[1]]))){
    background <- lapply(seq(length(background)), function(x){list(background[[x]])})
  } 
  for (j in 1:length(xy)) {
    pres.1km <- xy[[j]]
    polys <- list()
    for (i in 1:nrow(pres.1km)) {
      discbuff <- disc(radius = r, centre = c(pres.1km[i, 
                                                       1], pres.1km[i, 2]))
      discpoly <- Polygon(rbind(cbind(discbuff$bdry[[1]]$x, 
                                      y = discbuff$bdry[[1]]$y), c(discbuff$bdry[[1]]$x[1], 
                                                                   y = discbuff$bdry[[1]]$y[1])))
      polys <- c(polys, discpoly)
    }
    spolys <- list()
    for (i in 1:length(polys)) {
      spolybuff <- Polygons(list(polys[[i]]), ID = i)
      spolys <- c(spolys, spolybuff)
      spol <- SpatialPolygons(spolys)
       proj4string(spol) <- projection
    }
    polybuffs[[j]] <- spol
  }
  aa <- list()
  for (j in 1:length(background)) {
    message("generating pseudo-absences for species ", 
                j, " out of ", length(background))
    aus <- list()
    coords.l <- background[[j]]
    polpol <- polybuffs[[j]]
    pr <- xy[[j]]
    for (i in 1:length(coords.l)) {
#       print(paste("b =", i, "out of", as.character(length(coords.l))))
      coords <- coords.l[[i]]
      sp.coords <- SpatialPoints(coords)
       proj4string(sp.coords) <- projection
      a <- over(sp.coords, polpol)
      abs.bg <- coords[which(is.na(a)), 1:2]
      
      if (kmeans == TRUE) {
        abs.aus<-cbind(abs.bg, rep(0, nrow(abs.bg)))
        abs.bio<-biomat(data = abs.aus, varstack, projection)
        aus[[i]] <- tryCatch({kmeans(cbind(abs.bg, abs.bio[,-1]), centers = prev * nrow(pr))$centers[,1:2]},
                              error = function(err){aus[[i]] <- NULL})
        
      }
      else {
       
        aus[[i]] <- tryCatch({abs.bg[sample(1:nrow(abs.bg), size = prev * 
                                    nrow(pr)), ]}, error = function(err){aus[[i]] <- NULL})
 
      }
    }
    
    names(aus) <- names(coords.l)
    ind <- rep(NA,length(aus))
    for(n in 1:length(aus)){
      if(is.null(aus[[n]])){
        message("Background ", names(coords.l)[n], " is too small for sampling and will be ignored")
        ind[n] <- n
      }
    }
    as <- unname(na.omit(ind))
    if(length(as)!=0){
      aa[[j]] <- aus[-as]
    }else{
      aa[[j]] <- aus
    }
  }
  if (length(aa) == 1) {
    aa <- aa[[1]]
  } else {
    names(aa) <- names(xy)
  }
  return(aa)
}

