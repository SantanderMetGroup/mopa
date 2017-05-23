
#' @title Pseudo-absence data generation
#' @description Pseudo-absence data generation at random or by k-means clustering inside a single
#' background or a group of backgrounds (e.g. of different extent, \code{\link[mopa]{backgroundRadius}})
#' @param xy Data frame or list of data frames with coordinates (each row is a point), this is, 
#' presence data
#' @param background Matrix or list/s of matrixes with background coordinates in columns.
#' Object derived from function \code{\link[mopa]{backgroundGrid}}, \code{\link[mopa]{OCSVMprofiling}} 
#' or \code{\link[mopa]{backgroundRadius}}. 
#' @param realizations Integer. Number of realizations (default = 1).
#' @param exclusion.buffer value of the minimum distance to be kept between presence data and 
#' pseudo-absence data. Default is 0.0166
#' @param prevalence Proportion of presences against absences. Default is 0.5 (equal number of 
#' presences and absences)
#' @param kmeans Logical. If FALSE (default) pseudo-absences are generated at random. If TRUE
#' k-means clustering of the background is done and centroids are extracted as pseudo-absences.
#' @param varstack RasterStack of variables for to compute the k-means clustering. Used if \code{kmeans}
#' = TRUE.
#' 
#' 
#' @return  data frame or list/s of data frames 
#' 
#' 
#' @details Details. The application of this function could be preceded by the application
#' of functions \code{\link[mopa]{OCSVMprofiling}} and/or \code{\link[mopa]{backgroundRadius}}
#' in order to consider alternative methods for pseudo-absence data generation (see references).
#' 
#' @seealso \code{\link[mopa]{mopaTrain}}
#' 
#' 
#' @author M. Iturbide 
#' 
#' @examples
#' data(Oak_phylo2)
#' data(biostack)
#' projection(biostack$baseline) <- CRS("+proj=longlat +init=epsg:4326")
#' r <- biostack$baseline[[1]]
#' ## Background of the whole study area
#' bg <- backgroundGrid(r)
#' 
#' ## environmental profiling
#' bg.profiled <- OCSVMprofiling(xy = Oak_phylo2, varstack = biostack$baseline, 
#'                               background = bg$xy)
#' 
#' ## considering an unique background extent
#' RS_random <-pseudoAbsences(xy = Oak_phylo2, background = bg$xy, 
#'                            exclusion.buffer = 0.083*5, 
#'                            prevalence = -0.5, kmeans = FALSE)
#' RSEP_random <-pseudoAbsences(xy = Oak_phylo2, background = bg.profiled$absence, 
#'                              exclusion.buffer = 0.083*5, 
#'                              prevalence = -0.5, kmeans = FALSE)
#' 
#' ## different background extents
#' bg.extents <- backgroundRadius(xy = Oak_phylo2, background = bg$xy, 
#'                                start = 0.166, by = 0.083*20, 
#'                                unit = "decimal degrees")
#' 
#' ## inside different background extents
#' TS_random <-pseudoAbsences(xy = Oak_phylo2, background = bg.extents, 
#'                            exclusion.buffer = 0.083*10, 
#'                            prevalence = -0.5, kmeans = FALSE)
#' 
#' 
#' ## with k-means clustering
#' TS_kmeans <-pseudoAbsences(xy = Oak_phylo2, background = bg.extents, 
#'                            exclusion.buffer = 0.083*5, 
#'                            prevalence = -0.5, kmeans = TRUE, 
#'                            varstack = biostack$baseline)
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @export
#' 
#' @import sp
#' @importFrom spatstat disc
#' @importFrom stats kmeans na.omit
#' @importFrom gtools combinations
#' 
#' 



pseudoAbsences <- function (xy, background, realizations = 1, exclusion.buffer = 0.0166, prevalence = 0.5, 
                            kmeans = FALSE, varstack = NULL){
  
  if (any(c("data.frame", "matrix") == class(xy))) xy <- list(xy)
  if (any(c("data.frame", "matrix") == class(background))){
    background <- rep(list(background), length(xy))
    if(length(xy) > 1) message("The same background will be used for all presence datasets in xy")
  } 
  if(length(xy) != length(background)) stop("xy and background do not have the same length")
  if(any(c("matrix", "data.frame") == class(background[[1]]))){
    background <- lapply(seq(length(background)), function(x){list(background[[x]])})
  }
  spa <- list()
  for(j in 1:length(xy)){
    pa <- list()
    nm <- character()
    message("[", Sys.time(), "] Generating pseudo-absences for species ", j)
    xy1 <- xy[[j]]
    background1 <- background[[j]]
    for(i in 1:realizations){
      message(":::[", Sys.time(), "] Realization ", i)
      pa[[i]] <- pseudoAbsences0(xy1, background1, exclusion.buffer = exclusion.buffer, prevalence = prevalence,
                                 kmeans =  kmeans, varstack = varstack)
      if(i < 10){
        nm[i] <- paste0("0", i)
      }else{
        nm[i] <- as.character(i)
      }
    }
    names(pa) <- paste0("PA", nm)
    spa[[j]] <- pa
  }
  if(is.null(names(xy))) names(xy) <- paste0("species", 1:length(xy))
  names(spa) <- names(xy)
  return(spa) 
}

#end



#' @title Pseudo-absences internal
#' @description Pseudo-absence data generation at random or by k-means clustering inside a single
#' background or a group of backgrounds (e.g. of different extent, \code{\link[mopa]{backgroundRadius}})
#' @param xy Data frame or list of data frames with coordinates (each row is a point), this is, 
#' presence data
#' @param background Matrix or list/s of matrixes with background coordinates in columns.
#' Object derived from function \code{\link[mopa]{backgroundGrid}}, \code{\link[mopa]{OCSVMprofiling}} 
#' or \code{\link[mopa]{backgroundRadius}}. 
#' @param exclusion.buffer value of the minimum distance to be kept between presence data and 
#' pseudo-absence data. Default is 0.0166
#' @param prevalence Proportion of presences against absences. Default is 0.5 (equal number of 
#' presences and absences)
#' @param kmeans Logical. If FALSE (default) pseudo-absences are generated at random. If TRUE
#' k-means clustering of the background is done and centroids are extracted as pseudo-absences.
#' @param varstack RasterStack of variables for to compute the k-means clustering. Used if \code{kmeans}
#' = TRUE.
#' 
#' 
#' @author M. Iturbide 
#' 
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @keywords internal
#' @import sp
#' @importFrom spatstat disc
#' @importFrom stats kmeans na.omit
#' 



pseudoAbsences0 <- function(xy, background, exclusion.buffer = 0.0166, prevalence = 0.5, 
                            kmeans = FALSE, varstack = NULL){
  polybuffs <- list()
  r <- exclusion.buffer
  prev <- (1 - prevalence) * 2
  pr <- xy
  polys <- list()
  for (i in 1:nrow(pr)) {
    discbuff <- disc(radius = r, centre = c(pr[i, 
                                                     1], pr[i, 2]))
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
    # proj4string(spol) <- projection
  }
  
  
    aus <- list()
    coords.l <- background
    for (i in 1:length(coords.l)) {
      #       print(paste("b =", i, "out of", as.character(length(coords.l))))
      coords <- coords.l[[i]]
      sp.coords <- SpatialPoints(coords)
      # proj4string(sp.coords) <- projection
      a <- over(sp.coords, spol)
      abs.bg <- coords[which(is.na(a)), 1:2]
      
      if (kmeans == TRUE) {
        abs.aus <- cbind(abs.bg, rep(0, nrow(abs.bg)))
        if(length(abs.aus) != 0){
          abs.bio<-biomat(data = abs.aus, varstack)#, projection)
          aus[[i]] <- kmeans(cbind(abs.bg, abs.bio[,-1]), centers = prev * nrow(pr))$centers[,1:2]
        }else{
          aus[[i]] <- NULL
        }
      }else {
        
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
      aa <- aus[-as]
    }else{
      aa <- aus
    }
  
  ab <- bindPresAbs(xy, aa)
  return(ab)
}


#end



#' @title Bind presences and absences  
#' @description Binds presence and absence data for each background extension 
#' 
#' @param presences Data frame or list of data frames with coordinates for presence data 
#' (each row is a point)
#' @param absences Object returned by function \code{\link[mopa]{pseudoAbsences}}. 
#' List/s of data frames with coordinates for absence data (each row is a point)
#' 
#' @return  List/s of matrixes with xy coordinates for presence/pseudo-absence data.
#' Each matrix correspond to a different background extent
#' 
#' 
#'
#' @keywords internal
#' 
#' 
#' @author M. Iturbide




bindPresAbs <- function (presences, absences){
  pres <- presences
  prau<-list()
  pr <- cbind(pres, rep(1, nrow(pres)))
  names(pr)<-c("x", "y", "v")
  au <- absences
    for (j in 1:length(au)){
      aj <- cbind(as.data.frame(au[[j]]), rep(0,nrow(au[[j]])))
      names(aj)<-names(pr)
      prau[[j]]<-rbind(pr, aj)
    }
  names(prau)<-names(au)
  presaus <-prau
  rm(aj, pr, au, prau)
  return(presaus)
}

#end
