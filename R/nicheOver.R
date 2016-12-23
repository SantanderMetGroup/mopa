#' @title Niche overlap
#' @description Compute niche overlap among rasters in a RasterStack
#' 
#' @param stack RasterStack   
#' @param metric Metric for niche overlap. Options are "D" and "O" (see details).
#' 
#' @return Matrix of overlap values for metric D or O.
#' 
#' 
#' @details Niche overlap measures the similarity of the environmental ranges occupied by each constructed model
#' via operating the difference between two vectors of probability distributions (p), where p x,i and p y,i are 
#' the normalized suitability scores for biological entity or model X and Y in grid cell i.
#' 
#'  D: Schoener’s statistic for niche overlap (Warren et al., 2008).
#'  O: Pianka index (Pianka, 1973).
#' 
#' 
#' 
#' @author M. Iturbide 
#' 
#' @references Warren DL, Glor RE, Turelli M, Funk D (2008) Environmental niche equivalency versus conservatism:
#' Quantitative approaches to niche evolution. Evolution, 62, 2868–2883. doi:10.1111/j.1558-5646.2008. 00482.x.
#' 
#' Pianka ER (1973) The Structure of Lizard Communities. Annual Review of Ecology and Systematics, 4, 
#' 53–74. doi:10.1146/annurev.es.04.110173.000413.
#' 
#' @export
#' @importFrom gtools combinations   


nicheOver<-function(stack, metric = c("D", "O")){
  metric <- match.arg(metric, choices = c("D", "O"))
  n.sps <- nlayers(stack)
  ax.names <- names(stack)
  #   D<-matrix(NA,n.sps, n.sps)
  #   I<-matrix(NA,n.sps, n.sps)
  #   P<-matrix(NA,n.sps, n.sps)
  ov <- matrix(1,n.sps, n.sps)
  com <- gtools::combinations(n.sps, 2)
  for (i in 1:nrow(com)){
    p1 <- stack[[com[i,1]]]@data@values
    p2 <- stack[[com[i,2]]]@data@values
    p1 <- p1[-which(is.na(p1))]
    p2 <- p2[-which(is.na(p2))]
    #dp1<-density(p1, kernel="gaussian")#,from=xmin,to=xmax,n=length(p1),cut=0) # calculate the density of occurrences in a vector of R pixels along the score gradient
    # using a gaussian kernel density function, with R bins.
    #dp2<-density(p2, kernel="gaussian")#,from=xmin,to=xmax,n=length(p2),cut=0) 
    
    #x1<-dp1$x   										# breaks on score gradient
    
    #x2<-dp2$x   										# breaks on score gradient
    p1s<-(p1/sum(p1))
    p2s<-(p2/sum(p2))
    #     dp1s<-(dp1$y/sum(dp1$y))
    #     dp2s<-(dp2$y/sum(dp2$y))
    if(metric == "D"){
      ov[com[i,2],com[i,1]] <- 1-(0.5*(sum(abs(p1s-p2s)))) # overlap metric D
    }else if(metric == "O"){
      ov[com[i,2],com[i,1]] <- sum(p1s*p2s)/(sqrt(sum(p1s^2) * sum(p2s^2))) # overlap metric P
    }else if(metric == "I"){
      ov[com[i,2],com[i,1]] <- 1-(0.5*(sqrt(sum((sqrt(p1s)-sqrt(p2s))^2)))) 
    }
    ov[com[i,1],com[i,2]] <- NA
    
    #     D[com[i,1],com[i,2]]<- 1-(0.5*(sum(abs(p1s-p2s)))) # overlap metric D
    #     D[com[i,2],com[i,1]]<- 1-(0.5*(sum(abs(dp1s-dp2s)))) # overlap metric D kernel
    #     
    #     I[com[i,1],com[i,2]] <- 1-(0.5*(sqrt(sum((sqrt(p1s)-sqrt(p2s))^2))))  # overlap metric I
    #     I[com[i,2],com[i,1]] <- 1-(0.5*(sqrt(sum((sqrt(dp1s)-sqrt(dp2s))^2))))  # overlap metric I kernel
    #     
    #     P[com[i,1],com[i,2]] <- sum(p1s*p2s)/(sqrt(sum(p1s^2) * sum(p2s^2))) #overlap metrik P
    #     P[com[i,2],com[i,1]] <- sum(dp1s*dp2s)/(sqrt(sum(dp1s^2) * sum(dp2s^2))) #overlap metrik P kernel
    #     
  }
  
  #   ov[which(is.na(ov) == TRUE)] <- 1
  #   D[which(is.na(D) == TRUE)]<-1
  #   I[which(is.na(I) == TRUE)]<-1
  #   P[which(is.na(P) == TRUE)]<-1
  #   
  dimnames(ov) <- list(ax.names,ax.names)
  attr(ov, "metric") <- metric
  #   dimnames(D)<-list(ax.names,ax.names)
  #   dimnames(I)<-list(ax.names,ax.names)
  #   dimnames(P)<-list(ax.names,ax.names)
  #   
  #   overlap <- list("D"=D, "I"=I, "P" = P)
  #   return(overlap)
  return(ov)
}

