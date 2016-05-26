

#' @title Non-linear model fitting for extracting an index of threshold background extents 
#' @description Given a matrix as returned by function loadTextValues
#' extracts an index of the extents at which the AUC (or other test statistic) overcomes the Vm coefficient 
#' obtained from the non-linear model fitting. Non-linear asymptotic functions are tested:
#' Michaelis-Menten, 3-parameter asymptotic exponential (exponential 3) and 2-parameter 
#' asymptotic exponential (exponential 2). For details about these functions, see the 
#' "non-linear regression" section of THE R BOOK by Michael J. Crawley.
#' 
#' @param testmat Object returned by function loadTestValues (matrix with extents 
#' in columns and species in rows)
#' @param diagrams Logical. If TRUE diagrams of the fitted models are returned
#'   
#' @return Named integer. Index of the extents at which the AUC, kappa or tss 
#' statistics obtain the maximum score. The residual sum of squares is printed in the console and 
#' the index corresponding to the best fit is returned.
#' 
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
#' auc_mars <-loadTestValues(data = presaus, "auc", "mars") 
#' ind <- indextent(testmat = auc_mars, diagrams = TRUE)
#' }
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @export
#' @import lattice


indextent <- function (testmat, diagrams = TRUE) {
  e <- list()
  dat<-list()
  
  message(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
  
  for (i in 1:nrow(testmat)) {
    
    y<-na.omit(testmat[i,])
    x<-as.integer(sub("km", x=names(y), replacement= ""))[1:length(y)]
    
    
    group<-rep(rownames(testmat)[i], length(y))
    
    
    out <- nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), alg = "plinear")
    # simpler model to get better starting values
    
    a <- coef(out)[2]
    b <- coef(out)[2]-0.5
    c <- coef(out)[1]
    
    
    stmic <- list(Vm=max(predict(nls(y ~ 1*(1-exp(-c * x)), start = list(c=0.001), alg ="plinear"))), K=min(x))
    micmen <- nls(y~Vm*x/(K+x), start= stmic) 
    
    
    asym <- nls(y ~ a*(1-exp(-c * x)), start = list(a = a, c=c))
    
    
    asym3 <- nls(y ~ a - b*exp(-c * x), start = list(a = a, b=b, c=c))
    
    resid<-numeric(length = 3)
    names(resid)<-c("Michaelis Menten", "exponential3", "exponential2")
    resid[1] <- sum((predict(micmen)-y)^2)
    resid[2] <- sum((predict(asym3)-y)^2)
    resid[3] <- sum((predict(asym)-y)^2)
    resid <- sort(resid, decreasing = F)
    message (paste ("residual sum of squares for species ", rownames(testmat)[i]))
    mess<-(paste(names(resid),"=", resid))
    message(mess[1])
    message(mess[2])
    message(mess[3])
    nls.fun <- names(resid[1])
    message(paste("best function = ", nls.fun))
    
    
    
    dat[[i]]<- data.frame(x,y, group)
    
    #--------------------------------------------------------------------------------------------
    if (nls.fun == "Michaelis Menten"){
      a<-coef(micmen)[1] 
      qs <- which(na.omit(testmat[i, ]) > a)
      
      if (length(qs) == 0){
        rm(a,qs)
        id <- which(names(resid) == "Michaelis Menten" )+1
        
        if (length(id) > length(resid)){
          message("IMPORTANT WARNING: The coefficient of the asymptote is not outperformed for any of the non linear functions implemented. Index corresponding to the maximum AUC value is returned")
          nls.fun <- "max.value"
          diagrams = FALSE
        }else{
          nls.fun <- names(resid[id])
          message(paste("WARNING: in species -->", rownames(testmat)[i],  "<-- the Michaelis Menten coefficient for the maximum achieved by the system is not outperformed, index corresponding to the", 
                        nls.fun, "non-linear function is returned. See ?indextent.", sep=" "))
        }
        
        
        
      }else{ e[[i]] <- qs[which(qs == min(qs))]}
      
    }
    ########################################################################-
    
    if (nls.fun == "exponential3"){
      
      a<-coef(asym3)[1] 
      qs <- which(na.omit(testmat[i, ]) > a)
      
      if (length(qs) == 0){
        
        rm(a,qs)
        id <- which(names(resid) == "exponential3" )+1
        
        if (length(id) > length(resid)){
          message("IMPORTANT WARNING: The coefficient of the asymptote is not outperformed for any of the non linear functions implemented. Index corresponding to the maximum AUC value is returned")
          nls.fun <- "max.value"
          diagrams = FALSE
        }else{
          nls.fun <- names(resid[id])
          message(paste("WARNING: in species -->", rownames(testmat)[i],  "<-- the exponential3 coefficient of the asymptote is not outperformed, index corresponding to", 
                      nls.fun, "non-linear function is returned. See ?indextent.", sep=" "))
        }
        
      }else{
        e[[i]] <- qs[which(qs == min(qs))]}
    }
    
    ########################################################################-
    
    if (nls.fun == "exponential2"){
      
      a<-coef(asym)[1] 
      qs <- which(na.omit(testmat[i, ]) > a)
      
      if (length(qs) == 0){
        
        rm(a,qs)
        id <- which(names(resid) == "exponential2" )+1
        
        if (length(id) > length(resid)){
          message("IMPORTANT WARNING: The coefficient of the asymptote is not outperformed for any of the non linear functions implemented. Index corresponding to the maximum AUC value is returned")
          nls.fun <- "max.value"
          diagrams = FALSE
        }else{
          nls.fun <- names(resid[id])
          message(paste("WARNING: in species -->", rownames(testmat)[i],  "<-- the exponential2 coefficient of the asymptote is not outperformed, index corresponding to", 
                      nls.fun, "non-linear function is returned. See ?indextent.", sep=""))
        }
        
      }else{
        e[[i]] <- qs[which(qs == min(qs))]}
    }
    
    message(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    
    ########################################################################-
    
    if (nls.fun == "max.value"){e[[i]] <- which(y == max(y))}
    
  } #loop end
  
  #---------------diagrams------------------------------------------------------------------------------------------
  
  if (diagrams == TRUE){
    
    dfa<-do.call("rbind", dat)
    
    
    pl<-xyplot( y ~ x|group, data = dfa, ylab="AUC", xlab="Background extent", 
                key = list(text = list(lab=c("Michaelis Menten", "exponential 3", "exponential 2")), lines = list(lwd=c(1,1,1), span=0.5, col=c("red", "black", "green")),
                           columns =1),
                
                
                panel=function(x, y){
                  ## add lines to the graph which denote means of x and y
                  panel.xyplot(x,y)
                  
                  
                  
                  ##############################exp2#################################------
                  
                  panel.lines(x, predict(nls(y ~ a * (1 - exp(-c * x)), 
                                             start = list(a = coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), alg = "plinear"))[2], 
                                                          
                                                          c=coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), alg = "plinear"))[1]))),
                              col="green", lwd = 1.5)
                  
                  
                  panel.abline(coef(nls(y ~ a * (1 - exp(-c * x)), 
                                        start = list(a = coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), alg = "plinear"))[2], 
                                                     
                                                     c=coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), alg = "plinear"))[1])))[1],
                               lty = 5, col = "green", lwd = 1.5)
                  
                  
                  panel.abline(v = min(x[which(y > coef(nls(y ~ a * (1 - exp(-c * x)), 
                                                            start = list(a = coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), alg = "plinear"))[2], 
                                                                         
                                                                         c=coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), alg = "plinear"))[1])))[1])]),
                               lty = 5, col = "green", lwd = 1.5)
                  
                  ###############################micmen##################################-----
                  
                  panel.lines(x, predict(nls(y ~ Vm * x/(K + x), start=
                                               list(Vm = max(predict(nls(y ~ 1 * (1 - exp(-c * x)), start=list(c=0.001), alg ="plinear"))), 
                                                    K = min(x)))),
                              col="red", lwd = 1.5)
                  
                  
                  panel.abline(coef(nls(y ~ Vm * x/(K + x), start=
                                          list(Vm = max(predict(nls(y ~ 1 * (1 - exp(-c * x)), start=list(c=0.001), alg ="plinear"))), 
                                               K = min(x))))[1],
                               lty = 5, col = "red", lwd = 1.5)
                  
                  
                  panel.abline(v = min(x[which(y > coef(nls(y ~ Vm * x/(K + x), start=
                                                              list(Vm = max(predict(nls(y ~ 1 * (1 - exp(-c * x)), start=list(c=0.001), alg ="plinear"))), 
                                                                   K = min(x))))[1])]),
                               lty = 5, col = "red", lwd = 1.5)
                  
                  
                  
                  ##############################exp3#################################------
                  
                  panel.lines(x, predict(nls(y ~ a - b*exp(-c * x), 
                                             start = list(a = coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), alg = "plinear"))[2], 
                                                          b=coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), alg = "plinear"))[2]-0.5, 
                                                          c=coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), alg = "plinear"))[1]))),
                              col="black", lwd = 1.5)
                  
                  
                  panel.abline(coef(nls(y ~ a - b*exp(-c * x), 
                                        start = list(a = coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), alg = "plinear"))[2], 
                                                     b=coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), alg = "plinear"))[2]-0.5, 
                                                     c=coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), alg = "plinear"))[1])))[1],
                               lty = 5, col = "black", lwd = 1.5)
                  
                  
                  panel.abline(v = min(x[which(y > coef(nls(y ~ a - b*exp(-c * x), 
                                                            start = list(a = coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), alg = "plinear"))[2], 
                                                                         b=coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), alg = "plinear"))[2]-0.5, 
                                                                         c=coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), alg = "plinear"))[1])))[1])]),
                               lty = 5, col = "black", lwd = 1.5)
                  
                  
                  
                })
    
    print(pl)
  }
  e <- unlist(e)
  attr(e, "algorithm") <- attr(testmat, "algorithm")
  attr(e, "species") <- attr(testmat, "species")
  attr(e, "extents") <- attr(testmat, "extents")
  attr(e, "source directory") <- attr(testmat, "source directory")
  return(e)
} #end
