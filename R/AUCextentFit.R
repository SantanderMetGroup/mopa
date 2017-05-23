

#' @title Non-linear model fitting for extracting an index of the minimum value in x that obtains a value in y above the asymptote
#' @description Given a matrix 
#' extracts an index of the extents at which the AUC (or other test statistic) overcomes the Vm coefficient 
#' obtained from the non-linear model fitting. Non-linear asymptotic functions are tested:
#' Michaelis-Menten, 3-parameter asymptotic exponential (exponential 3) and 2-parameter 
#' asymptotic exponential (exponential 2). For details about these functions, see the 
#' "non-linear regression" section of THE R BOOK by Michael J. Crawley.
#' 
#' @param testmat Matrix. Data in each row is fitted separatelly.
#' @param extents Vector of increasing values that is fitted to each row in \code{testmat}.
#' @param diagrams Logical (default is FALSE). If TRUE diagrams of the fitted models are returned.
#' 
#' @details This function is internally used by function \code{\link[mopa]{mopaTrain}}, to fit AUC 
#' scores to the background extents at which were obtained (see reference below).
#'   
#' @return Named integer. Index of the minimum value in \code{extents} at which the data in \code{testmat}
#' is above the asymptote. The residual sum of squares is printed in the console and 
#' the index corresponding to the best nls fit is returned.
#' 
#' 
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
#' @export
#' @import lattice
#' @importFrom stats nls coef na.omit 


AUCextentFit <- function (testmat, extents, diagrams = FALSE) {
  e <- list()
  dat<-list()
  if(is.null(rownames(testmat))) rownames(testmat) <- paste("species", 1:nrow(testmat))
  message(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
  
  for (i in 1:nrow(testmat)) {
    print(i)
    
    y <- na.omit(testmat[i,])
    x <- extents[i,][1:length(y)]
    
    
    group <- rep(rownames(testmat)[i], length(y))
    
    
    out <- nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), algorithm = "plinear")
    # simpler model to get better starting values
    
    a <- coef(out)[2]
    b <- coef(out)[2]-0.5
    c <- coef(out)[1]
    
    
    stmic <- list(Vm=max(predict(nls(y ~ 1*(1-exp(-c * x)), start = list(c=0.001), algorithm = "plinear"))), K=min(x))
    micmen <- nls(y~Vm*x/(K+x), start= stmic) 
    
    
    asym <- nls(y ~ a*(1-exp(-c * x)), start = list(a = a, c=c))
    
    
    asym3 <- tryCatch({nls(y ~ a - b*exp(-c * x), start = list(a = a, b=b, c=c))}, 
                      error = function(err){NULL})
    
    resid<-numeric(length = 3)
    names(resid)<-c("Michaelis Menten", "exponential3", "exponential2")
    resid[1] <- sum((predict(micmen)-y)^2)
    resid[2] <- tryCatch({sum((predict(asym3)-y)^2)}, error = function(err){NA})
    resid[3] <- sum((predict(asym)-y)^2)
    resid <- sort(resid, decreasing = F)
    message (paste ("residual sum of squares for species ", rownames(testmat)[i]))
    mess <- (paste(names(resid),"=", resid))
    message(mess[1])
    message(mess[2])
    if(is.null(asym3)){
      message("Fitting of exponential3 cannot be done")
    }else{
      message(mess[3])
    }
    nls.fun <- names(resid)[1]
    message(paste("best function = ", nls.fun))
    
    
    
    dat[[i]]<- data.frame(x,y, group)
    
    #--------------------------------------------------------------------------------------------
    coefs <- c("Michaelis Menten" = unname(coef(micmen)[1]), 
               "exponential3" = unname(coef(asym3)[1]), 
               "exponential2" = unname(coef(asym)[1]))
    rep <- 1
    if(any(testmat[i, ] > min(coefs))){
      while(!any(testmat[i, ] > coefs[which(names(coefs) == nls.fun)], na.rm = TRUE)){
        rep <- rep + 1
        nls.fun <- names(resid)[rep]
      }
      if(rep != 1){
        message("IMPORTANT WARNING: The coefficient of the asymptote is not outperformed for function/s ", 
                names(resid)[1:(rep-1)], ". Index corresponding to", names(resid)[rep], " is returned.")
      } 
    }else{
        nls.fun <- "max.value"
        message("IMPORTANT WARNING: The coefficient of the asymptote is not outperformed for any of the 
                non linear functions implemented. Index corresponding to the maximum AUC value is returned")
    }
    
    
    if (nls.fun == "Michaelis Menten"){
      a <- coef(micmen)[1]
    }else if(nls.fun == "exponential3"){
      a <- coef(asym3)[1] 
    }else if(nls.fun == "exponential2"){
      a<-coef(asym)[1] 
    }
    qs <- which(na.omit(testmat[i, ]) > a)
    e[[i]] <- qs[which(qs == min(qs))]
    
    
    message(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    
    ########################################################################-
    
    if (nls.fun == "max.value"){e[[i]] <- which(y == max(y))}
    
  } #loop end
  
  #---------------diagrams------------------------------------------------------------------------------------------
  
  if (diagrams == TRUE){
    
    dfa <- do.call("rbind", dat)
    
    
    pl<-xyplot( y ~ x|group, data = dfa, ylab="AUC", xlab="Background radius", 
                key = list(text = list(lab=c("Michaelis Menten", "exponential 3", "exponential 2")), lines = list(lwd=c(1,1,1), span=0.5, col=c("red", "black", "green")),
                           columns =1),
                
                
                panel=function(x, y){
                  ## add lines to the graph which denote means of x and y
                  panel.xyplot(x,y)
                  
                  
                  
                  ##############################exp2#################################------
                  
                  panel.lines(x, predict(nls(y ~ a * (1 - exp(-c * x)), 
                                             start = list(a = coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), algorithm = "plinear"))[2], 
                                                          
                                                          c=coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), algorithm = "plinear"))[1]))),
                              col="green", lwd = 1.5)
                  
                  
                  panel.abline(coef(nls(y ~ a * (1 - exp(-c * x)), 
                                        start = list(a = coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), algorithm = "plinear"))[2], 
                                                     
                                                     c=coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), algorithm = "plinear"))[1])))[1],
                               lty = 5, col = "green", lwd = 1.5)
                  
                  
                  panel.abline(v = min(x[which(y > coef(nls(y ~ a * (1 - exp(-c * x)), 
                                                            start = list(a = coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), algorithm = "plinear"))[2], 
                                                                         
                                                                         c=coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), algorithm = "plinear"))[1])))[1])]),
                               lty = 5, col = "green", lwd = 1.5)
                  
                  ###############################micmen##################################-----
                  
                  panel.lines(x, predict(nls(y ~ Vm * x/(K + x), start=
                                               list(Vm = max(predict(nls(y ~ 1 * (1 - exp(-c * x)), start=list(c=0.001), algorithm ="plinear"))), 
                                                    K = min(x)))),
                              col="red", lwd = 1.5)
                  
                  
                  panel.abline(coef(nls(y ~ Vm * x/(K + x), start=
                                          list(Vm = max(predict(nls(y ~ 1 * (1 - exp(-c * x)), start=list(c=0.001), algorithm ="plinear"))), 
                                               K = min(x))))[1],
                               lty = 5, col = "red", lwd = 1.5)
                  
                  
                  panel.abline(v = min(x[which(y > coef(nls(y ~ Vm * x/(K + x), start=
                                                              list(Vm = max(predict(nls(y ~ 1 * (1 - exp(-c * x)), start=list(c=0.001), algorithm ="plinear"))), 
                                                                   K = min(x))))[1])]),
                               lty = 5, col = "red", lwd = 1.5)
                  
                  
                  
                  ##############################exp3#################################------
                  if(!is.null(asym3)){
                    panel.lines(x, predict(nls(y ~ a - b*exp(-c * x), 
                                               start = list(a = coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), algorithm = "plinear"))[2], 
                                                            b=coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), algorithm = "plinear"))[2]-0.5, 
                                                            c=coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), algorithm = "plinear"))[1]))),
                                col="black", lwd = 1.5)
                    
                    
                    panel.abline(coef(nls(y ~ a - b*exp(-c * x), 
                                          start = list(a = coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), algorithm = "plinear"))[2], 
                                                       b=coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), algorithm = "plinear"))[2]-0.5, 
                                                       c=coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), algorithm = "plinear"))[1])))[1],
                                 lty = 5, col = "black", lwd = 1.5)
                    
                    
                    panel.abline(v = min(x[which(y > coef(nls(y ~ a - b*exp(-c * x), 
                                                              start = list(a = coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), algorithm = "plinear"))[2], 
                                                                           b=coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), algorithm = "plinear"))[2]-0.5, 
                                                                           c=coef( nls (y ~ 1*(1-exp(-c * x)), start = list(c=0.001), algorithm = "plinear"))[1])))[1])]),
                                 lty = 5, col = "black", lwd = 1.5)
                    
                  } 
                    
                  
                }
    )
    
    print(pl)
  }
  e <- unlist(e)
  attr(e, "algorithm") <- attr(testmat, "algorithm")
  attr(e, "species") <- attr(testmat, "species")
  attr(e, "extents") <- attr(testmat, "extents")
  attr(e, "source directory") <- attr(testmat, "source directory")
  return(e)
} #end
