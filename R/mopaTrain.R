#' @title Easy species distribution modeling and cross validation
#' @description Species distribution modeling and k-fold cross validation 
#' for a set of presence/absence data per species, also considering different background 
#' extents (optional). Algorithms supported are "glm", "svm", "maxent", "mars", "rf", "cart.rpart" 
#' and "cart.tree" 
#' 
#' @param y Object returned by function \code{\link[mopa]{pseudoAbsences}} or data frame or list/s of data 
#' frames with coordinates
#'  in the first two columns and presence/absence (1=presence, 0=absence) in the third column. 
#' @param x RasterStack ot list of RasterStacks of variables for modeling, a.k.a baseline environment/climatology
#' @param k Integer. Number of folds for cross validation. Default is 10
#' @param algorithm Character string of the algorithms for modeling. Options are the following: 
#' "glm", "svm", "maxent", "mars", "rf", "cart.rpart" and "cart.tree" (see details)
#' @param algorithm.args Further arguments to be passed to the selected algorithm for modeling (functions involved 
#' are described in details). 
#' @param weighting Logical for model fitting with weighted presence/absences. Applicable for algorithms "glm", "mars", 
#' "rf", cart.tree and "cart.rpart". Default is FALSE.
#' The processing time is considerably increased if weighting option is selected when the  
#' "mars" algorithm (see \code{\link[earth]{earth}} is applied.
#' @param threshold Cut value between 0 and 1 to calculate the confusion matrix. Default is NULL (see Details).
#' @param diagrams Logical. Only applied if \code{x} contains data for different background extents 
#' (see \code{\link[mopa]{backgroundRadius}} and \code{\link[mopa]{pseudoAbsences}}). Should diagrams of 
#' AUC extent fitting be printed? default is FALSE. 
#' @param tuneRF.args list of arguments from function \code{\link[randomForest]{tuneRF}}. Only used when algorihm = "rf"
#' 
#'  
#' 
#' @return A list of six components is returned for each species in \code{x}:
#' \itemize{
#'  \item{model }{fitted model using all data for training}
#'  \item{auc }{AUC statistic in the cross validation}
#'  \item{kappa }{kappa statistic in the cross validation}
#'  \item{tss }{true skill statistic in the cross validation }
#'  \item{fold.models }{fitted models of each data partition for cross validation}
#'  \item{ObsPred }{cross model prediction (e.g. for further assessment of model accuracy)}
#'  }
#' 
#' @details This function calculates the AUC with the function \code{\link[PresenceAbsence]{auc}} from package 
#' \pkg{PresenceAbsence}. \strong{Note:} Package \pkg{SDMTools} must be detached. 
#' 
#' If \code{threshold} is not specified the value that maximisez the TSS (true skill statistic) is 
#' used to calculate the confusion matrix.
#' 
#' 
#' 
#' If \code{y} contains data for different background extents (see \code{\link[mopa]{backgroundRadius}} and
#' \code{\link[mopa]{pseudoAbsences}}), \code{\link[mopa]{mopaTrain}} performs the species distribution modeling for 
#' each different background extent, and fits obtained AUCs (corresponding to different background extents) 
#' to three non linear models (Michaelis-Menten, exponential2 and exponential3). 
#' The model that scores the lowest error is automatically selected to extract the Vm coefficient (equation 1 in 
#' Iturbide et al., 2015). Then, the minimum extent at which the AUC surpasses the Vm value is selected 
#' as the threshold extent (see Figure 3 in Iturbide et al., 2015), being the corresponding fitted SDM the 
#' one returned by \code{mopaFitting}. If argument \code{diagrams} is set to TRUE, A fitted model plot 
#' (as in Fig. 3 in Iturbide et al., 2015) is printed in the plotting environment.
#'
#'   
#' \code{mopaTrain} uses the algorithm implementations of the following functions and R packages: 
#' \itemize{
#' \item{"mars" }{function \code{\link[earth]{earth}} from package \pkg{earth}} 
#' \item{"rf" }{function \code{\link[ranger]{ranger}} from package \pkg{ranger}} 
#' \item{"maxent" }{function \code{\link[dismo]{maxent}} from package \pkg{dismo}} 
#' \item{"cart.rpart" }{function \code{\link[rpart]{rpart}} from package \pkg{rpart}}
#' \item{"svm" }{function \code{\link[e1071]{best.svm}} from package \pkg{e1071}}
#' \item{"cart.tree" }{function \code{\link[tree]{tree}} from package \pkg{tree}} 
#' \item{"glm" }{function \code{\link[stats]{glm}} from  package \pkg{stats}}
#' }
#' 
#' For example, when appying "glm", further arguments from function \code{\link[stats]{glm}} can be 
#' passed to \code{mopaTrain} by using \code{algorithm.args}.
#' 
#' @seealso \code{\link[mopa]{mopaPredict}}, \code{\link[mopa]{pseudoAbsences}}, \code{\link[mopa]{backgroundGrid}}, 
#' \code{\link[mopa]{OCSVMprofiling}}, \code{\link[mopa]{backgroundRadius}}, \code{\link[mopa]{extractFromModel}}
#' @author M. Iturbide 
#' 
#' @examples
#' data(Oak_phylo2)
#' data(biostack)
#' r <- biostack$baseline[[1]]
#' ## Create background grid
#' bg <- backgroundGrid(r)
#' ## Generate pseudo-absences
#' RS_random <-pseudoAbsences(xy = Oak_phylo2, background = bg$xy, 
#'                            exclusion.buffer = 0.083*5, prevalence = -0.5, kmeans = FALSE)
#' ## Model training
#' fittedRS <- mopaTrain(y = RS_random, x = biostack$baseline, 
#'                       k = 10, algorithm = "glm", weighting = TRUE)
#' ## Extract fitted models
#' mods <- extractFromModel(models = fittedRS, value = "model")
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @export
#' 
#' @import raster
#' @import sp
#' @importFrom sampling strata  


mopaTrain <- function(y, 
                      x, 
                      k = 10, 
                      algorithm = c("glm", "svm", "maxent", "mars", "rf", "cart.rpart", "cart.tree"), 
                      algorithm.args = NULL,
                      weighting = FALSE,
                      threshold = NULL,
                      diagrams = FALSE,
                      tuneRF.args = NULL){
  if(class(x) != "list") x <- list(x)
  smfit <- list()
  for(j in 1:length(y)){
    message("[", Sys.time(), "] Modeling species ", j)
    yy <- y[[j]]
    mfit <- list()
    for(i in 1:length(yy)){
      message(":::[", Sys.time(), "] realization ", i)
      fit <- list()
        for(h in 1:length(algorithm)){
          fitbase <- list()
          nm <- character()
          for(l in 1:length(x)){
          fitbase[[l]] <- mopaTrain0(y = yy[[i]], 
                              x = x[[l]], 
                              k = k, 
                              algorithm = algorithm[h], 
                              algorithm.args = algorithm.args,
                              weighting = weighting,
                              threshold = threshold,
                              diagrams = diagrams,
                              tuneRF.args = tuneRF.args,
                              plotnames = paste0(names(y)[j], " PAralization ", i))
          if(l < 10){
            nm[l] <- paste0("0", l)
          }else{
            nm[l] <- as.character(l)
          }
          }
          if(is.null(names(x))) names(x) <- paste0("baseClim", nm)
          names(fitbase) <- names(x)
          fit[[h]] <- fitbase
    }
      names(fit) <- algorithm
      mfit[[i]] <- fit
    }
    names(mfit) <- names(yy)
    smfit[[j]] <- mfit
  }
  names(smfit) <- names(y)
  return(smfit)
}


#end

#' @title Easy species distribution modeling and cross validation
#' @description Species distribution modeling and k-fold cross validation 
#' for a set of presence/absence data per species, also considering different background 
#' extents (optional). Algorithms supported are "glm", "svm", "maxent", "mars", "rf", "cart.rpart" 
#' and "cart.tree" 
#' 
#' @param x Object returned by function \code{link[mopa]{pseudoAbsences}} or list/s of data frames with coordinates
#'  in the first two columns and presence/absence (1=presence, 0=absence) in the third column. 
#' @param y RasterStack of variables for modelling
#' 
#' @param k Integer. Number of folds for cross validation. Default is 10
#' @param algorithm Any character of the following: "glm", "svm", "maxent", "mars", "rf", "cart.rpart" 
#' or "cart.tree"
#' @param algorithm.args Further arguments to be passed to the selected algorithm for modeling (functions involved 
#' are described in details)
#' @param weighting Logical for "glm", "mars" and "rf" fitting with weighted presence/absences-s. Default is FALSE.
#' @param threshold Cut value between 0 and 1 to calculate the confusion matrix. Default is NULL (see Details).
#' @param diagrams logical. Only applied if \code{x} contains data for different background extents 
#' (see \code{\link[mopa]{backgroundRadius}} and \code{\link[mopa]{pseudoAbsences}}). Should diagrams of 
#' AUC extent fitting be printed? default is FALSE. 
#' @param plotnames names to be printed in the diagrams
#' @param tuneRF.args list of arguments from function \code{\link[randomForest]{tuneRF}}. Only used when algorihm = "rf"
#' 
#'  
#' 
#' @return A list of six components is returned for each species in \code{x}:
#' 
#'  \item{model }{fitted model using all data for training}
#'  \item{auc }{AUC statistic in the cross validation}
#'  \item{kappa }{kappa statistic in the cross validation}
#'  \item{tss }{true skill statistic in the cross validation }
#'  \item{fold.models }{fitted model with partitioned data}
#'  \item{ObsPred }{cross model prediction}
#'  
#' 
#' 
#' 
#' 
#' @seealso \code{\link[mopa]{mopaPredict}}, \code{\link[mopa]{pseudoAbsences}}, \code{\link[mopa]{backgroundGrid}}, 
#' \code{\link[mopa]{OCSVMprofiling}}, \code{\link[mopa]{backgroundRadius}}
#' 
#' @author M. Iturbide 
#' 
#' 
#' @references Iturbide, M., Bedia, J., Herrera, S., del Hierro, O., Pinto, M., Gutierrez, J.M., 2015. 
#' A framework for species distribution modelling with improved pseudo-absence generation. Ecological 
#' Modelling. DOI:10.1016/j.ecolmodel.2015.05.018.
#' 
#' @keywords internal
#' 
#' @import raster
#' @import sp
#' @importFrom sampling strata  



mopaTrain0 <- function(y, 
                        x, 
                        k = 10, 
                        algorithm = c("glm", "svm", "maxent", "mars", "rf", "cart.rpart", "cart.tree"), 
                        algorithm.args = NULL,
                        weighting = FALSE,
                        threshold = NULL,
                        diagrams = FALSE,
                        tuneRF.args = NULL,
                        plotnames = "unnamed"){
  algorithm <- match.arg(algorithm, choices = c("glm", "svm", "maxent", "mars", "rf", "cart.rpart", "cart.tree"))
  data <- y
  biostack <- x
  if (class(data[[1]]) != "list"){
    data <- list(data)
    names(data) <- plotnames
  }
  extents <- array(dim = c(length(data), max(unlist(lapply(data, length)))), dimnames = list(c(names(data))))
  aucmat <- array(dim = c(length(data), max(unlist(lapply(data, length)))), dimnames = list(c(names(data))))
  for (i in 1:length(data)){
    sp_01 <- data[[i]]
    if(is.null(names(sp_01))) names(sp_01) <- "maximum extent"
    extnames <- names(sp_01)
    extents[i, 1:length(extnames)] <- as.integer(sub(unlist(extnames), pattern = "km", replacement = ""))
    for(j in 1:length(sp_01)){
      #print(paste("running model for species", i, "considering pseudo-absences inside the extent of", names (sp_01)[j]))
      sp.bio <- biomat(sp_01[[j]], biostack)
      x <- kfold(k, df = sp.bio)
      xx <- leaveOneOut(x)
      # mod <- tryCatch({modelo(kdata = xx, data=sp.bio, algorithm = algorithm, weighting = weighting, threshold = threshold)},
      #                 error = function(err){xxx = list(rep(NA, k), NA, NA)})
      aucmat[i, j] <- modelo(kdata = xx, data=sp.bio, algorithm = algorithm, algorithm.args = algorithm.args,
                             weighting = weighting, threshold = threshold, tuneRF.args = tuneRF.args)$auc
      rm(xx, x, sp.bio)
    }
  }
  ind <- rep(1, length(data))
  if(ncol(aucmat) > 1)  ind <- AUCextentFit(aucmat, extents = extents, diagrams = diagrams)
  mod <- list()
  for (i in 1:length(data)){
    sp_01 <- data[[i]]
    sp.bio <- biomat(sp_01[[ind[i]]], biostack)
    x <- kfold(k, df = sp.bio)
    xx <- leaveOneOut(x)
    # mod <- tryCatch({modelo(kdata = xx, data=sp.bio, algorithm = algorithm, weighting = weighting, threshold = threshold)},
    #                 error = function(err){xxx = list(rep(NA, k), NA, NA)})
    mod[[i]] <- modelo(kdata = xx, data=sp.bio, algorithm = algorithm, algorithm.args = algorithm.args,
                       weighting = weighting, threshold = threshold, tuneRF.args = tuneRF.args)
    mod[[i]]$extent <- extents[[i]][ind[i]]
    rm(xx, x, sp.bio)
  }
  if(length(mod) == 1) mod <- mod[[1]]
  return(mod)
}

#end




#' @title Species distribution modeling and cross validation
#' @description Species distribution modeling with k-fold cross validation. 
#' Algorithms supported are "glm", "svm", "maxent", "mars", "rf", "cart.rpart" 
#' and "cart.tree" 
#' @param kdata Object returned by function leaveOneOut
#' @param data Object returned by function biomat. 2D matrix with the dependent variable 
#' (presence/absence) in the first column and the independent variables in the rest 
#' (extracted from varstack) 
#' @param algorithm Any character of the following: \code{"glm"}, "svm", "maxent", "mars", "rf", "cart.rpart" 
#' or "cart.tree"
#' @param algorithm.args Further arguments to be passed to the selected algorithm for modeling (functions involved 
#' are described in details)
#' @param weighting Logical for model fitting with weighted presence/absences-s. Applicable for algorithms "glm", "mars", 
#' "rf" and "cart.rpart". Default is FALSE.
#' @param threshold Cut value between 0 and 1 to calculate the confusion matrix. Default is 0.5.
#' @param tuneRF.args list of arguments from function \code{\link[randomForest]{tuneRF}}. Only used when algorihm = "rf"
#'
#'  
#' 
#' @return  A list with six components:
#'  \item{model }{fitted model using all data for training}
#'  \item{auc }{AUC statistic in the cross validation}
#'  \item{kappa }{kappa statistic in the cross validation}
#'  \item{tss }{true skill statistic in the cross validation }
#'  \item{fold.models }{fitted model with partitioned data}
#'  \item{ObsPred }{cross model prediction}
#' 
#' @details This function calculates the AUC with the function "auc" from package 
#' "PresenceAbsence". Package SDMTools must be detached.
#' 
#' 
#' 
#' @author M. Iturbide 
#' 
#'
#' @keywords internal
#' @import PresenceAbsence
#' @importFrom e1071 svm
#' @importFrom sampling strata  
#' @importFrom e1071 best.svm 
#' @importFrom dismo maxent
#' @importFrom earth earth
#' @importFrom rpart rpart
#' @importFrom tree tree
#' @importFrom randomForest randomForest tuneRF
#' @importFrom ranger ranger
#' @importFrom stats glm binomial na.omit 



modelo <- function(kdata, data, 
                   algorithm = c("glm","svm","maxent","mars","rf","cart.rpart","cart.tree"), 
                   algorithm.args = NULL,
                   weighting = FALSE, threshold = NULL,
                   tuneRF.args = NULL){
  mod <- list()
  pmod <- list()
  algorithm <- as.character(algorithm)
  orig<-list()
  alltrn<-na.omit(data)
  colnames(alltrn)<-c("V1", colnames(alltrn)[2:ncol(alltrn)])
  if(weighting){
    Force.weights <- TRUE
    w.factor <- round(length(which(alltrn$V1 == 0))/length(which(alltrn$V1 == 1))*10)
    all.weights <- alltrn$V1*w.factor
    all.weights[(which(alltrn$V1 == 0))] <- 10
    all.weights.rf <- rep(length(which(alltrn$V1 == 1)), 2)
  }else{
    Force.weights <- FALSE
    all.weights <- rep(10, length(alltrn$V1))
    all.weights.rf <- nrow(alltrn)
  }
  nm <- character()
  for(i in 1:length(kdata)) {
    if(i < 10){
      nm[i] <- paste0("0", i)
    }else{
      nm[i] <- as.character(i)
    }
    length(mod) <- i
    length(pmod) <- i
    dftst <- na.omit(kdata[[i]]$test)
    colnames(dftst)<-c("V1", colnames(dftst)[2:ncol(dftst)])   
    dftrn <- na.omit(kdata[[i]]$train)
    colnames(dftrn)<-c("V1", colnames(dftrn)[2:ncol(dftrn)])  
    orig[[i]]<-na.omit(dftst[,1])
    if(weighting){
      k.w.factor <- round(length(which(dftrn$V1 == 0))/length(which(dftrn$V1 == 1))*10)
      weights <- dftrn$V1*k.w.factor
      weights[(which(dftrn$V1 == 0))] <- 10
      pres <- length(which(dftrn$V1 == 1))
      aus <- length(which(dftrn$V1 == 0))
      if(pres>aus){
        weights.rf <- rep(aus, 2)
      }else{
        weights.rf <- rep(pres, 2) 
      }
    }else{
      weights <- rep(10, length(dftrn$V1))
      weights.rf <- nrow(dftrn)
    }
    # Training
    if(algorithm == "glm") {
      mod[[i]] <- do.call(glm, c(list(formula = V1 ~., data = dftrn, family = binomial(link="logit"), weights = weights), algorithm.args))
    } else if(algorithm == "svm"){
      mod[[i]] <- do.call(best.svm, c(list(x = V1 ~., data = dftrn), algorithm.args))
    } else if(algorithm == "maxent"){
      mod[[i]] <- do.call(maxent, c(list(x = dftrn[ ,-1], p=dftrn[,1]), algorithm.args))
    } else if(algorithm == "mars"){
      mod[[i]] <- do.call(earth, c(list(x = dftrn[ ,-1], y = dftrn[,1], weights = weights), algorithm.args))
    } else if(algorithm=="cart.rpart"){
      mod[[i]] <- do.call(rpart, c(list(V1 ~., data=dftrn, weights = weights), algorithm.args))
    } else if(algorithm == "cart.tree"){
      mod[[i]]<- do.call(tree, c(list(V1 ~., data=dftrn, weights = weights), algorithm.args))
    } else if(algorithm == "rf"){
      # dftrn$V1 <- as.factor(dftrn$V1)
      # mod[[i]]<- randomForest(V1 ~., data=dftrn, strata = dftrn$V1, sampsize = weights.rf)
      if("mtry" %in% names(algorithm.args) == FALSE){
        setbefore <- TRUE
        suppressMessages(
          mtry1 <- do.call(tuneRF, c(list(dftrn[,-1], dftrn[,1], plot = F), tuneRF.args))
        )
        mtry <- mtry1[which(mtry1[,2] == min(mtry1[,2])),1]
        algorithm.args$mtry <- mtry
      }
      mod[[i]]<- do.call(ranger, c(list(V1 ~., data=dftrn, case.weights = weights), algorithm.args))
    }
    # Predictions
    if (algorithm == "cart.rpart") {
      pmod[[i]] <- predict(mod[[i]], dftst[,-1])
    }else if (algorithm=="cart.tree"){
      pmod[[i]] <- predict(mod[[i]], dftst[,-1])
    }else if(algorithm == "rf"){
      pmod[[i]] <- predict(mod[[i]], dftst[ ,-1])$predictions
    }else{
      pmod[[i]] <- predict(mod[[i]], dftst[ ,-1], type="response")
    }
    
  }
  
  
  # allTraining
  if(algorithm == "glm") {
    allmod <- do.call(glm, c(list(V1 ~., data=alltrn, family=binomial(link="logit"), weights = all.weights), algorithm.args))
  } else if(algorithm == "svm"){
    allmod <- do.call(best.svm, c(list(V1 ~., data=alltrn), algorithm.args))
  } else if(algorithm == "maxent"){
    allmod <- do.call(maxent, c(list(x = alltrn[ ,-1], p=alltrn[,1]), algorithm.args))
  } else if(algorithm == "mars"){
    allmod <- do.call(earth, c(list(x = alltrn[ ,-1], y = alltrn[,1], weights = all.weights), algorithm.args))
  } else if(algorithm=="cart.rpart"){
    allmod <- do.call(rpart, c(list(V1 ~., data=alltrn, weights = all.weights), algorithm.args))
  } else if(algorithm == "cart.tree"){
    allmod <- do.call(tree, c(list(V1 ~., data=alltrn, weights = all.weights), algorithm.args))
  } else if(algorithm == "rf"){
    # alltrn$V1 <- as.factor(alltrn$V1)
    # allmod <- randomForest(V1 ~., data=alltrn, strata = alltrn$V1, sampsize = all.weights.rf)
    if(setbefore){
      message(":::Applying tuneRF to obtain parameter mytr")
      suppressMessages(
        mtry1 <- tuneRF(alltrn[,-1], alltrn[,1], plot = F)
      )
      mtry <- mtry1[which(mtry1[,2] == min(mtry1[,2])),1]
      algorithm.args$mtry <- mtry
    }
    allmod <- do.call(ranger, c(list(V1 ~., data=alltrn, case.weights = all.weights), algorithm.args))
  }
  
  
  ori <- unlist(orig)
  p <- unlist(pmod)
  dat <- as.data.frame(cbind("id" = as.integer(c(1:length(ori))), "obs" = ori, "pred" = p))
  if(is.null(threshold)) threshold <- cutTSS(Obs = dat$obs, Fit = dat$pred)$CutOff 
  mod.auc <- auc(dat, st.dev = FALSE)
  mod.kappa <- Kappa(cmx(dat, threshold = threshold), st.dev = FALSE)
  mod.tss <- sensitivity(cmx(dat, threshold = threshold),st.dev = FALSE) + 
    specificity (cmx(dat, threshold = threshold),st.dev = FALSE) - 1
  rm(ori, p)
  names(mod) <- paste0("fold", nm)
  return(list("model"=allmod, "auc" = mod.auc, "kappa"=mod.kappa, "tss"= mod.tss ,"fold.models" = mod, "ObsPred" = dat)) # behar dira ere modeloak
}

#end




#' @title Cut value of the max TSS
#' @description Cut value of the max TSS
#' @param Obs Observed values
#' @param Fit fitted values
#'    
#' @return threshold of the max TSS
#' @keywords internal


cutTSS <-  function(Obs, Fit){
  if(sum(Obs)==0) stop("\n The observed data only contains 0")
  tab <- as.data.frame(matrix(0, nrow=101, ncol=2))  ### 
  
  if(length(unique(Fit))==1){
    Misc<-table(as.vector(Fit) >= as.numeric(unique(Fit)), Obs) 
    a <- TSS.Stat(Misc)
    TP <- Misc[4]
    TN <- Misc[1]
    ca0 <- (TN * 100)/sum(Misc[,1])
    ca1 <- (TP * 100)/sum(Misc[,2])
    if(is.na(ca0)) ca0<-0
    if(is.na(ca1)) ca1<-0
    invisible(list(TSS=0, CutOff=unique(Fit), TP=TP, se=ca1, TN=TN, sp=ca0))
  }
  else{
    Quant <- quantile(Fit)
    for(j in 0:100){
      Seuil <- Quant[1] + (j*((Quant[5] - Quant[1])/100))
      Misc<-table(Fit >= Seuil, Obs)
      a <- TSS.Stat(Misc)
      if(!is.na(a)) if(a > 0) {tab[j+1, 1] <- Seuil; tab[j+1, 2] <- a}
      rm(Misc, Seuil)
    }
    
    t <- max(tab[,2],na.rm=TRUE)
    seuil <- tab[tab[,2]==t,1]   
    if(t > 0) {
      Misc<-table(Fit >= seuil[1], Obs)
      TP <- Misc[4]
      TN <- Misc[1]
      ca0 <- (TN * 100)/sum(Misc[,1])
      ca1 <- (TP * 100)/sum(Misc[,2])
      invisible(list(TSS = t, CutOff = seuil[1], TP = TP, se = ca1, TN = TN, sp = ca0))
    }
    else {
      invisible(list(TSS = 0, CutOff = 0, TP = 0, se = 0, TN = 0, sp = 0))
    }
  }
}

#end



#' @title Internarl function for cutTSS 
#' @description  Internarl function for cutTSS
#' @param data data generated in function cutTSS
#' @keywords internal


TSS.Stat <-
  function(data)
  {
    if(dim(data)[1]==1){
      if(row.names(data)[1]=="FALSE") data <- rbind(data, c(0,0))
      else {
        a<-data
        data<-c(0,0)
        data<-rbind(data, a)
        
      }
    }
    n <- sum(data)
    a <- data[1,1]
    b <- data[1,2]
    c <- data[2,1]
    d <- data[2,2]
    sens <- a/(a+c)
    spec <- d/(b+d)
    K <- (sens + spec) - 1        #TSS
    return(K)
  }

#end
