
#' @title Species distribution modelling and cross validation
#' @description Species distribution modelling with k-fold cross validation. 
#' Algorithms supported are "glm", "svm", "maxent", "mars", "randomForest", "cart.rpart" 
#' and "cart.tree" 
#' @param kdata Object returned by function leaveOneOut
#' @param data Object returned by function biomat. 2D matrix with the dependent variable 
#' (presence/absence) in the first column and the independent variables in the rest 
#' (extracted from varstack) 
#' @param algorithm Any character of the following: \code{"glm"}, "svm", "maxent", "mars", "randomForest", "cart.rpart" 
#' or "cart.tree"
#' @param threshold Cut value between 0 and 1 to calculate the confusion matrix. Default is 0.5.
#'
#'  
#' 
#' @return  A list with six components:
#'  \item{allmod }{fitted model using all data for training}
#'  \item{auc }{AUC statistic in the cross validation}
#'  \item{kappa }{kappa statistic in the cross validation}
#'  \item{tss }{true skill statistic in the cross validation }
#'  \item{mod }{fitted model with partitioned data}
#'  \item{p }{cross model prediction}
#' 
#' @details This function calculates the AUC with the function "auc" from package 
#' "PresenceAbsence". Package SDMTools must be detached.
#' 
#' 
#' 
#' @author M. Iturbide 
#' 
#'
#' @export
#' @keywords internal
#' @import PresenceAbsence
#' @importFrom e1071 svm
#' @importFrom sampling strata  
#' @importFrom e1071 best.svm 
#' @importFrom dismo maxent
#' @importFrom earth earth
#' @importFrom rpart rpart
#' @importFrom tree tree
#' @importFrom randomForest randomForest



modelo<-function(kdata, data, algorithm = c("glm","svm","maxent","mars","randomForest","cart.rpart","cart.tree"), threshold = 0.5){
  mod <- list()
  pmod <- list()
  algorithm <- as.character(algorithm)
  orig<-list()
  alltrn<-na.omit(data)
  colnames(alltrn)<-c("V1", colnames(alltrn)[2:ncol(alltrn)])
                      
  for(i in 1:length(kdata)) {
    length(mod) <- i
    length(pmod) <- i
    dftst <- na.omit(kdata[[i]]$test)
      colnames(dftst)<-c("V1", colnames(dftst)[2:ncol(dftst)])   
    dftrn <- na.omit(kdata[[i]]$train)
      colnames(dftrn)<-c("V1", colnames(dftrn)[2:ncol(dftrn)])  
    orig[[i]]<-na.omit(dftst[,1])
    
    # Training
    if(algorithm == "glm") {
      mod[[i]] <- glm(V1 ~., data=dftrn, family=binomial(link="logit"))
    } else if(algorithm == "svm"){
      mod[[i]] <- best.svm(V1 ~., data=dftrn)
    } else if(algorithm == "maxent"){
      mod[[i]] <- maxent(x = dftrn[ ,-1], p=dftrn[,1])
    } else if(algorithm == "mars"){
      mod[[i]] <- earth(x = dftrn[ ,-1], y = dftrn[,1])
    } else if(algorithm=="cart.rpart"){
      mod[[i]] <- rpart(V1 ~., data=dftrn)
    } else if(algorithm == "cart.tree"){
      mod[[i]]<- tree(V1 ~., data=dftrn)
    } else if(algorithm == "randomForest"){
      mod[[i]]<- randomForest(V1 ~., data=dftrn)
    }
    # Predictions
    if (algorithm == "cart.rpart") {
      pmod[[i]] <- predict(mod[[i]], dftst[,-1])
    }else if (algorithm=="cart.tree"){
      pmod[[i]] <- predict(mod[[i]], dftst[,-1])
    }else{
      pmod[[i]] <- predict(mod[[i]], dftst[ ,-1], type="response")
    }
    
  }
  
  
  # allTraining
  if(algorithm == "glm") {
    allmod <- glm(V1 ~., data=alltrn, family=binomial(link="logit"))
  } else if(algorithm == "svm"){
    allmod <- best.svm(V1 ~., data=alltrn)
  } else if(algorithm == "maxent"){
    allmod <- maxent(x = alltrn[ ,-1], p=alltrn[,1])
  } else if(algorithm == "mars"){
    allmod <- earth(x = alltrn[ ,-1], y = alltrn[,1])
  } else if(algorithm=="cart.rpart"){
    allmod <- rpart(V1 ~., data=alltrn)
  } else if(algorithm == "cart.tree"){
    allmod <- tree(V1 ~., data=alltrn)
  } else if(algorithm == "randomForest"){
    allmod <- randomForest(V1 ~., data=alltrn)
  }
  
  
  ori<-unlist(orig)
  p<-unlist(pmod)
  p.auc <- p
  p.auc[which(p>=threshold)] <- 1
  p.auc[which(p.auc!=1)] <- 0
  dat<-as.data.frame(cbind(as.integer(c(1:length(ori))), ori, p))
  dat.auc <- as.data.frame(cbind(as.integer(c(1:length(ori))), ori, p.auc))
  mod.auc <- auc(dat.auc, st.dev = FALSE)
  mod.kappa <- Kappa(cmx(dat, threshold = threshold), st.dev = FALSE)
  mod.tss <- sensitivity(cmx(dat, threshold = threshold),st.dev = FALSE) + 
    specificity (cmx(dat, threshold = threshold),st.dev = FALSE) - 1
  rm(ori, p)
  return(list("allmod"=allmod, "auc" = mod.auc, "kappa"=mod.kappa, "tss"= mod.tss ,"mod" = mod,"p" = dat)) # behar dira ere modeloak
}
