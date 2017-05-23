
#' @title Leave out a different subset for test each fold
#' @description Given a data frame splitted in k subsets, this function binds 
#' subsets for training data and leaves one subset out for test k times, leaving 
#' out a different subset each fold. 
#' 
#' 
#' @param kf.data Object returned by function kfold or list of data frames
#' 
#' @return  List of length k. Each slot contains another list of two components:
#'  \item{train }{data frame of training data}
#'  \item{test }{data frame for test}
#' 
#' 
#' @details Details.
#' 
#' 
#' 
#' @author M. Iturbide 
#' 
#' @keywords internal
#'
#' @export
#' @keywords internal
#' @importFrom sampling strata


leaveOneOut<-function(kf.data){
  x<-kf.data
  k<-length(kf.data)
  tt.list <- list()
  for (i in 1:length(x)) {
    test.x <- x[[i]]
    train.x <- data.frame()
    for (j in c(1:k)[-i]) {
      train.x <- rbind.data.frame(train.x, x[[j]])
    }
    length(tt.list) <- i
    tt.list [[i]]<- list("test"=test.x, "train"=train.x)
  }
  return(tt.list)
}
