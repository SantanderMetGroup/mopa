#' @title Stratified random partitioning into subsets 
#' @description Internal function for stratified random partitioning into subsets of presence/absence
#' 
#' @param k Integer. Number of subsets. Default is 10
#' @param df Data frame with the variable for stratification in the first column 
#' @return  List with k data frames
#' 
#' 
#' @author M. Iturbide 
#' 
#' 
#' @keywords internal
#' 
#' @importFrom sampling strata



kfold <- function(k = 10, df) {
  df <- df
  k <- k  
  sz <- floor(nrow(df)/k)
  strat <- list()
  counter <- 1
  while (sum(counter) < k) {
    #     print(sum(counter))
    st <- strata(df, size = sz, method = "srswor")
    length(strat) <- sum(counter)
    strat[[sum(counter)]] <- df[st[,1], ]
    df <- df[-st[,1],]
    counter <- c(counter, 1)
  }
  length(strat) <- k
  strat[[k]] <- df
  names(strat) <- paste("fold", 1:k, sep="")
  return(strat)
}
