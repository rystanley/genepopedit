# Function for inserting rows
#' @title Insert values at location
#' @description Function to insert values at specified positons within a dataframe and-or vector
#' @param Vec A data vector where values will be inserted.
#' @param breaks vector of positions where values will be inserted.
#' Here the output will have length(vec)+length(breaks) rows with the values
#' specified in 'newVal' positioned at the 'breaks' positions
#' @param newVal values to be inserted
#' @rdname insert_vals
#' @export

##
insert_vals <- function(Vec,breaks,newVal){
  break.space <- 1:(length(breaks))
  breaks <- breaks+break.space-1 #To space out the insertion points.
  newvec <- rep(NA,length(Vec)+length(breaks)) #Pre-allocate memory by creating final dataframe.
  for(i in 1:length(breaks)){newvec[breaks[i]]=newVal} #Insert added rows into new dataframe>
  x <- 1:length(newvec)
  x <- x[-(breaks)] #Finding the rows of the new dataframe that will receive old rows
  for(i in 1:length(Vec)){newvec[x[i]]=Vec[i]}
  return(newvec)}
