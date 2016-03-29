# Function for identifying common string sequence
#' @title Identify common string sequence
#' @description A function which returns the largest common sequence among a vector of strings.
#' @param Vec A character vector.
#' @rdname common_string
#' @export

##
common_string <- function(Vec){

  Vec.split <- strsplit(Vec, '')
  Vec.split <- lapply(Vec.split, `length<-`, max(nchar(Vec)))
  Vec.mat <- do.call(rbind, Vec.split)
  Vec.common.length <- which.max(apply(Vec.mat, 2, function(col) !length(unique(col)) == 1)) - 1
  return(substr(Vec[1], 1, Vec.common.length))

}#end function
