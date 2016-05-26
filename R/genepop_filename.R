# Get File name specified by GenePop variable
#' @title genepop_filename
#' @description Gets the file name that is specified by GenePop variable
#' @param X what to check
#' @rdname genepop_filename
#' @export
#' @importFrom stringr str_split


genepop_filename <- function(X){

  GeneSplit <- unlist(stringr::str_split(string = X, pattern = "/"))
  GeneSplit <- GeneSplit[grep(x = GeneSplit, pattern = ".txt")]

  return(GeneSplit)

}
