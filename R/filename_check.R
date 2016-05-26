# Get File name specified by GenePop variable
#' @title filename_check
#' @description Gets the file name that is specified by GenePop variable. This function is used by various functions of genepopedit and is applied to the Genepop variable.
#' @param X what to check
#' @rdname filename_check
#' @export
#' @importFrom stringr str_split


filename_check <- function(X){

  GeneSplit <- unlist(stringr::str_split(string = X, pattern = "/"))
  GeneSplit <- GeneSplit[grep(x = GeneSplit, pattern = ".txt")]

  return(GeneSplit)

}
