# Check that file paths end in "/"
#' @title path_ending
#' @description Checks that the paths end in "/" as required
#' @param path path to check
#' @rdname path_ending
#' @export
#' @importFrom stringi stri_sub

path_ending <- function(path){

  if(stringi::stri_sub(path, -1) != "/"){
    path <- paste0(path, "/")
          }
  if(stringi::stri_sub(path, -1) == "/"){path <- path}

  return(path)

}
