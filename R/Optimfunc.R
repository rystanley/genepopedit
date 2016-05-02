# Decision matrix for selecting a panel of unlinked loci by fst
#' @title Decision matrix for selecting a panel of unlinked loci by fst
#' @description Decision matrix for selecting unlinked loci while maximizing fst
#' @param Matrix ordered by Loci high (1) to low (n loci) by fst in rows and the number of columns corresponds to the linked loci.
#' @rdname Optimfunc
#' @export

Optimfunc <- function(x)
{
  subdat <- as.matrix(x)
  returned <- NULL
  while(nrow(subdat)>0)
  {
    highest <- subdat[1,which.min(subdat[1,])] #loci returned to the panel
    otherdat <- subdat[1,-which.min(subdat[1,])] #loci which will be eliminated from the panel due to being linked to the top loci by fst
    for(i in otherdat){subdat[subdat==i] <- NA}
    subdat <- subdat[apply(subdat,1,function(x){!highest%in%x}),]
    returned <- c(returned,highest)
  }
  return(returned)
}
