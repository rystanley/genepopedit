# Decision matrix for selecting a panel of unlinked loci by Fst
#' @title Decision matrix for selecting a panel of unlinked loci by Fst
#' @description Decision matrix for selecting unlinked loci while maximizing Fst. Function is used by genepop_toploci.
#' @param Matrix ordered by Loci high (1) to low (n loci) by fst in rows and the number of columns corresponds to the linked loci.
#' @rdname Optimfunc
#' @export

Optimfunc <- function(x)
{
  subdat <- as.matrix(x)
  returned <- NULL
  nloci <- length(subdat)
  #while(nrow(subdat)>0)
  CheckProgress <- txtProgressBar(min = 0, max = nrow(subdat), style = 3)
  while(class(subdat)!="integer")
  {
    highest <- subdat[1,which.min(subdat[1,])] #loci returned to the panel
    otherdat <- subdat[1,-which.min(subdat[1,])] #loci which will be eliminated from the panel due to being linked to the top loci by fst
    otherdat <- otherdat[!is.na(otherdat)]
    if(length(!is.na(otherdat))>0)
      {
    for(i in otherdat){subdat[subdat==i] <- NA}
    subdat <- subdat[apply(subdat,1,function(x){!highest%in%x}),]
    }

    returned <- c(returned,highest)
    if(class(subdat)!="integer" & sum(!is.na(subdat[1,]))<=1){subdat <- subdat[-1,]}
    if(is.null(nrow(subdat))){returned <- c(returned,subdat[1])}

    # print(nrow(subdat))
    setTxtProgressBar(CheckProgress, nrow(subdat))
  }
  return(returned)
}
