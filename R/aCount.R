# BGC Allele conversion
#' @title Function used when converting allele formats for Bayesian Genomic Clines admixed file.
#' @description BGC admixed allele conversion
#' @param x flattend allele dataframe
#' @rdname aCount
#' @export


aCount <- function(x){
  if(length(unique(as.numeric(x[c("allele1","allele2")])))==1 &
     sum(as.numeric(x[c("allele1","allele2")])!= -9)>1 &
     unique(as.numeric(x["allele1"])==as.numeric(x["alleleMin"]))){return(c(1,1))} # same alleles (small)

  if(length(unique(as.numeric(x[c("allele1","allele2")])))==1 &
     sum(as.numeric(x[c("allele1","allele2")])!= -9)>1 &
     unique(as.numeric(x["allele1"])==as.numeric(x["alleleMax"]))){return(c(2,2))} # same alleles (small)

  if(as.numeric(x[,"allele1"])>as.numeric(x[,"allele2"])){return(c(2,0))}
  if(as.numeric(x[,"allele1"])<as.numeric(x[,"allele2"])){return(c(0,2))}
  if(unique(as.numeric(x[c("allele1","allele2")]))== -9){return(c(0,0))}
}
