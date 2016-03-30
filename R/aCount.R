# BGC Allele conversion - this function is built specifically for genepop_bgc in genepopedit
#' @title Function used when converting allele formats for Bayesian Genomic Clines admixed file.
#' @description BGC admixed allele conversion
#' @param x flattened allele dataframe
#' @rdname aCount
#' @export


aCount <- function(x){
  #homozygous minor
  if(length(unique(as.numeric(x[c("allele1","allele2")])))==1 &
     sum(as.numeric(x[c("allele1","allele2")])!= -9)>1 &
     unique(as.numeric(x["allele1"])==as.numeric(x["alleleMinor"]))){return(c(0,2))} # same alleles (small)
  #homozygous major
  if(length(unique(as.numeric(x[c("allele1","allele2")])))==1 &
     sum(as.numeric(x[c("allele1","allele2")])!= -9)>1 &
     unique(as.numeric(x["allele1"])==as.numeric(x["alleleMajor"]))){return(c(2,0))} # same alleles (small)
  #heterozygous
  if(length(unique(as.numeric(x[c("allele1","allele2")])))==2 &
     sum(as.numeric(x[c("allele1","allele2")])!= -9)>1) {return(c(1,1))}
  #Missing data
  if(length(which(as.numeric(x[c("allele1","allele2")])== -9))>1){return(c(-9,-9))}
  #Missing data -- contingency for all missing data for a given loci and population
  if(length(which(as.numeric(x[c("allele1","allele2")])== 999))>1){return(c(-9,-9))}
}
