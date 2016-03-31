# BGC Allele conversion - this function is built specifically for genepop_bgc in genepopedit
#' @title Function to convert allele in character format to bgc zygosity format.
#' @description sub-function used by genepop_bgc
#' @param Vector of characters representing alleles
#' @param allele_length the number of characters in an allele (default: 6 e.g. 130120)
#' @rdname majorminor
#' @export

majorminor<- function(Vec,allele_length=6){
  Vec <- as.character(Vec)
  firstAllele <-  as.data.frame(sapply(Vec,function(x)as.numeric(as.character(substring(x,1,(allele_length/2))))),stringsAsFactors = F)
  secondAllele <-  as.data.frame(sapply(Vec,function(x)as.numeric(as.character(substring(x,(allele_length/2)+1,allele_length)))),stringsAsFactors = F)
  x=c(firstAllele[,1],secondAllele[,1])
  AlleleMajor <- as.numeric(names(which(table(x)==max(table(x)))))
  AlleleMinor <- as.numeric(names(which(table(x)==min(table(x)))))
  if(length(AlleleMajor)>1){AlleleMajor <- AlleleMajor[1];AlleleMinor=AlleleMinor[2]}

  Vec[is.na(Vec)]="-9 -9" #missing data
  Vec=gsub(paste0(AlleleMajor,AlleleMajor),"2 0",Vec) #homozygous major
  Vec=gsub(paste0(AlleleMajor,AlleleMinor),"1 1",Vec) #heterozygous
  Vec=gsub(paste0(AlleleMinor,AlleleMajor),"1 1",Vec) #heterozygous
  Vec=gsub(paste0(AlleleMinor,AlleleMinor),"0 2",Vec) #homozygous minor
  return(Vec)
}
