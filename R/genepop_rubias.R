# Genepop -> Rubias
#' @title Convert Genepop to rubias format.
#' @description Function to convert Genepop to rubias format for individual assignment and mixture analyses.
#' @param genepop the genepop data to be manipulated. This is read in as a complete file path.
#' This will the standard genepop format with a the first n+1 rows corresponding the the n loci names,
#' or a single comma deliminated row of loci names followed by the locus data. Populations are
#' separated by "Pop". Each individual ID is linked to the locus data by " ,  " (space,space space) and is read in as
#' as a single row (character).
#' @param sampletype Specifies whether your rubias input file is a reference (baseline) or mixture file. Can only be one of "reference" or "mixture"
#' @param repunits A vector from the workspace that specifies which individuals are assigned which reporting unit. The length of this vector
#' should be equal to the number of individuals in the original genepop file. Can also use genepop_aggregate to generate this vector.
#' @param path file path to directory where the rubias input file will be saved.
#' @rdname genepop_rubias
#' @importFrom utils write.table
#' @importFrom utils read.table
#' @export

genepop_rubias<-function(genepop,sampletype,repunits,path){

  if(!((sampletype == "reference") || (sampletype == "mixture")))
    stop("sampletype must be one of 'reference' or 'mixture'. Function stopped.",call. = F)

  #flatten a genepop file
  rubiasformat<-genepopedit::genepop_flatten(genepop)

   rubiasformat$repunit<-repunits
  #make everything characters
  for(i in 1:ncol(rubiasformat)){
    rubiasformat[, i] = as.character(rubiasformat[, i])
  }

  #Add the columns and col names rubias wants
  colnames(rubiasformat)[1]<-"indiv"
  colnames(rubiasformat)[2]<-"collection"
  locusnames<-colnames(rubiasformat[,-c(1:3)])
   rubiasformat$sample_type<-rep(sampletype,length(rubiasformat$indiv))

   #Reorganize the columns to preferred rubias order
  rubiasformat1<-rubiasformat[c("sample_type","repunit","collection","indiv",locusnames)]

  #Now split each locus into two columns with one allele in each
  alleleEx <- max(sapply(rubiasformat1[, 5], FUN = function(x) {
    nchar(as.character(x[!is.na(x)]))
  }))

  firstAllele <- as.data.frame(sapply(rubiasformat1[, -c(1:4)], function(x) as.character(substring(x, 1, alleleEx/2))))
  secondAllele <- as.data.frame(sapply(rubiasformat1[, -c(1:4)], function(x) as.character(substring(x, (alleleEx/2) + 1, alleleEx))))

  #Add _1 to the second allele name for each locus
  colnames(secondAllele) <- paste0(colnames(secondAllele), "_1")
  splitloci <- cbind(firstAllele, secondAllele)
  indx <- rbind(names(firstAllele), names(secondAllele))
  splitloci <- splitloci[, indx]
  rubiasinput <- as.data.frame(sapply(splitloci[, -c(1:4)], function(x) as.character(as.factor(x))), stringsAsFactors = FALSE)
  rubiasinput1<-cbind(rubiasformat1[,1:4],rubiasinput)

  #Make missing data NA
  rubiasinput1[rubiasinput1 == "000"] <- NA

  #write the file
  write.table(rubiasinput1,file=paste0(path,"RubiasInput.txt"),quote=F, row.names = F)

}
