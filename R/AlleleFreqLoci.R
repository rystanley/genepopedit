# Function to calculate the frequency of a given allele
#' @title Calculate allele frequency for a given allele
#' @description sub-function used by genepop_allelefreq
#' @param x Vector of characters representing alleles
#' @param major is the major allele for a given loci
#' @rdname AlleleFreqLoci
#' @export

AlleleFreqLoci <- function(x,major)
{

  #allele coding length
  alleleEx <- max(sapply(x,FUN=function(x){nchar(as.character(x[!is.na(x)]))})) #presumed allele length

  #check to make sure the allele length is a even number
  if(!alleleEx %% 2 ==0){stop(paste("The length of each allele is assumed to be equal (e.g. loci - 001001 with 001 for each allele), but a max loci length of", alleleEx, "was detected. Please check data."))}

    if(alleleEx==2) # four character locus - 2 digit allele code
  {
    x1 <- as.numeric(substring(x,1,1)) # get first allele
    x2 <- as.numeric(substring(x,2,2)) # get second allele
  }

  if(alleleEx==4) # four character locus - 2 digit allele code
  {
    x1 <- as.numeric(substring(x,1,2)) # get first allele
    x2 <- as.numeric(substring(x,3,4)) # get second allele
  }

  if(alleleEx==6) # six character locus - 3 digit allele code
  {
    x1 <- as.numeric(substring(x,1,3)) # get first allele
    x2 <- as.numeric(substring(x,4,6)) # get second allele
  }

  x3 <- c(x1,x2) #combine together in one string
  if(length(unique(x3))>2){x4 <- x3[-which(x3==0)]} # delete the zeros (missing values)
  if(length(unique(x3))<3){x4 <- x3}
  # warning if data is coded wrong (too many alleles or 0 isn't the reference to )
  if(length(unique(x4))>2){print(paste("More than two alleles present. ",
                                       length(unique(x4))," unique alleles in vector:",
                                       paste(unique(x4),collapse=","),sep=""))}

  return(length(which(x4==major))/length(x4)) # calculate frequency of each
  }
