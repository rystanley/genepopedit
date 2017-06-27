# Function to calculate allele frequencies from genotype data.
#' @title Calculate allele frequencies - sub-function of genepop_allelefreq()
#' @description sub-function used by genepop_allelefreq
#' @param x Vector of characters representing alleles
#' @rdname frequencies
#' @export

frequencies <- function(x)
{
  x=x$Value
  xtest <- x[1:length(x)]
  AlleleLength <- max(nchar(as.character(xtest)),na.rm=T)

  if(AlleleLength==2) # four character locus - 2 digit allele code
  {
    x1 <- as.numeric(substring(x,1,1)) # get first allele
    x2 <- as.numeric(substring(x,2,2)) # get second allele
  }

  if(AlleleLength==4) # four character locus - 2 digit allele code
  {
    x1 <- as.numeric(substring(x,1,2)) # get first allele
    x2 <- as.numeric(substring(x,3,4)) # get second allele
  }

  if(AlleleLength==6) # six character locus - 3 digit allele code
  {
    x1 <- as.numeric(substring(x,1,3)) # get first allele
    x2 <- as.numeric(substring(x,4,6)) # get second allele
  }

  x3 <- c(x1,x2) #combine together in one string
  #If all zeros return NA
  if(sum(x3==0)==length(x3)){return(data.frame(allele=NA,Freq=NA,stringsAsFactors = F))}

  #If not all zeros return allele frequencies for each locus
  if(!sum(x3==0)==length(x3)){
      x4 <- x3[!x3==0] # delete the zeros (missing values)
      x5 <- as.data.frame(table(x4)/sum(length(x4)),stringsAsFactors = F) # calculate frequency of each
      colnames(x5)=c("allele","Freq")

      return(x5[order(x5$Freq,decreasing = T),])
  }

}

