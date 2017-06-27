# Function to calculate allele frequencies
#' @title Calculate the major allele
#' @description sub-function used by genepop_allelefreq
#' @param x Vector of characters representing alleles
#' @rdname AlleleFreq
#' @export

AlleleFreq <- function(x)
{

  xvec <- x[1:(floor(length(x)/2))] #how long are the characters. For speed this only completes the task using half the vector.
  AlleleLength <- max(nchar(as.character(xvec)),na.rm=T)

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
  if(length(unique(x3))>2){x4 <- x3[!x3==0]} # delete the zeros (missing values)
  if(length(unique(x3))<3){x4 <- x3}

 #  warning if data is coded wrong (too many alleles or 0 isn't the reference to )
   if(length(unique(x4))>2){print(paste("More than two alleles present. ",
                                        length(unique(x4))," unique alleles in vector:",
                                        paste(unique(x4),collapse=","),sep=""))}

  x5 <- as.data.frame(table(x4)/sum(length(x4))) # calculate frequency of each
  colnames(x5)=c("allele","Freq")

  #For missing values
  x6 <- as.data.frame(table(x3)/sum(length(x3))) # calculate frequency of each

  # if it is 50 - 50 return the first allele
  if(x5[1,2]==0.5){x5[1,2]=x5[1,2]+0.01} # add 0.01% to break the tie

  return(as.character(x5[which(x5$Freq==max(x5$Freq)),"allele"]))
}
