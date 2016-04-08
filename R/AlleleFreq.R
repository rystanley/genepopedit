# Function to calculate allele frequencies
#' @title Calculate allele frequencies
#' @description sub-function used by genepop_allelefreq
#' @param x Vector of characters representing alleles
#' @param Freq vector defining which allele frequency is returned (Major, Minor or Missing)
#' @rdname AlleleFreq
#' @export

AlleleFreq <- function(x,Freq="Major")
{
  # Freq vector of 2 ("Major" or "Minor" or "missing") allele frequencie or percent missing

  xtest <- x[1:10]
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
  if(length(unique(x3))>2){x4 <- x3[-which(x3==0)]} # delete the zeros (missing values)
  if(length(unique(x3))<3){x4 <- x3}
  # warning if data is coded wrong (too many alleles or 0 isn't the reference to )
  if(length(unique(x4))>2){print(paste("More than two alleles present. ",
                                       length(unique(x4))," unique alleles in vector:",
                                       paste(unique(x4),collapse=","),sep=""))}

  x5 <- as.data.frame(table(x4)/sum(length(x4))) # calculate frequency of each
  colnames(x5)=c("allele","Freq")

  #For missing values
  x6 <- as.data.frame(table(x3)/sum(length(x3))) # calculate frequency of each

  # if it is 50 - 50 return the first allele
  if(x5[1,2]==0.5){x5[1,2]=x5[1,2]+0.01} # add 0.01% to break the tie

  if(Freq=="Major"){return(x5[which(x5$Freq==max(x5$Freq)),])}
  if(Freq=="Minor"){return(x5[which(x5$Freq==min(x5$Freq)),])}
  if(Freq=="Missing"){return(x6[which(x6[,1]==0),])}

}
