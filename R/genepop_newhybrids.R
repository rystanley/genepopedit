# Genepop -> New Hybrids
#' @title Convert Genepop to New Hybrids format.
#' @description Function to convert Genepop to NewHybrids
#' @param GenePop the genepop data to be manipulated. This can be either a file path
#' or a dataframe read in with tab seperation, header=FALSE , quote="", and stringsAsFactors=FALSE.
#' This will the standard genepop format with a the first n+1 rows corresponding the the n loci names,
#' or a single commma deliminated row of loci names followed by the locus data. Populations are
#' seperated by "Pop". Each individual ID is linked the the locus data by "  , " and is read in as
#' as single row (character).
#' @param path the filepath and filename of output.
#' @rdname genepop_newhybrids
#' @importFrom tidyr separate
#' @export


genepop_newhybrids <- function(GenePop,path=NULL){


  #Check to see if Genepop is a file path or dataframe
  if(is.character(GenePop)){
    GenePop <- read.table(GenePop,
                          header = FALSE, sep = "\t",
                          quote = "", stringsAsFactors = FALSE)
  }

## check if loci names are read in as one large character vector (1 row)
header <- GenePop[1,]
if(length(gregexpr(',', header, fixed=F)[[1]])>1){
  lociheader <- strsplit(header,",")
  lociheader <- gsub(" ","",unlist(lociheader))
  #remove the first column of loci names
  GenePop <- as.vector(GenePop)
  GenePop <- GenePop[-1,]
  GenePop <- c(lociheader,GenePop)
  GenePop <- data.frame(GenePop,stringsAsFactors = FALSE)
}

## Stacks version information
stacks.version <- GenePop[1,] # this could be blank or any other source. First row is ignored by GenePop

#Remove first label of the stacks version
GenePop <- as.vector(GenePop)
GenePop <- GenePop[-1,]

#Add an index column to Genepop and format as a dataframe
GenePop <- data.frame(data=GenePop,ind=1:length(GenePop))
GenePop$data <- as.character(GenePop$data)

#ID the rows which flag the Populations
Pops  <-  which(GenePop$data == "Pop" | GenePop$data =="pop" | GenePop$data == "POP")
npops  <-  1:length(Pops)

## Seperate the data into the column headers and the rest
ColumnData <- GenePop[1:(Pops[1]-1),"data"]
snpData <- GenePop[Pops[1]:NROW(GenePop),]

#Get a datafile with just the snp data no pops
tempPops <- which(snpData$data=="Pop"| snpData$data =="pop" | snpData$data == "POP") ## Changed because we allowed
## alternate spelling on line 48, so had to change this so it would identify properly and not make an empty DF
snpData <- snpData[-tempPops,]

#Seperate the snpdata
#First we pull out the population data which follows
#"TEXT ,  "
temp <- tidyr::separate(snpData,data,into=c("Pops","snps"),sep=",")
temp$snps <- substring(temp$snps,3) # delete the extra spaces at the beginning
temp2 <- as.data.frame(do.call(rbind, strsplit(temp$snps," "))) #split characters by spaces

#Contingency to see if R read in the top line as the "stacks version"
if (length(temp2)!=length(ColumnData)){colnames(temp2) <- c(stacks.version,ColumnData)}
if (length(temp2)==length(ColumnData)){colnames(temp2) <- ColumnData}
if (length(temp2)!=length(ColumnData)){stacks.version="No STACKS version specified"}

#paste the loci and row individual count together.
temp3=paste(1:nrow(temp),temp$snps,sep=" ")

#compile the output
Output=c(paste0("NumIndivs ",nrow(temp)),
         paste0("NumLoci ",length(temp2)),
         paste0("Digits ",nchar(as.character(temp2[1,1]))/2),
         "Format Lumped",
         paste0("LocusNames ",do.call(paste,c(as.list(colnames(temp2))))),
         temp3)

write.table(Output,path,col.names=FALSE,row.names=FALSE,quote=FALSE)

}