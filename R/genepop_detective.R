# GenePop detective
#' @title Explore Genepop data structure.
#' @description Funtion returns Genepop file meta-data.
#' @param GenePop the genepop file to be manipulated. This will the standard
#' genepop format with a the first n+1 rows corresponding the the n loci names
#' followed by the locus data. Populations are seperated by "Pop".
#' Each individual ID is linked the the locus data by "  , " and is read in as
#' as single row (character)
#' e.g.
#' Stacks Ver 1.0
#' 1
#' 2
#' 3
#' Pop
#' Pop01_01  , 120120 110110 110110
#' Pop01_02  , 100100 110110 110110
#' Pop
#' Pop02_01  , 120120 110110 110110
#' ...
#' @param variable data to be returned
#' Four options \code{default = "Pops"}
#' "Pops" = vector of population names
#' "PopNum" = dataframe of population names and counts
#' "Inds" = vector of sample IDs
#' "Loci" = vector of Loci
#' @rdname genepop_detective
#' @importFrom tidyr separate
#' @export
#'

genepop_detective <- function(GenePop,variable="Pops"){

## Stacks version information
stacks.version <- GenePop[1,] #this could be blank or any other source. First row is ignored by GenePop

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

## Get the population names (prior to the _ in the Sample ID)
NamePops <- temp[,1] # Sample names of each
NamePops <- gsub(" ","",NamePops) #get rid of space
NameExtract <- substr(NamePops,1,regexpr("_",NamePops)-1)

PopNum <- data.frame(table(NameExtract))
colnames(PopNum)[1] <- "Population"

#return data vector of interest
if(variable=="Pops"){return(unique(NameExtract))}
if(variable=="PopNum"){return(PopNum)}
if(variable=="Inds"){return(NamePops)}
if(variable=="Loci"){return(names(temp2))}

}
