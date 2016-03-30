# Genepop -> FSTAT
#' @title Convert Genepop to FSTAT format.
#' @description Function to convert Genepop to FSTAT
#' @param GenePop the genepop data to be manipulated. This can be either a file path
#' or a dataframe read in with tab seperation, header=FALSE , quote="", and stringsAsFactors=FALSE.
#' This will be the standard genepop format with the first n+1 rows corresponding to the n loci names,
#' or a single commma deliminated row of loci names followed by the locus data. Populations are
#' seperated by "Pop". Each individual ID is linked to the locus data by " ,  " (space,space space) and is read in as
#' as a single row (character).
#' @param path the filepath and filename of output.
#' @param addworkspace logical statement defining whether the converted data
#' should be saved as a file specified in the path (default) arguement or whether it should be returned to the workspace
#' if returned to the workspace the object will be called "Output_fstat".
#' @rdname genepop_fstat
#' @importFrom tidyr separate
#' @export


##
genepop_fstat <- function(GenePop,path=NULL,addworkspace=FALSE)
  {

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
    ColumnData <- gsub("\r","",ColumnData)#remove any hidden carriage returns
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

#convert the snp data into character format to get rid of factor levels
    temp2[] <- lapply(temp2, as.character)

#get the allele values for FSTAT summary header
    firstAllele <-  as.data.frame(sapply(temp2,function(x)as.numeric(as.character(substring(x,1,3)))))
    secondAllele <-  as.data.frame(sapply(temp2,function(x)as.numeric(as.character(substring(x,4,6)))))

## replace Genepop blanks ("000000") with fstat blanks ("0")
    temp2[temp2=="000000"]="     0"

## population column
    numPop=rep(1:length(table(NameExtract)),
               times=table(factor(NameExtract,levels=unique(NameExtract))))

## Add loci back together as one character vector
   Loci <- do.call(paste,c(temp2[,], sep=" "))

## Add numeric population grouping and spacing
   Loci <- paste(numPop,Loci,sep="   ") #seperated by three spaces

## Add the loci names
   Loci <- c(names(temp2),Loci)

## Add the header information
    nPops <- length(unique(numPop))
    nLoci <- length(temp2)
    allelemax <- max(c(max(as.matrix(firstAllele)),max(as.matrix(secondAllele))))
    allelecallsize <- nchar(temp2[1,1])/2

    FStat_header <- paste(nPops,nLoci,allelemax,allelecallsize,sep=" ")

#compile into fstat format
    Output <- c(FStat_header,Loci)

    #return to the workspace
    if(addworkspace){
      Output_fstat<<-as.data.frame(Output)
    }
    # Save the file
    if(!addworkspace){
      write.table(Output,path,col.names=FALSE,row.names=FALSE,quote=FALSE)
    }
} #End function
