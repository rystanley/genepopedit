# Genepop -> FSTAT
#' @title Convert Genepop to FSTAT format.
#' @description Function to convert Genepop to FSTAT
#' @param GenePop the genepop data to be manipulated. This can be either a file path
#' or a dataframe read in with tab separation, header=FALSE , quote="", and stringsAsFactors=FALSE.
#' This will be the standard genepop format with the first n+1 rows corresponding to the n loci names,
#' or a single comma delimited row of loci names followed by the locus data. Populations are
#' separated by "Pop". Each individual ID is linked to the locus data by " ,  " (space,space space) and is read in as
#' as a single row (character).
#' @param path the filepath and filename of output.
#' @param addworkspace logical statement defining whether the converted data
#' should be saved as a file specified in the path (default) argument or whether it should be returned to the workspace
#' if returned to the workspace the object will be called "Output_fstat".
#' @rdname genepop_fstat
#' @importFrom data.table fread
#' @importFrom utils write.table
#' @export


##
genepop_fstat <- function(GenePop,path=NULL,addworkspace=FALSE)
  {

  #Check to see if GenePop is a data.frame from the workspace and convert to data.table
  if(is.data.frame(GenePop)){GenePop <- as.data.table(GenePop)}

  #Check to see if Genepop is a file path or dataframe
  if(is.character(GenePop)){
    GenePop <- data.table::fread(GenePop,
                                 header = FALSE, sep = "\t",
                                 stringsAsFactors = FALSE)
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
    GenePop <- as.data.table(GenePop,stringsAsFactors = FALSE)
  }

  ## Stacks version information
  stacks.version <- GenePop[1,] #this could be blank or any other source. First row is ignored by GenePop

  #Remove first label of the stacks version
  GenePop <- GenePop[-1,]
  colnames(GenePop) <- "data"

  #ID the rows which flag the Populations
  Pops  <-  which(GenePop$data == "Pop" | GenePop$data =="pop" | GenePop$data == "POP")
  npops  <-  1:length(Pops)

  ## separate the data into the column headers and the rest
  ColumnData <- GenePop$data[1:(Pops[1]-1)]
  ColumnData <- gsub("\r","",ColumnData)#remove any hidden carriage returns
  snpData <- GenePop[Pops[1]:NROW(GenePop),]

  #Get a datafile with just the snp data no pops
  tempPops <- which(snpData$data=="Pop"| snpData$data =="pop" | snpData$data == "POP") ## Changed because we allowed
  snpData <- snpData[-tempPops,]

  #separate the snpdata
  temp <- as.data.frame(do.call(rbind, strsplit(snpData$data," ")))

  #data format check
  if(unique(temp[,2])!="," | !length(which(temp[,3]==""))>1){
    stop("Genepop sampleID delimiter not in proper format. Ensure sampleIDs are separated from loci by ' ,  ' (space comma space space). Function stopped.",call. = FALSE)
  }

  temp2 <- temp[,4:length(temp)] #split characters by spaces

  #Contingency to see if R read in the top line as the "stacks version"
  if (length(temp2)!=length(ColumnData)){colnames(temp2) <- c(stacks.version,ColumnData)}
  if (length(temp2)==length(ColumnData)){colnames(temp2) <- ColumnData}
  if (length(temp2)!=length(ColumnData)){stacks.version="No STACKS version specified"}

  #stacks version character
  stacks.version <- as.character(stacks.version)

  ## Get the population names (prior to the _ in the Sample ID)
  NamePops <- temp[,1] # Sample names of each
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
   Loci <- paste(numPop,Loci,sep="   ") #separated by three spaces

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
      Output_fstat <- as.data.frame(Output)
    }
    # Save the file
    if(!addworkspace){
      utils::write.table(Output,path,col.names=FALSE,row.names=FALSE,quote=FALSE)
    }
} #End function
