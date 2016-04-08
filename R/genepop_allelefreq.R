# GenePop allelefreq
#' @title Explore population specific allele frequencies.
#' @description Function returns population derived allele frequencies.
#' @param GenePop the genepop data to be manipulated. This can be either a file path
#' or a dataframe read in with tab seperation, header=FALSE , quote="", and stringsAsFactors=FALSE.
#' This will be the standard genepop format with the first n+1 rows corresponding to the n loci names,
#' or a single commma delimited row of loci names followed by the locus data. Populations are
#' seperated by "Pop". Each individual ID is linked to the locus data by " ,  " (space space,space) and is read in as
#' as a single row (character).
#' @param popgroup if specified (Default: NULL) popgroup is a dataframe or path to a csv.
#' This dataframe contains two columns. Column 1 corresponds to the population names. These names
#' should match the individual IDs (e.g. BON_01 ,  110110 would be 'BON'). The next column
#' has the group. If groupings are the same as populations then leave as NULL (Default). Populations will be defined in this function using the population grouping format in the GENEPOP file.
#' @rdname genepop_detective
#' @importFrom data.table fread as.data.table
#' @export
#'

genepop_allelefreq <- function(GenePop,popgroup=NULL){

  #Check to see if GenePop is a data.frame from the workspace
  if(is.data.frame(GenePop)){GenePop <- data.table::as.data.table(GenePop)}

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
    GenePop <- data.frame(GenePop,stringsAsFactors = FALSE)
  }

  ## Stacks version information
  stacks.version <- GenePop[1,] #this could be blank or any other source. First row is ignored by GenePop

  #Remove first label of the stacks version
  GenePop <- GenePop[-1,]
  colnames(GenePop) <- "data"

  #ID the rows which flag the Populations
  Pops  <-  which(GenePop$data == "Pop" | GenePop$data =="pop" | GenePop$data == "POP")
  npops  <-  1:length(Pops)

  ## Seperate the data into the column headers and the rest
  ColumnData <- GenePop$data[1:(Pops[1]-1)]
  ColumnData <- gsub("\r","",ColumnData)#remove any hidden carriage returns
  snpData <- GenePop[Pops[1]:NROW(GenePop),]

  #Get a datafile with just the snp data no pops
  tempPops <- which(snpData$data=="Pop"| snpData$data =="pop" | snpData$data == "POP") ## Changed because we allowed
  snpData <- snpData[-tempPops,]

  #Seperate the snpdata
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

  ## Get the population names (prior to the _ in the Sample ID)
  NamePops <- temp[,1] # Sample names of each
  NameExtract <- substr(NamePops,1,regexpr("_",NamePops)-1)

    PopNum <- data.frame(table(NameExtract))
    colnames(PopNum)[1] <- "Population"

    #convert the snp data into character format to get rid of factor levels

    #get the allele values summary header
    firstAllele <-  as.data.frame(sapply(temp2,function(x)as.numeric(as.character(substring(x,1,(nchar(alleleEx)/2))))))
    secondAllele <-  as.data.frame(sapply(temp2,function(x)as.numeric(as.character(substring(x,(nchar(alleleEx)/2)+1,nchar(alleleEx))))))

    ## population grouping variables
    pPops <- NULL
    for (i in 2:length(tempPops)){
      pPops <- c(pPops,tempPops[i]-(i-1))
    }

    pPops <- c(1,pPops)

    #count populatons
    lPops <- NULL
    for(i in 1:(length(pPops)-1)){
      if(i==1){lPops=pPops[i+1]-1}
      if(i < (length(pPops))){lPops <- c(lPops,pPops[i+1]-pPops[i])}else{
        lPops=c(lPops,nrow(temp2)-(pPops[i]-1))
      }
    }

    PopGroupVec <- rep(1:length(lPops),times=lPops)
    PopGroupVec <- PopGroupVec[order(PopGroupVec)]

    #Allele Frequency table
    AlleleDat <- rbind(firstAllele,secondAllele)
    AlleleDat[AlleleDat==0]=NA #replace 0's with NA
    AlleleDat$Pops <- rep(PopGroupVec,times=2) # add population tag

    temp2$Pops <- PopGroupVec

    ## Identify major alleles for each loci
    mallele <- temp2[,-length(temp2)]%>%summarise_each(funs(AlleleFreq(.)))
    majorfreqs <- as.vector(mallele[1,]);colnames(majorfreqs)=NULL
    mallele <- data.frame(loci=names(mallele),major=mallele[1,],stringsAsFactors = F)

    transposeData <- reshape2::melt(temp2,id.vars="Pops")

}