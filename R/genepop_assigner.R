# Genepop -> assigner
#' @title Convert Genepop to assigner format.
#' @description Function to convert Genepop to assigner
#' @param GenePop the genepop data to be manipulated. This can be either a file path
#' or a dataframe read in with tab separation, header=FALSE , quote="", and stringsAsFactors=FALSE.
#' This will be the standard genepop format with a the first n+1 rows corresponding to the n loci names,
#' or a single comma delimited row of loci names followed by the locus data. Populations are
#' separated by "Pop". Each individual ID is linked to the locus data by " ,  " (space,space space) and is read in as
#' as a single row (character).
#' @param popgroup if specified (Default: NULL) popgroup is a dataframe or path to a csv.
#' This dataframe contains two columns. Column 1 corresponds to the population names. These names
#' should match the individual IDs (e.g. BON_01 ,  110110 would be 'BON'). The next column
#' has the group. These values must be numeric. If groupings are the same as populations then leave as NULL (Default).
#' @param path the filepath and filename of output.
#' @rdname genepop_assigner
#' @importFrom data.table fread
#' @importFrom utils write.table
#' @importFrom utils read.csv
#' @export


genepop_assigner <- function(GenePop,popgroup=NULL,path=NULL){

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

##Ouput

if(!is.null(popgroup)) #if popgroup isn't NULL
{
  if(is.character(popgroup)){popgroup <- utils::read.csv(popgroup,header=T)} #if it is a path then read it in

  if(length(intersect(unique(NameExtract),popgroup[,1]))!=length(unique(NameExtract))){
        message("Popuation levels missing form popgroups input. asssignr groups now set to default population levels")
        groupvec <- NameExtract
        for (i in 1:length(unique(NameExtract))) # replace with numbers
        {
          groupvec[which(groupvec==unique(NameExtract)[i])] = i
        }
    }

  groupvec=NameExtract
  for (i in 1:nrow(popgroup))
  {
    groupvec[which(groupvec==popgroup[i,1])]=rep(popgroup[i,2],length(groupvec[which(groupvec==popgroup[i,1])]))
  }

}

if(is.null(popgroup)) #if popgroup isn NULL
{
  groupvec <- NameExtract
      for (i in 1:length(unique(NameExtract))) # replace with numbers
        {
         groupvec[which(groupvec==unique(NameExtract)[i])] = i
      }

}

Output <- cbind(groupvec,NamePops,temp2) #dataframe with individual groups, IDs,  & Loci
Output <- apply(Output,2,as.character)
Output <- data.frame(Output,stringsAsFactors = FALSE)
colnames(Output)=c("POP_ID","INDIVIDUALS",colnames(temp2)) #add the column headers
Output$POP_ID <- as.numeric(Output$POP_ID)
Output <- Output[order(Output$POP_ID,decreasing = FALSE),]
utils::write.table(Output,path,row.names=FALSE,col.names = TRUE,quote=FALSE,sep="\t")

}
