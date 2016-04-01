# Update the sample IDs of a genepop file
#' @title Add seperation ("_") between population name and sample number
#' @description Function to add seperation ("_") between population name and sample number
#' @param GenePop the genepop data to be manipulated. This can be either a file path
#' or a dataframe read in with tab seperation, header=FALSE , quote="", and stringsAsFactors=FALSE.
#' This will be the standard genepop format with the first n+1 rows corresponding to the n loci names,
#' or a single commma delimited row of loci names followed by the locus data. Populations are
#' seperated by "Pop". Each individual ID is linked to the locus data by "   , " (space,space space) and is read in as
#' as a single row (character).
#' @param path the filepath and filename of output.
#' @rdname genepop_ID
#' @importFrom data.table fread as.data.table
#' @export


genepop_ID <- function(GenePop,path){

  #Check to see if GenePop is a data.frame from the workspace and convert to data.table
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
    GenePop <- data.table::as.data.table(GenePop,stringsAsFactors = FALSE)
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

  #Use the 'Pop' original seperators to separate the unique population names
  PopLengths=NULL
  for (i in 1:(length(tempPops)-1)){
    pl=length((tempPops[i]+1):(tempPops[i+1]-1))
    PopLengths=c(PopLengths,pl)}
  #add the length of the last pop which is the only one not bound by a "pop" label
  PopLengths=c(PopLengths,length((tempPops[length(tempPops)]+1):nrow(snpData)))
  popvector=rep(1:length(PopLengths),times=PopLengths) #vector to differentiate the pops based on the locaton of the "Pop" labels.

  snpData <- snpData[-tempPops,] #remove the pop labels.

  #Seperate the snpdata
  temp <- as.data.frame(do.call(rbind, strsplit(snpData$data," ")))
  temp2 <- temp[,4:length(temp)] #split characters by spaces

  #Contingency to see if R read in the top line as the "stacks version"
  if (length(temp2)!=length(ColumnData)){colnames(temp2) <- c(stacks.version,ColumnData)}
  if (length(temp2)==length(ColumnData)){colnames(temp2) <- ColumnData}
  if (length(temp2)!=length(ColumnData)){stacks.version="No STACKS version specified"}

  #stacks version character
  stacks.version <- as.character(stacks.version)

  ## Get the population names (prior to the _ in the Sample ID)
  SampleID <- as.character(temp[,1]) # Sample names of each
  NamePops <- as.character(temp[,1]) # Sample names of each

  #add the "_" to differentiate the populations from the sample number in each sample ID.
  #This asssumes there is a common string among the sample IDs between each "Pop" label in the Genepopfile.
  for(i in unique(popvector)){
    commonname <- common_string(SampleID[which(popvector==i)])
    NamePops[which(popvector==i)] <- paste0(commonname,"_")
    SampleID[which(popvector==i)] <- gsub(commonname,paste0(commonname,"_"),SampleID[which(popvector==i)])
  }

  NameExtract <- substr(NamePops,1,regexpr("_",NamePops)-1)

    #the number of individuals for all popualtions but the last (Pop tagged to the end)
    PopLengths <- table(factor(NamePops, levels=unique(NamePops)))[-length(table(NamePops))]

    if(length(table(NameExtract))==2){PopPosition = PopLengths+1}

    if(length(table(NameExtract))>2){
      PopPosition <- c(PopLengths[1]+1,rep(NA,(length(PopLengths)-1)))
      for (i in 2:length(PopLengths)){
        PopPosition[i] <- PopLengths[i]+PopPosition[i-1]
      }
    }

    #Now stitch the data together
    # paste together the Loci as one long integer seperated for each loci by a space
    Loci <- do.call(paste,c(temp2[,], sep=" "))

    #Paste these to the Loci
    Loci <- paste(SampleID,Loci,sep=" ,  ")

    #Insert the value of "Pop" which partitions the data among populations #only if more than one population
    if(length(table(NamePops))!=1){Loci <- insert_vals(Vec=Loci,breaks=PopPosition,newVal="Pop")}

    #Add the first "Pop" label
    Loci <- c("Pop",Loci)

    Output <- c(as.character(stacks.version),names(temp2),Loci)

    # Save the file
    write.table(Output,path,col.names=FALSE,row.names=FALSE,quote=FALSE)

} #End function
