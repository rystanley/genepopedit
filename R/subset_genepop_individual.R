# Subset samples from Genepop file
#' @title Genepop remove or keep specific sample IDs
#' @description Function for the manipulation of genopop format SNP datasets
#' @param genepop the genepop data to be manipulated. This can be either a file path
#' or a dataframe read in with tab separation, header=FALSE , quote="", and stringsAsFactors=FALSE.
#' This will be the standard genepop format with the first n+1 rows corresponding to the n loci names,
#' or a single comma delimited row of loci names followed by the locus data. Populations are
#' separated by "Pop". Each individual ID is linked to the locus data by " ,  " (space,space space) and is read in as
#' as a single row (character).
#' @param indiv vector sample IDs of interest.
#' These can be either the order by which they occur or the exact name of the loci
#' indiv <- \code{c("Pop01_01","Pop03_15","Pop16_02")} would individuals with these sample names.
#' @param keep logical whether to delete sample IDs specified by indiv (default: TRUE) or delete all other IDs.
#' @param path the filepath and filename of output.
#' @rdname subset_genepop_individual
#' @importFrom data.table fread as.data.table
#' @importFrom utils write.table
#' @export

##
subset_genepop_individual <- function(genepop,indiv=NULL,keep=FALSE,path){

  #Check to see if genepop is a data.frame from the workspace and convert to data.table
  if(is.data.frame(genepop)){genepop <- data.table::as.data.table(genepop)}

  #Check to see if genepop is a file path or dataframe
  if(is.character(genepop)){
    genepop <- data.table::fread(genepop,
                                 header = FALSE, sep = "\t",
                                 stringsAsFactors = FALSE)
  }

  ## check if loci names are read in as one large character vector (1 row)
  header <- genepop[1,]
  if(length(gregexpr(',', header, fixed=F)[[1]])>1){
    lociheader <- strsplit(header,",")
    lociheader <- gsub(" ","",unlist(lociheader))
    #remove the first column of loci names
    genepop <- as.vector(genepop)
    genepop <- genepop[-1,]
    genepop <- c(lociheader,genepop)
    genepop <- data.table::as.data.table(genepop,stringsAsFactors = FALSE)
  }

  ## Stacks version information
  stacks.version <- genepop[1,] #this could be blank or any other source. First row is ignored by genepop

  #Remove first label of the stacks version
  genepop <- genepop[-1,]
  colnames(genepop) <- "data"

  #ID the rows which flag the Populations
  Pops  <-  which(genepop$data == "Pop" | genepop$data =="pop" | genepop$data == "POP")
  npops  <-  1:length(Pops)

  ## separate the data into the column headers and the rest
  ColumnData <- genepop$data[1:(Pops[1]-1)]
  ColumnData <- gsub("\r","",ColumnData)#remove any hidden carriage returns
  snpData <- genepop[Pops[1]:NROW(genepop),]

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
    NamePops <- as.character(temp[,1]) # Sample names of each
    NameExtract <- substr(NamePops,1,regexpr("_",NamePops)-1)

## now subset the individuals

    if(keep){
        temp <- temp[which(NamePops %in% indiv),]
        temp2 <- temp2[which(NamePops %in% indiv),]
        NameExtract2 <- NameExtract[which(NamePops %in% indiv)]
    }

    if(!keep){
      temp <- temp[which(NamePops %in% setdiff(NamePops,indiv)),]
      temp2 <- temp2[which(NamePops %in% setdiff(NamePops,indiv)),]
      NameExtract2 <- NameExtract[which(NamePops %in% setdiff(NamePops,indiv))]
    }

    #Now recompile the GenePop format

    #Create a new vector with the new population names
    NameExtract3 <- NameExtract2
    for (i in 1:length(unique(NameExtract2))){
      NameExtract3[which(NameExtract3==unique(NameExtract3)[i])]=i
    }

    #the number of individuals for all populations but the last (Pop tagged to the end)
    PopLengths <- table(factor(NameExtract3, levels=unique(NameExtract3)))[-length(table(NameExtract3))]

    if(length(table(NameExtract3))==2){PopPosition = PopLengths+1}

    if(length(table(NameExtract3))>2){
      PopPosition <- c(PopLengths[1]+1,rep(NA,(length(PopLengths)-1)))
      for (i in 2:length(PopLengths)){
        PopPosition[i] <- PopLengths[i]+PopPosition[i-1]
      }
    }

    #Now stitch the data together
    # paste together the Loci as one long integer separated for each loci by a space
    Loci <- do.call(paste,c(temp2[,], sep=" "))

    PopVec <- paste0(as.character(temp[,1])," ,  ")

    #Paste these to the Loci
    Loci <- paste(PopVec,Loci,sep="")

    #Insert the value of "Pop" which partitions the data among populations #only if more than one population
    if(length(table(NameExtract2))!=1){Loci <- insert_vals(Vec=Loci,breaks=PopPosition,newVal="Pop")}

    #Add the first "Pop" label
    Loci <- c("Pop",Loci)

    #Derive the genepop output
    Output <- c(stacks.version,names(temp2),Loci)

    # Save the file
    utils::write.table(Output,path,col.names=FALSE,row.names=FALSE,quote=FALSE)

} #End function
