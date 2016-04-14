# Subset Genepop Aggregate
#' @title Genepop subset, combine, and reorder populations
#' @description Function to cluster populations together and remove specific loci
#' @param GenePop the genepop data to be manipulated. This can be either a file path
#' or a dataframe read in with tab separation, header=FALSE , quote="", and stringsAsFactors=FALSE.
#' This will be the standard genepop format with the first n+1 rows corresponding to the n loci names,
#' or a single comma delimited row of loci names followed by the locus data. Populations are
#' separated by "Pop". Each individual ID is linked to the locus data by " ,  " (space,space space) and is read in as
#' as a single row (character).
#' @param subs the loci names of interest or a vector
#' subs <- c("190-56","145_21",456_12") would return loci with these defined names.
#' @param keep logical vector which defines whether you want to remove the loci or keep them.
#' The default is to keep them keep <- TRUE assuming you are removing neutral markers
#' and only keeping the subs
#' @param agPopFrame a dataframe or path to a csv.
#' This dataframe contains two columns: Column 1 corresponds to the population names.
#' Here we consider the alpha-numeric characters before the first underscore '_' to be the population name.
#' so that IDs are "Population_sample#" (e.g. Aqua23_04 = Population Aqua23, individual 4).
#' These names can be obtained using the genepop_detective function.
#' The next column has grouping variables. If you don't want to change the grouping just repeat original name.
#' If the input is a dataframe object from the workspace it must be a data.frame object and therefore will have headers.
#' e.g. data.frame(Opop=c("AAA","BBB","CCC","DDD","EEE","FFF","GGG"),AgPop=c("Pop1","Pop2","CCC","Pop2","EEE","Pop2","GGG"))
#' In this case AAA/BBB & DDD/FFF would be clustered together between population flags in genepop.
#' @param path to the output .txt file.
#' @rdname subset_genepop_aggregate
#' @importFrom data.table fread as.data.table
#' @export

subset_genepop_aggregate <- function(GenePop,subs=NULL,keep=TRUE,agPopFrame,path){

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
## Now subset out the the data according to the specified loci and whether or not you want to keep them.

      if(!keep)# neutral
          {
            if(length(subs)>0){reqCols <- temp2[,-which(names(temp2)%in%subs)]}
            if(length(subs)==0){reqCols <- temp2}
          }

        if(keep)# outliers or loci under divergent selection
            {
            if(length(subs)>0){reqCols <- temp2[,c(subs)]}
            if(length(subs)==0){reqCols <- temp2}
            }


## Now subset the Populations

    if(!is.data.frame(agPopFrame)) #if it isn't a dataframe then read in the path
    {
      agPopFrame <- read.csv(agPopFrame,header=T)
    }

    sPop <- as.character(agPopFrame[,1]) # these are the populations of interest

    # is a population subset required
        reqCols <- reqCols[which(NameExtract %in% sPop),]
        temp <- temp[which(NameExtract %in% sPop),]
        temp2 <- temp2[which(NameExtract %in% sPop),]
        NamePops <- NamePops[which(NameExtract %in% sPop)]

    #Now recompile the GenePop format
    NameExtract2 <- NameExtract[which(NameExtract %in% sPop)]
    #Create a new vector with the new population names
    for (i in sPop){
      NameExtract2[which(NameExtract2 == i)]=as.character(agPopFrame[which(agPopFrame[,1]==i),2])
    }

    #reorder variables according to the new population groupings
    reqCols <- reqCols[order(NameExtract2),]
    temp <- temp[order(NameExtract2),]
    temp2 <- temp2[order(NameExtract2),]
    NamePops <- NamePops[order(NameExtract2)]

    NameExtract3 <- NameExtract2[order(NameExtract2)]
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
    Loci <- do.call(paste,c(reqCols[,], sep=" "))

    #Paste these to the Loci
    Loci <- paste(NamePops," ,  ",Loci,sep="")

    #Insert the value of "Pop" which partitions the data among populations #only if more than one population
    if(length(table(NameExtract2))!=1){Loci <- insert_vals(Vec=Loci,breaks=PopPosition,newVal="Pop")}

    #Add the first "Pop" label
    Loci <- c("Pop",Loci)

      if(!keep)
      {
        if(length(subs)==0){Output <- c(stacks.version,names(temp2),Loci)}
        if(length(subs)>0){Output <- c(stacks.version,names(temp2)[-which(names(temp2)%in%c(subs))],Loci)}
      }

      if(keep)
      {
        if(length(subs)==0){Output <- c(stacks.version,names(temp2),Loci)}
        if(length(subs)>0){Output <- c(stacks.version,subs,Loci)}
      }

    # Save the file
    write.table(Output,path,col.names=FALSE,row.names=FALSE,quote=FALSE)

} #End function
