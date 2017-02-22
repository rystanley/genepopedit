# Subset Genepop Rename
#' @title Genepop subset and rename populations
#' @description Function for the manipulation of genopop format SNP datasets and renaming of populations
#' @param genepop the genepop data to be manipulated. This can be either a file path
#' or a dataframe read in with tab separation, header=FALSE , quote="", and stringsAsFactors=FALSE.
#' This will be the standard genepop format with the first n+1 rows corresponding to the n loci names,
#' or a single comma delimited row of loci names followed by the locus data. Populations are
#' separated by "Pop". Each individual ID is linked to the locus data by " ,  " (space,space space) and is read in as
#' as a single row (character).
#' @param nameframe a dataframe or path to a csv and is required for this function. The first column of the dataframe defines the original value
#' and the second corresponds to the change. If populations (default: meta="populations") are the metadata to be changed
#' the names should match the individual IDs (e.g. BON_01 ,  110110 120120 = 'BON'). The next column
#' has the new names that you want to change or rename. If you don't want to change the name then just repeat
#' from column one in that row. This function assumes that each population will have a unique name. Names in this case are
#' comprised of alpha characters and not numbers.
#' (e.g. Pop01_01 and Pop02_01 would each be considered 'Pop' for the population name)
#' e.g. data.frame(Opop=c("BON","BRA","EDN","CRA","MAL","TRY"),Rename=c("BON","BON","EDN","CRA","BON","CRA")).
#' @param renumber is a logical (default=FALSE) defining whether you want to change the sample unique identity
#' i.e. sample number - BON_01 where 01 is the unique qauntity. If multiple populations are combined this will
#' prevent two samples from having the same name.
#' @param meta which metadata will be renamed with nameframe. Default is "Pop" for populations, alternative is "Ind" for individuals.
#' @param path the filepath and filename of output.
#' @rdname subset_genepop_rename
#' @importFrom data.table fread as.data.table
#' @importFrom utils read.csv
#' @importFrom utils write.table
#' @export

subset_genepop_rename <- function(genepop,nameframe,renumber=FALSE,meta="Pop",path){

  if(length(which(meta %in% c("Pop","Ind")))==0){
    stop("Parameter 'meta' must be defined as Pop or Ind for renaming either the population or individual ID. Function stopped.",call. = FALSE)
  }

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
  NamePops <- temp[,1] # Sample names of each
  NameExtract <- substr(NamePops,1,regexpr("_",NamePops)-1)

## Now change the population names

  if(meta == "Pop"){

    if(!is.data.frame(nameframe)) #if it isn't a dataframe then read in the path
    {
      nameframe <- utils::read.csv(nameframe,header=T)
    }

    sPop <- as.character(nameframe[,1]) # these are the populations of interest

    #Now recompile the GenePop format

    NameExtract2 <- NameExtract
    #Create a new vector with the new population names
    for (i in sPop){
      NameExtract2[which(NameExtract2 == i)]=as.character(nameframe[which(nameframe[,1]==i),2])
    }

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

    #Grab the Population tags that each individual had following the format ID_,__
    popvec1 <- unlist(strsplit(as.character(temp[,1]),split = "_"))
    popvec2 <- popvec1[!popvec1%in%unique(NameExtract)]

    if(renumber)
    {
      popvec2=NULL
      for (i in 1:length(table(NameExtract2)))
      {
        popvec2=c(popvec2,sub('^(.)$', '0\\1', 1:table(NameExtract2)[which(names(table(NameExtract2))==unique(NameExtract2)[i])]))
      }
    }

    PopVec <- paste0(NameExtract2,"_",popvec2," ,  ")

    #Paste these to the Loci
    Loci <- paste(PopVec,Loci,sep="")

    #Insert the value of "Pop" which partitions the data among populations #only if more than one population
    if(length(table(NameExtract2))!=1){Loci <- insert_vals(Vec=Loci,breaks=PopPosition,newVal="Pop")}

  }# end if (meta=="Pop")

  if(meta=="Ind")
    {
    if(!is.data.frame(nameframe)) #if it isn't a dataframe then read in the path
    {
      nameframe <- utils::read.csv(nameframe,header=T)
    }

    nameframe[]=lapply(nameframe,as.character)
    NamePops=as.character(NamePops)

    #get the row indicies in the right order.
    nameframe$ind=sapply(nameframe[,1],FUN=function(x){which(NamePops == x)})

    NamePops[nameframe$ind] <- nameframe[,2] #replace values

    NameExtract <- substr(NamePops,1,regexpr("_",NamePops)-1)

    #Now stitch the data together
    # paste together the Loci as one long integer separated for each loci by a space
    Loci <- do.call(paste,c(temp2[,], sep=" "))

    #Paste these to the Loci
    Loci <- paste(NamePops,Loci,sep=" ,  ")

    PopLengths <- table(factor(NameExtract, levels=unique(NameExtract)))[-length(table(NameExtract))]

    if(length(table(NameExtract))==2){PopPosition = PopLengths+1}

    if(length(table(NameExtract))>2){
      PopPosition <- c(PopLengths[1]+1,rep(NA,(length(PopLengths)-1)))
      for (i in 2:length(PopLengths)){
        PopPosition[i] <- PopLengths[i]+PopPosition[i-1]
      }
    }

    #Insert the value of "Pop" which partitions the data among populations #only if more than one population
    if(length(table(NameExtract))!=1){Loci <- insert_vals(Vec=Loci,breaks=PopPosition,newVal="Pop")}

  } # end of if (meta == "Ind")


   #Add the first "Pop" label
    Loci <- c("Pop",Loci)

    Output <- c(stacks.version,names(temp2),Loci)

    # Save the file
    utils::write.table(Output,path,col.names=FALSE,row.names=FALSE,quote=FALSE)

} #End function
