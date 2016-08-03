# Subset Genepop
#' @title Genepop subset loci and populations
#' @description Function to subset loci and populations
#' @param genepop the genepop data to be manipulated. This can be either a file path
#' or a dataframe read in with tab separation, header=FALSE , quote="", and stringsAsFactors=FALSE.
#' This will be the standard genepop format with the first n+1 rows corresponding to the n loci names,
#' or a single comma delimited row of loci names followed by the locus data. Populations are
#' separated by "Pop". Each individual ID is linked to the locus data by " ,  " (space,space space) and is read in as
#' as a single row (character).
#' @param subs the loci names of interest or a vector which corresponds to the order of which
#' they appear in the genepop file.
#' These can be either the order by which they occur or the exact name of the loci
#' e.g. subs <-c(1,2,3,4) would return the first 4 loci
#' & subs <- c("190_56","145_21",456_12") would return loci with defined names.
#' @param keep logical vector which defines whether you want to remove the loci or keep them.
#' The default is to keep them keep <- TRUE assuming you are removing neutral markers
#' and only keeping the subs.
#' @param spop is the populations of interest. Note these are specified in the order which they appear in the
#'  original genepop file. i.e. first pop = 1 second pop = 2
#'  Examples: numeric - spop <- c(1,3,4,7) or
#'  the population ID (alpha-numeric code before the underscore). Here we assume conventional
#'  naming of "Population_sample#" e.g. (Aqua01_05: population Aqua01 & sample #5).
#'            text- spop <- c("Aqua01", "GRR","GHR","TRS").
#' @param path the filepath and filename of output.
#' @rdname subset_genepop
#' @importFrom data.table fread as.data.table
#' @importFrom utils write.table
#' @export


##
subset_genepop <- function(genepop,subs=NULL,keep=TRUE,spop=NULL,path)
  {

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

  ## Now add the population tags using npops (number of populations and Pops for the inter differences)
    tPops <- c(Pops,NROW(genepop))
    PopIDs <- NULL
    for (i in 2:length(tPops)){
      hold <- tPops[i]-tPops[i-1]-1
      if(i==length(tPops)){hold=hold+1}
      pophold <- rep(npops[i-1],hold)
      PopIDs <- c(PopIDs,pophold)
    }

    temp2$Pop <- PopIDs;rm(hold,pophold,tPops,PopIDs)

## Now subset out the the data according to the specified loci and whether or not you want to keep them.

    if(is.numeric(subs))
    { #column number instead of name depending on the output from Outlier detection

      if(!keep) # neutral
      {
        if(length(subs)>0){reqCols <- temp2[,-subs]}
        if(length(subs)==0){reqCols <- temp2}
      }


      if(keep) # outliers or loci under divergent selection
      {
        PopInd=which(names(temp2)=="Pop")
        if(length(subs)>0){reqCols <- temp2[,c(subs,PopInd)]}
        if(length(subs)==0){reqCols <- temp2}
      }

    }

    if(!is.numeric(subs))
    { #column name

      if(!keep)# neutral
      {
        if(length(subs)>0){reqCols <- temp2[,-which(names(temp2)%in%subs)]}
        if(length(subs)==0){reqCols <- temp2}
      }

      if(keep)# outliers or loci under divergent selection
      {
        if(length(subs)>0){reqCols <- temp2[,c(subs,"Pop")]}
        if(length(subs)==0){reqCols <- temp2}
      }
    }

## Now subset the rows
    # is a population subset required
    if(length(spop)>0){

      if(sum(is.numeric(spop))>0){ # if the subsetted populations are numeric
      ind <- which(reqCols$Pop %in% spop) # index where the populations are
      reqCols <- reqCols[ind,]
      temp <- temp[ind,]
      temp2 <- temp2[ind,]
      NamePops <- NamePops[ind]
      }

      if(sum(is.numeric(spop))==0){ # if the subsetted populations are character indexes
        reqCols <- reqCols[which(NameExtract %in% spop),]
        temp <- temp[which(NameExtract %in% spop),]
        temp2 <- temp2[which(NameExtract %in% spop),]
        NamePops <- NamePops[which(NameExtract %in%spop)]
      }


    } # end of subset population if statement

    reqCols <- reqCols[,-length(reqCols)] # delete the "Pop" data * last column

#Now recompile the GenePop format

    #the number of individuals for all populations but the last (Pop tagged to the end)
    if(length(table(temp2$Pop))>1){PopLengths <- table(temp2$Pop)[-length(table(temp2$Pop))]} else {PopLengths=1}

    if(length(table(temp2$Pop))==2){PopPosition = PopLengths+1}

    if(length(table(temp2$Pop))>2){
          PopPosition <- c(PopLengths[1]+1,rep(NA,(length(PopLengths)-1)))
          for (i in 2:length(PopLengths)){
            PopPosition[i] <- PopLengths[i]+PopPosition[i-1]
          }
    }

    # paste together the Loci as one long integer separated for each loci by a space
    Loci <- do.call(paste,c(reqCols[,], sep=" "))

    #Grab the Population tags that each invididual had following the format ID_,__
    PopVec <- paste0(NamePops," ,  ")

    #Paste these to the Loci
    Loci <- paste(PopVec,Loci,sep="")

    #Insert the value of "Pop" which partitions the data among populations #only if more than one population
    if(length(table(temp2$Pop))!=1){Loci <- insert_vals(Vec=Loci,breaks=PopPosition,newVal="Pop")}

    #Add the first "Pop" label
    Loci <- c("Pop",Loci)

    ## Add the column labels and the stacks version

    if(is.numeric(subs))
      { #Column numbers
        if(!keep)
        {
          PopInd=which(names(temp2)=="Pop")
          if(length(subs)==0){Output <- c(stacks.version,names(temp2)[-PopInd],Loci)}
          if(length(subs)>0){Output <- c(stacks.version,names(temp2)[-c(subs,PopInd)],Loci)}
        }

        if(keep)
        {
          PopInd=which(names(temp2)=="Pop")
          if(length(subs)==0){Output <- c(stacks.version,names(temp2)[-PopInd],Loci)}
          if(length(subs)>0){Output <- c(stacks.version,names(reqCols),Loci)}
        }
    }

    if(!is.numeric(subs))
    { # column names
      if(!keep)
      {
        PopInd=which(names(temp2)=="Pop")
        if(length(subs)==0){Output <- c(stacks.version,names(temp2)[-PopInd],Loci)}
        if(length(subs)>0){Output <- c(stacks.version,names(temp2)[-which(names(temp2)%in%c(subs,"Pop"))],Loci)}
      }

      if(keep)
      {
        PopInd=which(names(temp2)=="Pop")
        if(length(subs)==0){Output <- c(stacks.version,names(temp2)[-PopInd],Loci)}
        if(length(subs)>0){Output <- c(stacks.version,subs,Loci)}
      }
    }

    # Save the file
    utils::write.table(Output,path,col.names=FALSE,row.names=FALSE,quote=FALSE)

} #End function
