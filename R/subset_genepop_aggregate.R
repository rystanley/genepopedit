# Subset Genepop Aggregate
#' @title Genepop subset and combine/rename populations
#' @description Function for the manipulation of genopop format SNP datasets and renaming of populations
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
#' @param subs he loci names of interest or a vector which corresponds the the order of which
#' they appear in the genepop file.
#' These can be either the order by which they occur or the exact name of the loci
#' e.g. subs <- c(1,2,3,4) would return the first 4 loci &
#' subs <- c("190-56","145_21",456_12") would return loci with defined names.
#' @param keep logical vector which defines whether you want to remove the loci or keep them.
#' the default is to keep them (keep <- TRUE) assuming you are removing neutral markers
#' and only keeping the subs
#' @param dirname directory where the output file will be saved.
#' @param agPopFrame a dataframe or path to a csv. This must be specified in this function
#' this dataframe contains two columns. Column 1 corresponds to the population names. These names
#' should match the individual IDs (e.g. BON01  , 110110 120120 -- would be 'BON'). The next column
#' has the new names that you want to replace. If you don't want to change the name then just repeat
#' from column one in that row. **Note if this is a reference to a path then there should be column headers
#' though the headers do not need to match this example. If the input is a dataframe object from the
#' workspace it must be a data.frame object and therefore will have headers.
#' e.g.
#' Opop   AgPop
#' BON    BON
#' BRA    BON
#' EDN    EDN
#' CRA    CRA
#' MAL    BON
#' TRY    CRA
#'
#' @rdname subset_genepop_aggregate
#' @importFrom tidyr separate
#' @importFrom  stringr str_extract
#' @export

subset_genepop_aggregate <- function(GenePop,subs=NULL,keep=TRUE,dirname,agPopFrame){

## Stacks version information
    stacks.version <- GenePop[1,] # this could be blank or any other source. First row is ignored by GenePop

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
    if (length(temp2)!=length(ColumnData)){stacks.version="No STACKS version specified"}
    if (length(temp2)==length(ColumnData)){colnames(temp2) <- ColumnData}

## Get the Alpha names from the
    NamePops=temp[,1] # Names of each
    NameExtract=stringr::str_extract(NamePops, "[A-Z]+" ) # extract the text from the individuals names to denote population

## Now add the population tags using npops (number of populations and Pops for the inter differences)
     tPops <- c(Pops,NROW(GenePop))
      PopIDs <- NULL
          for (i in 2:length(tPops)){
            hold <- tPops[i]-tPops[i-1]-1
            if(i==length(tPops)){hold=hold+1}
            pophold <- rep(npops[i-1],hold)
            PopIDs <- c(PopIDs,pophold)
          }

    temp2$Pop <- PopIDs;rm(hold,pophold,tPops)

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

## Now subset the Populations

    if(!is.data.frame(agPopFrame)) #if it isn't a dataframe then read in the path
    {
      agPopFrame <- read.csv(agPopFrame,header=T)
    }

    sPop <- as.character(agPopFrame[,1]) # these are the populations of interest


    # is a population subset required
      if(sum(is.numeric(sPop))>0){ # if the subsetted populations are numeric
      ind <- which(reqCols$Pop %in% sPop) # index where the populations are
      reqCols <- reqCols[ind,]
      temp <- temp[ind,]
      temp2 <- temp2[ind,]
      }

      if(sum(is.numeric(sPop))==0){ # if the subsetted populations are character indexes
        reqCols <- reqCols[which(NameExtract %in% sPop),]
        temp <- temp[which(NameExtract %in% sPop),]
        temp2 <- temp2[which(NameExtract %in% sPop),]
      }

    #Now recompile the GenePop format

    NameExtract2 <- NameExtract[which(NameExtract %in% sPop)]
    #Create a new vector with the new population names
    for (i in sPop){
      NameExtract2[which(NameExtract2 == i)]=as.character(agPopFrame[which(agPopFrame[,1]==i),2])
    }

    NameExtract3 <- NameExtract2
    for (i in 1:length(unique(NameExtract2))){
      NameExtract3[which(NameExtract3==unique(NameExtract3)[i])]=i
    }

    #the number of individuals for all popualtions but the last (Pop tagged to the end)
    PopLengths <- table(factor(NameExtract3, levels=unique(NameExtract3)))[-length(table(NameExtract3))]

    if(length(table(NameExtract3))==2){PopPosition = PopLengths+1}

    if(length(table(NameExtract3))>2){
      PopPosition <- c(PopLengths[1]+1,rep(NA,(length(PopLengths)-1)))
      for (i in 2:length(PopLengths)){
        PopPosition[i] <- PopLengths[i]+PopPosition[i-1]
      }
    }

    #Now stitch the data together
    # paste together the Loci as one long integer seperated for each loci by a space
    reqCols <- reqCols[,-length(reqCols)]
    Loci <- do.call(paste,c(reqCols[,], sep=" "))

    #Grab the Population tags that each invididual had following the format ID_,__
    popvec1 <- unlist(strsplit(gsub(pattern="_",replacement="",temp[,1]),split = "[^0-9]+"))
    popvec2 <- popvec1[which(popvec1 != "")]

    PopVec <- paste0(NameExtract2,"_",popvec2," ,  ")

    #Paste these to the Loci
    Loci <- paste(PopVec,Loci,sep="")

    #Insert the value of "Pop" which partitions the data among populations #only if more than one population
    if(length(table(NameExtract2))!=1){Loci <- insert_vals(Vec=Loci,breaks=PopPosition,newVal="Pop")}

    #Add the first "Pop" label
    Loci <- c("Pop",Loci)

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
        if(length(subs)==0){Output <- c(stacks.version,names(temp2)[-length(names(temp2))],Loci)}
        if(length(subs)>0){Output <- c(stacks.version,names(temp2)[-which(names(temp2)%in%c(subs,"Pop"))],Loci)}
      }

      if(keep)
      {
        if(length(subs)==0){Output <- c(stacks.version,names(temp2)[-length(names(temp2))],Loci)}
        if(length(subs)>0){Output <- c(stacks.version,subs,Loci)}
      }
    }

    # Save the file
    write.table(Output,dirname,col.names=FALSE,row.names=FALSE,quote=FALSE)

} #End function
