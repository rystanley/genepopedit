# Subset samples from Genepop I
#' @title Genepop remove or keep specific sample IDs
#' @description Function for the manipulation of genopop format SNP datasets
#' @param GenePop = the genepop data to be manipulated. This can be either a file path
#' or a dataframe read in with tab seperation, header=FALSE , quote="", and stringsAsFactors=FALSE.
#' This will be the standard genepop format with the first n+1 rows corresponding to the n loci names,
#' or a single commma delimited row of loci names followed by the locus data. Populations are
#' seperated by "Pop". Each individual ID is linked to the locus data by "  , " and is read in as
#' as single row (character).
#' @param indiv vector sample IDs of interest.
#' These can be either the order by which they occur or the exact name of the loci
#' indiv <- \code{c("Pop01_01","Pop03_15","Pop16_02")} would individuals with these sample names.
#' @param keep logical whether to delete sample IDs specified by indiv (default: TRUE) or delete all other IDs.
#' @param path the filepath and filename of output
#' @rdname subset_genepop_individual
#' @importFrom tidyr separate
#' @export

##
subset_genepop_individual <- function(GenePop,indiv=NULL,keep=FALSE,path){

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

    #list of Individuals
    Individuals <- temp$Pops
    Individuals <- gsub(" ","",temp$Pops)

    #Contingency to see if R read in the top line as the "stacks version"
    if (length(temp2)!=length(ColumnData)){colnames(temp2) <- c(stacks.version,ColumnData)}
    if (length(temp2)!=length(ColumnData)){stacks.version="No STACKS version specified"}
    if (length(temp2)==length(ColumnData)){colnames(temp2) <- ColumnData}

    ## Get the population names (prior to the _ in the Sample ID)
    NamePops <- temp[,1] # Sample names of each
    NamePops <- gsub(" ","",NamePops) #get rid of space
    NameExtract <- substr(NamePops,1,regexpr("_",NamePops)-1)

## now subset the individuals

    if(keep){
        temp <- temp[which(Individuals %in% indiv),]
        temp2 <- temp2[which(Individuals %in% indiv),]
        NameExtract2 <- NameExtract[which(Individuals %in% indiv)]
    }

    if(!keep){
      temp <- temp[which(Individuals %in% setdiff(Individuals,indiv)),]
      temp2 <- temp2[which(Individuals %in% setdiff(Individuals,indiv)),]
      NameExtract2 <- NameExtract[which(Individuals %in% setdiff(Individuals,indiv))]
    }

    #Now recompile the GenePop format

    #Create a new vector with the new population names
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
    Loci <- do.call(paste,c(temp2[,], sep=" "))

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

    #Derive the genepop output
    Output <- c(stacks.version,names(temp2),Loci)

    # Save the file
    write.table(Output,path,col.names=FALSE,row.names=FALSE,quote=FALSE)

} #End function
