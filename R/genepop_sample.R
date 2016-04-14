# Subset Genepop Training
#' @title Randomly sample individuals from Genepop.
#' @description Stratified random sample of individuals from a Genepop file.
#' @param GenePop the genepop data to be manipulated. This can be either a file path
#' or a dataframe read in with tab separation, header=FALSE , quote="", and stringsAsFactors=FALSE.
#' This will be the standard genepop format with the first n+1 rows corresponding to the n loci names,
#' or a single comma delimited row of loci names followed by the locus data. Populations are
#' separated by "Pop". Each individual ID is linked to the locus data by " ,  " (space,space space) and is read in as
#' as a single row (character).
#' @param nsample vector or dataframe which defines the random stratified sampling.
#' If nsample is an integer (e.g. 5) then nsample individuals will be selected from each population.
#' If nsample is a fraction (0-0.9) then that percentage of individuals will be selected from each population.
#' If nsample is a dataframe then the number or fraction associated with each population will be sampled.
#' where the first column is the population and the second is the number to be sampled.
#' @rdname genepop_sample
#' @import magrittr
#' @importFrom data.table fread as.data.table
#' @importFrom dplyr filter sample_n sample_frac group_by ungroup
#' @export

##
genepop_sample <- function(GenePop,nsample){

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
  Individuals <- NamePops
  NameExtract <- substr(NamePops,1,regexpr("_",NamePops)-1)

    #create a dataframe which will be sampled
    Popdf <- data.frame(Pops=NameExtract,Indiv=NamePops)
    Popdf$Pops <- as.character(Popdf$Pops)
    Popdf$Indiv <- as.character(Popdf$Indiv)

    ## Random stratified sample

    #fixed number from each pop *note we use a NULL filter so we can call dplyr::
    if(!length(nsample)>1){
      if(as.integer(nsample)>0){
        sampledf <- as.data.frame(dplyr::filter(Popdf)%>%group_by(Pops)%>%sample_n(nsample)%>%ungroup())
        sub_individuals <- sampledf$Indiv
      }
    }

    #sample fixed proportions
   if(!length(nsample)>1){
     if(!as.integer(nsample)>0){
      if(nsample>0.9)
        {
          nsample <- 0.9#maximum sample of 0.9
          warning("Maximum proportion permitted is 90% nsample changed to 0.9")
        }
      sampledf  <- as.data.frame(dplyr::filter(Popdf)%>%group_by(Pops)%>%sample_frac(nsample)%>%ungroup())
      sub_individuals <- sampledf$Indiv
     }
   }

    #sample variable fixed amounts
    if(length(nsample)>1){
      if(as.integer(nsample[1,2])>0){
      sub_individuals=NULL
      for (i in 1:nrow(nsample))
      {
        temp <- filter(Popdf,Pops==nsample[i,1])
        temp <- as.data.frame(filter(temp)%>%group_by(Pops)%>%sample_n(nsample[i,2])%>%ungroup())
        sub_individuals=c(sub_individuals,temp$Indiv)
        }
      }
    }

    #sample variable fixed percentages
    if(length(nsample)>1){
      if(!as.integer(nsample[1,2])>0){
      sub_individuals=NULL
      for (i in 1:nrow(nsample))
        {
          if(nsample[i,2]>0.9)
          {
            nsample[i,2] <- 0.9#maximum sample of 0.9
            warning(paste("Maximum proportion permited is 90% nsample changed to 0.9 for population",nsample[i,1],sep=" "))
          }
          temp <- filter(Popdf,Pops==nsample[i,1])
          temp <- as.data.frame(filter(temp)%>%group_by(Pops)%>%sample_frac(nsample[i,2])%>%ungroup())
          sub_individuals=c(sub_individuals,temp$Indiv)
        }
      }
    }

    return(sub_individuals)

} #End function
