# Subset Genepop Training
#' @title Randomly sample individuals from Genepop.
#' @description Stratified random sample of individuals from a Genepop file.
#' @param GenePop the genepop data to be manipulated. This can be either a file path
#' or a dataframe read in with tab seperation, header=FALSE , quote="", and stringsAsFactors=FALSE.
#' This will the standard genepop format with a the first n+1 rows corresponding the the n loci names,
#' or a single commma deliminated row of loci names followed by the locus data. Populations are
#' seperated by "Pop". Each individual ID is linked the the locus data by "  , " and is read in as
#' as single row (character).
#' @param nsample vector or dataframe which defines the random stratified sampling
#' if nsample is an integer (e.g. 5) then nsample individuals will be selected from each population.
#' if nsample is a fraction (0-0.9) then that percentage of individuals will be selected from each population.
#' if nsample is a dataframe then the number or fraction associated with each population will be sampled.
#' where the first column is the population and the second is the number to be sampled.
#' @rdname genepop_sample
#' @import magrittr
#' @importFrom tidyr separate
#' @importFrom dplyr filter
#' @importFrom dplyr sample_n
#' @importFrom dplyr sample_frac
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr sample_frac
#' @export

##
genepop_sample <- function(GenePop,nsample){

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

    #list of Individuals
    Individuals <- temp$Pops
    Individuals <- gsub(" ","",temp$Pops)

    ## Get the population names (prior to the _ in the Sample ID)
    NamePops <- temp[,1] # Sample names of each
    NamePops <- gsub(" ","",NamePops) #get rid of space
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
          warning("Maximum proportion permited is 90% nsample changed to 0.9")
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
