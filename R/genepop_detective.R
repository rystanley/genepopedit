# GenePop detective
#' @title Explore Genepop data structure.
#' @description Function returns Genepop file meta-data.
#' @param genepop the genepop data to be manipulated. This can be either a file path
#' or a dataframe read in with tab separation, header=FALSE , quote="", and stringsAsFactors=FALSE.
#' This will be the standard genepop format with the first n+1 rows corresponding to the n loci names,
#' or a single comma delimited row of loci names followed by the locus data. Populations are
#' separated by "Pop". Each individual ID is linked to the locus data by " ,  " (space, space space) and is read in as
#' as a single row (character).
#' @param variable data to be returned
#' Four options \code{default = "Pops"}
#' "Pops" = vector of population names.
#' "PopNum" = dataframe of population names and counts.
#' "Inds" = vector of sample IDs.
#' "Loci" = vector of Loci.
#' "Allele" = vector of allele values.
#' @rdname genepop_detective
#' @importFrom data.table fread
#' @export
#'

genepop_detective <- function(genepop,variable="Pops"){

  #check the variable parameter
  if(length(which(variable %in% c("Pops","PopNum","Inds","Loci","All","Allele")))==0){
    stop("Parameter 'variable' defining metadata to be returned must be defined as one of Pop, Inds, Loci, PopNum, Allele, or All. Function stopped.",call. = FALSE)
  }

  #Check to see if genepop is a data.frame from the workspace
  if(is.data.frame(genepop)){genepop <- as.data.table(genepop)}

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
    genepop <- data.frame(genepop,stringsAsFactors = FALSE)
  }

  ## Stacks version information
  stacks.version <- genepop[1,] #this could be blank or any other source. First row is ignored by genepop

  #Remove first label of the stacks version
  genepop <- genepop[-1,]
  colnames(genepop) <- "data"

  #ID the rows which flag the Populations
  Pops  <-  which(genepop$data == "Pop" | genepop$data =="pop" | genepop$data == "POP")
  npops  <-  1:length(Pops)

  ## Seperate the data into the column headers and the rest
  ColumnData <- genepop$data[1:(Pops[1]-1)]
  ColumnData <- gsub("\r","",ColumnData)#remove any hidden carriage returns
  snpData <- genepop[Pops[1]:NROW(genepop),]

  #Get a datafile with just the snp data no pops
  tempPops <- which(snpData$data=="Pop"| snpData$data =="pop" | snpData$data == "POP") ## Changed because we allowed
  snpData <- snpData[-tempPops,]

  #Seperate the snpdata
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

  ## Get the population names (prior to the _ in the Sample ID)
  NamePops <- temp[,1] # Sample names of each
  NameExtract <- substr(NamePops,1,regexpr("_",NamePops)-1)

  PopNum <- data.frame(table(NameExtract))
  colnames(PopNum)[1] <- "Population"

  #convert the snp data into character format to get rid of factor levels
  #temp2[] <- lapply(temp2, as.character)

  if(variable=="Allele"){

    #allele coding length
    alleleEx <- max(sapply(temp2[,1],FUN=function(x){nchar(as.character(x[!is.na(x)]))})) #presumed allele length

    #check to make sure the allele length is a even number
    if(!alleleEx %% 2 ==0){stop(paste("The length of each allele is assumed to be equal (e.g. loci - 001001 with 001 for each allele), but a max loci length of", alleleEx, "was detected. Please check data."))}

    #get the allele values summary header
    firstAllele <-  as.data.frame(sapply(temp2,function(x)as.numeric(as.character(substring(x,1,(alleleEx/2))))))
    secondAllele <-  as.data.frame(sapply(temp2,function(x)as.numeric(as.character(substring(x,(alleleEx/2)+1,alleleEx)))))

    Allele <- unique(as.numeric(unique(unlist(firstAllele)))
                     ,as.numeric(unique(unlist(secondAllele))))

    Allele <- Allele[order(Allele)]} #sort the Allele values (NA or 0 will be first)

  if(variable=="All"){
    #create an list container for the data
    Output <- list()
    Output$Pops <- as.character(unique(NameExtract))
    Output$Loci=as.character(names(temp2))
    Output$Inds=as.character(NamePops)
    Output$ PopNum=PopNum
  }

  #return data vector of interest
  if(variable=="Pops"){return(as.character(unique(NameExtract)))}
  if(variable=="PopNum"){return(PopNum)}
  if(variable=="Inds"){return(as.character(NamePops))}
  if(variable=="Loci"){return(as.character(names(temp2)))}
  if(variable=="All"){return(Output)}
  if(variable=="Allele"){return(Allele)}

}
