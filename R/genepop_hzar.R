# GenePop -> HZAR
#' @title Convert Genepop to HZAR format.
#' @description Function to convert Genepop to the R package HZAR.
#' @param genepop the genepop data to be manipulated. This can be either a file path
#' or a dataframe read in with tab separation, header=FALSE , quote="", and stringsAsFactors=FALSE.
#' This will be the standard genepop format with the first n+1 rows corresponding to the n loci names,
#' or a single comma delimited row of loci names followed by the locus data. Populations are
#' separated by "Pop". Each individual ID is linked to the locus data by " ,  " (space,space space) and is read in as
#' as a single row (character).
#' @param distances A dataframe or path to a text file with your distances between populations. Should contain 2 columns -
#' Populations and Distances.There should be the same number of populations as in the Genepop file.
#' @param path the filepath and filename of output.
#' @rdname genepop_hzar
#' @import magrittr
#' @importFrom data.table fread as.data.table melt dcast
#' @importFrom dplyr filter summarise group_by ungroup summarise_each funs funs_
#' @export
#'

genepop_hzar<-function(genepop,distances,path){

  if(is.character(distances)){distances<-read.table(distances,header = TRUE)}

  colnames(distances)=c("Pop","Distance")

  if (length(distances[,1])!=length(genepop_detective(genepop))){
    stop("Number of distances does not match the number of populations present")
    }

  #First get allele frequencies for each population
  allelefreqs<-genepop_allelefreq(genepop,wide=TRUE)

  #Check to see if genepop is a data.frame from the workspace
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

  ## Get the population names (prior to the _ in the Sample ID)
  NamePops <- temp[,1] # Sample names of each
  NameExtract <- substr(NamePops,1,regexpr("_",NamePops)-1)

  PopNum <- data.frame(table(NameExtract))
  colnames(PopNum)[1] <- "Population"

  #allele coding length
  alleleEx <- max(sapply(temp2[,1],FUN=function(x){nchar(as.character(x[!is.na(x)]))})) #presumed allele length

  #check to make sure the allele length is a even number
  if(!alleleEx %% 2 ==0){stop(paste("The length of each allele is assumed to be equal (e.g. loci - 001001 with 001 for each allele), but a max loci length of", alleleEx, "was detected. Please check data."))}

  #get the allele values summary header
  firstAllele <-  as.data.frame(sapply(temp2,function(x)as.numeric(as.character(substring(x,1,alleleEx/2)))))
  secondAllele <-  as.data.frame(sapply(temp2,function(x)as.numeric(as.character(substring(x,(alleleEx/2)+1,alleleEx)))))

  #Add pop labels for grouping
  firstAllele$Pop <- NameExtract
  secondAllele$Pop <- NameExtract

  temp2$Pops <- NameExtract
  temp2[] <- lapply(temp2, as.character)

  ## Identify major alleles for each loci
  mallele <- temp2[,-length(temp2)]%>%summarise_each(funs(AlleleFreq(.)))
  majorfreqs <- as.vector(t(mallele[1,]))
  majordf <- data.frame(variable=names(temp2[-length(temp2)]),major=majorfreqs)

  zerocount=function(x){length(x)-length(which(x==0))}
  Count1 <- firstAllele%>%group_by(Pop)%>%summarise_each(funs(zerocount))%>%ungroup()%>%data.frame()
  Count2 <- secondAllele%>%group_by(Pop)%>%summarise_each(funs(zerocount))%>%ungroup()%>%data.frame()
  SumCount <- Count1[,-1]+Count2[,-1]
  colnames(SumCount)=paste("Count",colnames(SumCount),sep="_")
  SumCount$Pop <- Count1$Pop

  ## Add it all together
  #set factor levels to the same for reordering
  distances <- distances[order(distances$Distance),]
  distances$Pop <- factor(distances$Pop,levels=distances$Pop)
  allelefreqs$Pop <- factor(allelefreqs$Pop,levels=distances$Pop)
  SumCount$Pop <- factor(SumCount$Pop,levels=distances$Pop)

  #reorder allelefreq and SumCount
  allelefreqs <- allelefreqs[order(allelefreqs$Pop),]
  SumCount <- SumCount[order(SumCount$Pop),]

  Output <- cbind(distances,
                  allelefreqs[,-grep("Pop",colnames(allelefreqs))],
                  SumCount[,-grep("Pop",colnames(SumCount))])
  colnames(Output)[1:2]=c("Population","Distance")


  write.table(Output,file = path,quote = FALSE,row.names =FALSE, sep="\t")

}
