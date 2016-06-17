# Genepop -> GSI sim
#' @title Convert Genepop to GSI sim format.
#' @description Covert from GENEPOP to format required by gsi_sim. Note that output has SampleIDs formated as Population_Population_ID instead of conventioanl Population_ID to fit commong sampleID naming approach of genepopedit into the convention of gsi_sim.
#' @param GenePop the genepop data to be manipulated. This can be either a file path
#' or a dataframe read in with tab separation, header=FALSE , quote="", and stringsAsFactors=FALSE.
#' This will be the standard genepop format with the first n+1 rows corresponding to the n loci names,
#' or a single comma delimited row of loci names followed by the locus data. Populations are
#' separated by "Pop". Each individual ID is linked to the locus data by " ,  " (space, space space) and is read in as
#' as a single row (character).
#' @param path the filepath and filename of output.
#' @rdname genepop_GSIsim
#' @importFrom data.table fread as.data.table
#' @export
#'

##
genepop_GSIsim <- function(GenePop,path){

  #Check to see if GenePop is a data.frame from the workspace
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
    GenePop <- data.frame(GenePop,stringsAsFactors = FALSE)
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

  ## Get the population names (prior to the _ in the Sample ID)
  NamePops <- temp[,1] # Sample names of each
  NameExtract <- substr(NamePops,1,regexpr("_",NamePops)-1)
  UniquePops <- NameExtract[!duplicated(NameExtract)] #the unique values of each population for the gsisim label

    PopNum <- data.frame(table(NameExtract))
    colnames(PopNum)[1] <- "Population"

    #convert the snp data into character format to get rid of factor levels
    alleleEx <- max(unique(nchar(as.character(temp2[nrow(temp2),2]))),na.rm=T)

    #get the allele values summary header
    firstAllele <-  as.data.frame(sapply(temp2,function(x)as.character(substring(x,1,alleleEx/2))),stringsAsFactors = FALSE)
    secondAllele <-  as.data.frame(sapply(temp2,function(x)as.character(substring(x,(alleleEx/2)+1,alleleEx))),stringsAsFactors = FALSE)

    #Create temporary matrix which will be populated with interlaced even and odd columns
    #matrix dimensions
    rows.combined <- dim(firstAllele)[1]
    cols.combined <- dim(firstAllele)[2]*2 #because alleles for each SNP are seperated

    #create NULL frame to be populated
    AlleleFrame <- as.data.frame(matrix(NA, nrow=rows.combined, ncol=cols.combined))

    #populate the frame
    AlleleFrame[,seq(1,cols.combined,2)] <- firstAllele #first allele
    AlleleFrame[,seq(2,cols.combined,2)] <- secondAllele #second allele

    #The structure of gsi_sim has two underscores in front of each unique sample ID. Because in the case of
    #genepopedit all sample IDs are separeted from population ID's by a single _, the output of genepop_GSIsim will
    #have Population_Population_sampleID

    NamePops <- paste(NameExtract,NamePops,sep="_")

    AlleleFrame <- cbind(NamePops,AlleleFrame)

    ## population grouping variables
    pPops <- NULL
    for (i in 2:length(tempPops)){
      pPops <- c(pPops,tempPops[i]-(i-1))
    }

    pPops <- c(1,pPops)

    # paste together the Loci as one long integer separated for each loci by a space
    Loci <- do.call(paste,c(AlleleFrame[,], sep=" "))

    ##Add in the population seperation values
    if(length(npops)!=1){Loci <- insert_vals(Vec=Loci,breaks=pPops,newVal="Pop")}

    #Add the first Pop label if only one population
    if(length(npops)==1){Loci=c("Pop",Loci)}

    #create the seperation variables
    PopSeps=paste("Pop",UniquePops,sep=" ")

    #replace existing Pop labels with the modified "Pop + Population name" labels.
    Loci[which(Loci=="Pop")]=PopSeps

    #Add the individual count, loci count, Loci Names and allele calls.
    Output <- c(paste(nrow(temp2),length(names(temp2)),sep=" "),names(temp2),Loci)

    # Save the file
    write.table(Output,path,col.names=FALSE,row.names=FALSE,quote=FALSE)

}
