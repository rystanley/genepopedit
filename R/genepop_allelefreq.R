# GenePop allelefreq
#' @title Explore population specific allele frequencies.
#' @description Function returns population derived Major allele frequency (default) or full panel allele frequencies.
#' @param genepop the genepop data to be manipulated. This can be either a file path
#' or a dataframe read in with tab separation, header=FALSE , quote="", and stringsAsFactors=FALSE.
#' This will be the standard genepop format with the first n+1 rows corresponding to the n loci names,
#' or a single comma delimited row of loci names followed by the locus data. Populations are
#' separated by "Pop". Each individual ID is linked to the locus data by " ,  " (space, space space) and is read in as
#' as a single row (character).
#' @param popgroup population grouping using the "Pop" deliminiter (Default: NULL) or a dataframe or path to a csv. The grouping dataframe should have two columns, the first corresponding to the population name and the second to an aggregation vector of common groups. Each population can only be assigned to one group.
#' @param fullpanel (default: FALSE) a logical parameter specifying whether allele frequencies for all loci by population should be retured.
#' @param wide logical specifying whether the allele frequencies should be returned as long (default:FALSE) or wide (TRUE) format. Note that the wide format can be used as the input for alleleotype_genepop to simulate geneotypes.
#' @rdname genepop_allelefreq
#' @import magrittr
#' @importFrom data.table fread as.data.table melt dcast
#' @importFrom dplyr filter summarise group_by ungroup summarise_all funs funs_ do
#' @importFrom tidyr unnest
#' @export
#'

genepop_allelefreq <- function(genepop,popgroup=NULL,wide=FALSE,fullpanel=FALSE){

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

  ## population grouping variables
  pPops <- NULL
  for (i in 2:length(tempPops)){
    pPops <- c(pPops,tempPops[i]-(i-1))
  }

  pPops <- c(1,pPops)

  #count populatons
  lPops=as.vector(table(NameExtract)[unique(NameExtract)])
  PopGroupVec <- rep(1:length(lPops),times=lPops)
  PopGroupVec <- PopGroupVec[order(PopGroupVec)]


if(!fullpanel){

    if(is.null(popgroup)){
      temp2$Pops <- PopGroupVec
      temp2[] <- lapply(temp2, as.character)

      ## Identify major alleles for each loci
      mallele <- temp2[,-length(temp2)]%>%summarise_all(funs(AlleleFreq(.)))
      majorfreqs <- as.vector(t(mallele[1,]))
      majordf <- data.frame(variable=names(temp2[-length(temp2)]),major=majorfreqs)

      tData <- data.table::melt(temp2,id.vars="Pops") #transposed dataframe
      transposeData=merge(tData,majordf,by="variable") # add the major allele for a given locus

      Prec <- transposeData%>%
        group_by(Pops,variable)%>%
        summarise(prec=AlleleFreqLoci(value,unique(major)))%>%
        ungroup()%>%data.frame()

      PopLabels <- data.frame(Names=NameExtract,Pops=as.character(PopGroupVec))
      PopLabels <- PopLabels%>%
        group_by(Pops)%>%
        summarise(Name=paste(unique(Names),collapse="_"))%>%
        ungroup()%>%data.frame()
      PopLabels [] <- lapply(PopLabels , as.character)

      Output <- merge(Prec,PopLabels,by="Pops")
      Output <- Output[,c(4,2,3)]
      colnames(Output) <- c("Population","Loci","MAF")
      Output$Loci=as.character(Output$Loci)

    }

    if(!is.null(popgroup)){
      popgroup[] <- lapply(popgroup, as.character)
      colnames(popgroup)=c("Name","Group")
      temp3 <- temp2[which(NameExtract %in% popgroup[,1]),]
      temp3$Name <- NameExtract[which(NameExtract %in% popgroup[,1])]
      temp4 <- merge(temp3,popgroup,by="Name")

      PopLabels <- data.frame(PopNames=temp4$Name,Group=temp4$Group)
      PopLabels <- as.data.frame(PopLabels%>%
                                   group_by(Group)%>%
                                   summarise(Name=paste(unique(PopNames),collapse="_"))%>%ungroup())

      #temp5 <- merge(temp4[,-grep("Pops",names(temp4))],PopLabels,by="Group")
      temp5 <- merge(PopLabels,temp4[,-grep("Pops|Name",names(temp4))],by="Group")
      #temp5 <- temp5[,c(names(temp2),names(temp5)[length(temp5)])] # select columns
      #colnames(temp5)[length(temp5)]="Pops"
      temp5 <- temp5[,c(setdiff(names(temp2),c("Pops","Group")),"Name")] # select columns
      temp5[] <- lapply(temp5, as.character)

      ## Identify major alleles for each loci
      mallele <- temp5[,-length(temp5)]%>%summarise_all(funs(AlleleFreq(.)))
      majorfreqs <- as.vector(t(mallele[1,]))
      majordf <- data.frame(variable=names(temp5[-length(temp5)]),major=majorfreqs)

      tData <- data.table::melt(temp5,id.vars="Name") #transposed dataframe
      transposeData=merge(tData,majordf,by="variable") # add the major allele for a given locus

      Output <- transposeData%>%
        group_by(Name,variable)%>%
        summarise(prec=AlleleFreqLoci(value,unique(major)))%>%
        ungroup()%>%data.frame()

      colnames(Output)=c("Population","Loci","MAF")
      Output$Loci=as.character(Output$Loci)
    }

} #end if(!fullpanel)

### Full panel allele frequencies -------------------

  if(fullpanel & is.null(popgroup)){
    temp2$Pops <- PopGroupVec
    temp2[] <- lapply(temp2, as.character)

    temp3 <- reshape2::melt(temp2,id.vars="Pops")
    colnames(temp3) <- c("Pops","Loci","Value")

    temp3$Loci <- as.character(temp3$Loci)

    Output <- tidyr::unnest(temp3%>%group_by(Pops,Loci)%>%do(split=frequencies(.))%>%
                             ungroup()%>%data.frame(.,stringsAsFactors=F))
  }


  if(fullpanel & !is.null(popgroup)){
    popgroup[] <- lapply(popgroup, as.character)
    colnames(popgroup)=c("Name","Group")
    temp3 <- temp2[which(NameExtract %in% popgroup[,1]),]
    temp3$Name <- NameExtract[which(NameExtract %in% popgroup[,1])]
    temp4 <- merge(temp3,popgroup,by="Name")

    PopLabels <- data.frame(PopNames=temp4$Name,Group=temp4$Group)
    PopLabels <- as.data.frame(PopLabels%>%
                                 group_by(Group)%>%
                                 summarise(Name=paste(unique(PopNames),collapse="_"))%>%ungroup())

    #temp5 <- merge(temp4[,-grep("Pops",names(temp4))],PopLabels,by="Group")
    temp5 <- merge(PopLabels,temp4[,-grep("Pops|Name",names(temp4))],by="Group")
    #temp5 <- temp5[,c(names(temp2),names(temp5)[length(temp5)])] # select columns
    #colnames(temp5)[length(temp5)]="Pops"
    #temp5 <- temp5[,c(setdiff(names(temp2),c("Pops","Group")),"Name")] # select columns
    temp5[] <- lapply(temp5, as.character)
    temp5 <- temp5[,-grep("Group", names(temp5))]
    temp6 <- reshape2::melt(temp5,id.vars="Name")
    colnames(temp6) <- c("Pops","Loci","Value")

    temp6$Loci <- as.character(temp6$Loci)

    Output <- tidyr::unnest(temp6%>%group_by(Pops,Loci)%>%do(split=frequencies(.))%>%
                             ungroup()%>%data.frame(.,stringsAsFactors=F))
  }

  #If the output is to be in wide format or to be left (default) as long (transposed format). Note taht wide formate is accepted for genotype simulation by alleleotype_genepop().
  if(wide & fullpanel){print("The full allele frequency panel cannot be reported in wide format. Data will be returned in long format.")}


  if(wide & !fullpanel){
    Output[] <- lapply(Output,as.character)
    Output <- Output %>% data.table::dcast(.,Population~Loci, value.var="MAF") %>% data.frame()
    #data.table::dcast adds an "x" to the column names which must be removed for matching
    names(Output)[2:length(names(Output))]=gsub("X","",names(Output)[2:length(names(Output))])
    Output <- Output[,c("Population",names(temp2)[-length(temp2)])]
    colnames(Output)[1]="Pop"
  }

  return(Output)# return the allele frequency table

}
