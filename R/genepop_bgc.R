# Genepop -> Bayesian Genomic Clines (BGC)
#' @title Convert Genepop to Bayesian Genomic Clines (BGC) format.
#' @description Function to convert Genepop to BGC format.
#' @param genepop the genepop data to be manipulated. This can be either a file path
#' or a dataframe read in with tab separation, header=FALSE , quote="", and stringsAsFactors=FALSE.
#' This will the standard genepop format with a the first n+1 rows corresponding the the n loci names,
#' or a single comma deliminated row of loci names followed by the locus data. Populations are
#' separated by "Pop". Each individual ID is linked to the locus data by " ,  " (space,space space) and is read in as
#' as a single row (character).
#' @param popdef is a dataframe or path to a csv.
#' This dataframe contains two columns. Column 1 corresponds to the population names. These names
#' should match the individual IDs (e.g. BON_01 ,  110110 would be 'BON'). The next column
#' has the grouping classification corresponding to each population
#' defining parental 1 ("P1") parental 2 ("P2") and admixed ("Admixed") populations.
#' Note the classifications must be exactly as specified (caps sensitive). If populations are omitted from
#' this dataframe then they will be omitted from the output files.
#' @param fname collective name assigned to each of the output files for BGC.
#' e.g. "Lobster_analysis" would result in
#' "Lobster_analysis_P1.txt","Lobster_analysis_P2.txt", and "Lobster_analysis_Admixed.txt"
#' @param path file path to directory where the BGC files (3) will be saved.
#' @rdname genepop_bgc
#' @import magrittr
#' @import dplyr
#' @importFrom data.table fread
#' @importFrom utils write.table
#' @importFrom utils read.csv
#' @export

genepop_bgc <- function(genepop,popdef,fname,path){

  #Check to see if genepop is a data.frame from the workspace and convert to data.table
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
    genepop <- as.data.table(genepop,stringsAsFactors = FALSE)
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

  #convert the snp data into character format to get rid of factor levels
  temp2[] <- lapply(temp2, as.character)

  #allele coding length
  alleleEx <- max(sapply(temp2[,1],FUN=function(x){nchar(as.character(x[!is.na(x)]))})) #presumed allele length

  #check to make sure the allele length is a even number
  if(!alleleEx %% 2 ==0){stop(paste("The length of each allele is assumed to be equal (e.g. loci - 001001 with 001 for each allele), but a max loci length of", alleleEx, "was detected. Please check data."))}

  #get the allele values summary header
  firstAllele <-  as.data.frame(sapply(temp2,function(x)as.numeric(as.character(substring(x,1,(alleleEx/2))))))
  secondAllele <-  as.data.frame(sapply(temp2,function(x)as.numeric(as.character(substring(x,(alleleEx/2)+1,alleleEx)))))

  # switch from combined allele in one row to two allele each in their own row according to a given locus
  holdframe <- rbind(firstAllele,secondAllele)#create a dummy data frame
  holdframe[!1:nrow(holdframe) %% 2 == 0,] = firstAllele #odd rows
  holdframe[1:nrow(holdframe) %% 2 == 0,] = secondAllele #even rows

  holdframe[holdframe==0]= -9 # replace missing values with -9

  groupvec <- NameExtract
  for (i in 1:length(unique(NameExtract))) # replace with numbers
  {
    groupvec[which(groupvec==unique(NameExtract)[i])] = i
  }

  holdframe=cbind(rep(NamePops,each=2),rep(groupvec,each=2),rep(NameExtract,each=2),holdframe)
  colnames(holdframe)[1:3]=c("ID","PopID","Pop")

  #single digit format as required by BGC (1-4)
  holdframe[holdframe==100]=1
  holdframe[holdframe==110]=2
  holdframe[holdframe==120]=3
  holdframe[holdframe==130]=4

  if(is.character(popdef)){popdef <- utils::read.csv("popdef.csv",header=T)} #if popdef is a path then read it in

  #Extract the parental data and admixed data
  P1_raw <- holdframe[which(holdframe$Pop %in% popdef[which(popdef[,2]=="P1"),1]),]#Parental 1
  P2_raw <- holdframe[which(holdframe$Pop %in% popdef[which(popdef[,2]=="P2"),1]),]#Parental 2
  P3_raw <- holdframe[which(holdframe$Pop %in% popdef[which(popdef[,2]=="Admixed"),1]),]#Admixed

  # names of the snps
  snpnames <- colnames(temp2)

  #map to be used for missing alleles (if any present across loci)
  Allele_Map <- data.frame(SNP=snpnames,
                           Allele1=rep(999,length(snpnames)),
                           Allele2=rep(999,length(snpnames))) # 999 is a dummy placeholder

  for(i in 1:length(snpnames)){
    #unique alleles for a given snp (locus)
    alleleVals <- as.data.frame(table(as.character(c(P1_raw[,snpnames[i]],P2_raw[,snpnames[i]],P3_raw[,snpnames[i]]))))

    # if there is missing data (-9) delete it as a possibe allele
    if(length(which(alleleVals[,1]==(-9)))>0){
      alleleVals <- alleleVals[-which(alleleVals[,1]==(-9)),]
    }

    Allele_Map[i,"Allele1"]=as.character(alleleVals[1,1])
    Allele_Map[i,"Allele2"]=as.character(alleleVals[2,1])
  }

  #NULL vectors
  P1_BGC <- NULL
  P2_BGC <- NULL

  for(i in snpnames){
    # grab vector of alleles and delete replace missing values (-9) with NA
    P1_alleles <- P1_raw[,i];P1_alleles[which(P1_alleles==-9)]=NA
    P2_alleles <- P2_raw[,i];P2_alleles[which(P2_alleles==-9)]=NA

    #If the population only has one allele for a given locus then a zero and the allele have be be added
    if(length(table(P1_alleles))==1|sum(is.na(P1_alleles))==length(P1_alleles)){
      if(length(table(P1_alleles))==1){
      hold <- as.data.frame(table(P1_alleles))
      hold[,1] <- as.character(hold[,1])
      hold <- rbind(hold,c(setdiff(as.numeric(Allele_Map[which(Allele_Map$SNP==i),c("Allele1","Allele2")]),hold[1,1]),0)) #add in the extra value
      hold <- hold[order(hold[,1]),] #sort the right order from a conventional table output
      P1_alleles <- hold[,2]
      rm(hold)} else {P1_alleles <- c(0,0)}
    } else {P1_alleles <- as.character(as.data.frame(table(P1_alleles))[,2])}

    if(length(table(P2_alleles))==1 | sum(is.na(P2_alleles))==length(P2_alleles)){
      if(length(table(P2_alleles))==1){
      hold <- as.data.frame(table(P2_alleles))
      hold[,1] <- as.character(hold[,1])
      hold <- rbind(hold,c(setdiff(as.numeric(Allele_Map[which(Allele_Map$SNP==i),c("Allele1","Allele2")]),hold[1,1]),0)) #add in the extra value
      hold <- hold[order(hold[,1]),] #sort the right order from a conventional table output
      P2_alleles <- hold[,2]
      rm(hold)} else{P2_alleles <- c(0,0)}
    } else {P2_alleles <- as.character(as.data.frame(table(P2_alleles))[,2])}


    #for a given locus get the format for BGC
    P1_temp <- c(paste("locus_",i,sep=""),paste(P1_alleles[1],P1_alleles[2],sep=" "))
    P2_temp <- c(paste("locus_",i,sep=""),paste(P2_alleles[1],P2_alleles[2],sep=" "))

    #Combine output sequentially for each locus
    P1_BGC <- c(P1_BGC,P1_temp)
    P2_BGC <- c(P2_BGC,P2_temp)
  }

  ##Save output for BGC formatted for the parental populations ------------
  if(substring(path,nchar(path))!="/"){path=paste0(path,"/")}

  utils::write.table(x = P1_BGC,file=paste0(path,fname,"_Parental1_BGC.txt",sep=""),
              sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

  utils::write.table(x = P2_BGC,file=paste0(path,fname,"_Parental2_BGC.txt",sep=""),
              sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

  #Convert the admixed data to BGC format --------------

  #subset data for admixed populations
  missingfix<- function(x){ #create functions for apply loop
    hold=x
    hold[grep("000",hold)]=NA
    return(hold)}

  #Remove Alleles with missing data and replace with NA
  temp3 <- apply(temp2,2,missingfix)

  #convert to zygosity format (2 0 - homozygous major, 0 2 - homozygous minor, 1 1 - heterozygous, -9 -9 - missing data )
  temp4 <- apply(temp3,2,FUN = function(x) {majorminor(x,allele_length=alleleEx)})

  MixedStruct <- temp4[which(NameExtract %in% popdef[which(popdef[,2]=="Admixed"),1]),]
  MixedPops <- NameExtract[which(NameExtract %in% popdef[which(popdef[,2]=="Admixed"),1])]

  #the number of individuals for all populations but the last (Pop tagged to the end)
  PopLengths <- table(MixedPops)[-length(table(MixedPops))]

  if(length(table(MixedPops))==2){PopPosition = PopLengths+1}

  if(length(table(MixedPops))>2){
    PopPosition <- c(PopLengths[1]+1,rep(NA,(length(PopLengths)-1)))
    for (i in 2:length(PopLengths)){
      PopPosition[i] <- PopLengths[i]+PopPosition[i-1]
    }
  }

  #Insert the population labels
  if(length(table(MixedPops))!=1){
  temp5 <- apply(MixedStruct,2,function(x){insert_vals(x,breaks=PopPosition,
    newVal=paste0("pop_",unique(MixedPops)[2:length(unique(MixedPops))]))})} else {
    temp5 <- MixedStruct}

  temp5=as.data.frame(temp5,stringsAsFactors = FALSE)

  #Add the "locus_" and first "pop_" labels
  temp6=as.matrix(rbind(paste0("locus_",colnames(temp5)),
              rep(paste0("pop_",unique(MixedPops)[1]),length(temp5)),
              temp5))

  #redim as a single vector
  MixedData=as.vector(temp6)

  ##Save output for BGC formatted for the parental and mixed populations ------------
  utils::write.table(x = MixedData,file=paste(path,fname,"_Admixed_BGC.txt",sep=""),
              sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

} #end of function
