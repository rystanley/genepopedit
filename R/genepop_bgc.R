# Genepop -> Bayesian Genomic Clines (BGC)
#' @title Convert Genepop to Bayesian Genomic Clines (BGC) format.
#' @description Function to convert Genepop to BGC format.
#' @param GenePop the genepop data to be manipulated. This can be either a file path
#' or a dataframe read in with tab seperation, header=FALSE , quote="", and stringsAsFactors=FALSE.
#' This will the standard genepop format with a the first n+1 rows corresponding the the n loci names,
#' or a single commma deliminated row of loci names followed by the locus data. Populations are
#' seperated by "Pop". Each individual ID is linked to the locus data by " ,  " (space,space space) and is read in as
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
#' @importFrom tidyr separate
#' @importFrom dplyr filter do group_by ungroup
#' @export


genepop_bgc <- function(GenePop,popdef,fname,path){

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
  ColumnData <- gsub("\r","",ColumnData)#remove any hidden carriage returns
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
  if (length(temp2)==length(ColumnData)){colnames(temp2) <- ColumnData}
  if (length(temp2)!=length(ColumnData)){stacks.version="No STACKS version specified"}

  ## Get the population names (prior to the _ in the Sample ID)
  NamePops <- temp[,1] # Sample names of each
  NamePops <- gsub(" ","",NamePops) #get rid of space
  NameExtract <- substr(NamePops,1,regexpr("_",NamePops)-1)# extract values from before the "_" to denote populations

  #convert the snp data into character format to get rid of factor levels
  temp2[] <- lapply(temp2, as.character)

  alleleEx <- as.character(temp2[1,1]) #example allele

  #get the allele values summary header
  firstAllele <-  as.data.frame(sapply(temp2,function(x)as.numeric(as.character(substring(x,1,(nchar(alleleEx)/2))))))
  secondAllele <-  as.data.frame(sapply(temp2,function(x)as.numeric(as.character(substring(x,(nchar(alleleEx)/2)+1,nchar(alleleEx))))))

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

  if(is.character(popdef)){popdef <- read.csv("popdef.csv",header=T)} #if popdef is a path then read it in

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
        Admixed_BGC <- NULL

        for(i in snpnames){
          # grab vector of alleles and delete replace missing values (-9) with NA
          P1_alleles <- P1_raw[,i];P1_alleles[which(P1_alleles==-9)]=NA
          P2_alleles <- P2_raw[,i];P2_alleles[which(P2_alleles==-9)]=NA

          #If the population only has one allele for a given locus then a zero and the allele have be be added
          if(length(table(P1_alleles))==1){
            hold <- as.data.frame(table(P1_alleles))
            hold[,1] <- as.character(hold[,1])
            hold <- rbind(hold,c(setdiff(as.numeric(Allele_Map[which(Allele_Map$SNP==i),c("Allele1","Allele2")]),hold[1,1]),0)) #add in the extra value
            hold <- hold[order(hold[,1]),] #sort the right order from a conventional table output
            P1_alleles <- hold[,2]
            rm(hold)
          } else {P1_alleles <- as.character(as.data.frame(table(P1_alleles))[,2])}

          if(length(table(P2_alleles))==1){
            hold <- as.data.frame(table(P2_alleles))
            hold[,1] <- as.character(hold[,1])
            hold <- rbind(hold,c(setdiff(as.numeric(Allele_Map[which(Allele_Map$SNP==i),c("Allele1","Allele2")]),hold[1,1]),0)) #add in the extra value
            hold <- hold[order(hold[,1]),] #sort the right order from a conventional table output
            P2_alleles <- hold[,2]
            rm(hold)
          } else {P2_alleles <- as.character(as.data.frame(table(P2_alleles))[,2])}


          #for a given locus get the format for BGC
          P1_temp <- c(paste("locus_",i,sep=""),paste(P1_alleles[1],P1_alleles[2],sep=" "))
          P2_temp <- c(paste("locus_",i,sep=""),paste(P2_alleles[1],P2_alleles[2],sep=" "))

          #Combine output sequentially for each locus
          P1_BGC <- c(P1_BGC,P1_temp)
          P2_BGC <- c(P2_BGC,P2_temp)
        }


#Convert the admixed data to BGC format --------------

        MixedStruct=P3_raw

       #Data wrangle and construct the format for admixture populatons as per the instructions and table 2 in the BGC manual
          MixedData <- NULL
          for(col in names(MixedStruct)[4:length(MixedStruct)])
            {

            temp3 <- MixedStruct[,c("ID","Pop",col)]
            Locushold <- paste("locus_",col,sep="") # Start the locus data

              for (i in unique(MixedStruct$Pop)) # each populatoin
                {

                temp4 <- dplyr::filter(temp3,Pop==i) # subset for the population

                #define major and minor alleles(remove missing data (-9))

                # all missing data
                if(length(table(temp4[,3]))==1){if(unique(temp4[,3])==-9){
                  writeLines(paste0("Warning no data available for loci(",col,
                                    ") and population ",i,
                                    " coded as -9 -9 for all individuals"))}}
                if(length(table(temp4[,3]))==1){if(unique(temp4[,3])==-9){
                  temp4[,3]=999}}

                temp4a <- temp4
                temp4a[which(temp4[,3]==-9),3] <- NA

                AlleleMajor <- as.numeric(names(which(table(temp4a[,3])==max(table(temp4a[,3])))))
                AlleleMinor <- as.numeric(names(which(table(temp4a[,3])==min(table(temp4a[,3])))))

                #Reformat the data for one row for each individaul (ID, Pop, Allele1, Allele2)
                temp5 <- data.frame(ID=temp4[seq(1,nrow(temp4),2),"ID"],
                                    Pop=temp4[seq(1,nrow(temp4),2),"Pop"],
                                    allele1=temp4[seq(1,nrow(temp4),2),col],
                                    allele2=temp4[seq(2,nrow(temp4),2),col],
                                    alleleMajor=AlleleMajor,
                                    alleleMinor=AlleleMinor)

                #code for homozygous minor homozygous major and missing data.
                temp6 <- as.data.frame(temp5%>%group_by(ID)%>%do(col1=aCount(.)[1],col2=aCount(.)[2])%>%ungroup())
                temp6 <- temp6[ordered(temp5$ID),] #back to the original order

                temp7 <- paste(temp6[,"col1"],temp6[,"col2"],sep=" ")

                Locushold <- c(Locushold,paste("pop_",i,sep=""),temp7)

              } #end of population loop

              MixedData <- c(MixedData,Locushold) # add each successive locus
            } #end of locus loop


##Save output for BGC formated for the parental and mixed populations ------------
      if(substring(path,nchar(path))!="/"){path=paste0(path,"/")}

      write.table(x = P1_BGC,file=paste0(path,fname,"_Parental1_BGC.txt",sep=""),
                  sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

      write.table(x = P2_BGC,file=paste0(path,fname,"_Parental2_BGC.txt",sep=""),
                  sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

      write.table(x = MixedData,file=paste(path,fname,"_Admixed_BGC.txt",sep=""),
                  sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

} #end of function
