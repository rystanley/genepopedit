# Simulate genotypes from pooled DNA samples
#' @title Alleleotype -> GENEPOP format
#' @description Simulate individuals using pooled DNA major allele frequencies. Return as GENEPOP formatted data.
#' @param input path to input file (csv) containing major allele frequencies. First column is the SNP names and the remaining columsn are the population based allele frequencies from the pooled DNA sample.
#' @param numsim number of simulated individuals per population to be returned.
#' @param path the filepath and filename of output.
#' @import magrittr
#' @importFrom plyr round_any
#' @importFrom dplyr group_by do ungroup select mutate_each mutate
#' @importFrom data.table fread melt dcast as.data.table
#' @rdname alleleotype_genepop
#' @export

##
alleleotype_genepop <- function(input,numsim=100,path){

  ## Functions used for simulating individuals

    #create dataframes of alleles which will be used to create genotypes
      simind <- function(x,n=100){
        coeff <- n/100 #multiple of 100
        return(data.frame(alleleotype=c(rep("001",round(x*100)*coeff*2),
                                        rep("003",(round((1-x)*100)*coeff*2)))))}

    #Function used to shuffle output from 'simind' for each column (locus) to randomly simulate individual genotypes
      shufflefun <- function(x){sample(x,length(x),replace=FALSE)}

    ## The number of individuals to be sampled will be simulated at base 100 (for rounding)
  nsim <- plyr::round_any(numsim,100,f=ceiling)

  # Simulations are based on alleles so nsim and numsim need to be multiplied by two to generate snps
  numsim <- numsim*2

  #if inputdata is a path read in the data
    if(is.character(input)){
      df <- data.table::fread(input,stringsAsFactors = FALSE)
    }else{df <- data.table::as.data.table(input)}

  #create molten data for dplyr grouping
    df <- data.table::melt(df,id.vars="Pop")

  #create populate the alleles from alleleotype data prune to specified sample size (numsim)
    df2 <- df%>%group_by(Pop,variable)%>%
      do(simind(.$value,n=nsim))%>%
      sample_n(.,numsim)%>%ungroup()%>%
        as.data.table(stringsAsFactors=FALSE)

  #Simulate individuals by sampling without replacement the possible alleles. This is essentially shuffling the potential alleles randomly at each locus and then drawing individuals from this

    df2$newid="999"
    #set up a idvariable which helps dcast deal with repeat values see link [[1]]
    df2 <- df2%>%group_by(Pop,variable)%>%
      mutate(newid = paste(Pop,seq_along(newid)))%>%
      ungroup()%>%as.data.table(stringsAsFactors=FALSE)

    df2[] <- lapply(df2, as.character) # covert from factor to character

    #go from wide to long using the new ID variable which will be removed then the order of each column will be shifted
      df3 <- df2%>%dcast(.,Pop+newid~variable, value.var="alleleotype")%>%select(.,-newid)%>%
        group_by(Pop)%>%mutate_each(funs(shufflefun))%>%ungroup()%>%data.frame()

    ## split individuals up into single SNP calls
      firstAllele.row <- df3[!1:nrow(df3) %% 2 == 0,-1] #odd rows
      secondAllele.row <- df3[1:nrow(df3) %% 2 == 0,-1] #even ros

  #Create temporary matrix which will be populated with interlaced even and odd columns
    #matrix dimensions
      rows.combined <- dim(firstAllele.row)[1]
      cols.combined <- dim(firstAllele.row)[2]*2 #because alleles for each SNP are seperated

    #create NULL frame to be populated
      AlleleFrame <- as.data.frame(matrix(NA, nrow=rows.combined, ncol=cols.combined))

    #populate the frame
      AlleleFrame[,seq(1,cols.combined,2)] <- firstAllele.row #first allele
      AlleleFrame[,seq(2,cols.combined,2)] <- secondAllele.row #second allele

    #Get the indicies for first and second alleles of each SNP (consecutive pairs)
      firstAllele.col <- which(!1:cols.combined %% 2 == 0) #odd columns
      secondAllele.col <- which(1:cols.combined %% 2 == 0) #even columns

    #Create matrix which will be used in apply to specify which columns to be combined
      ColumnPasteFrame=rbind(firstAllele.col,secondAllele.col)

    #Combine the alleles for each SNP
      CombinationFrame <- apply(ColumnPasteFrame,2,FUN=function(x){apply(AlleleFrame[,x],1,paste,collapse="")})

    #Add the population labels
      AlleleFrame <- data.table::as.data.table(cbind(as.character(df3[!1:nrow(df3) %% 2 == 0,"Pop"]),
                                   CombinationFrame),
                                    stringsAsFactors=FALSE )#Add population
   #Fix column names
      colnames(AlleleFrame)=c("Pop",colnames(firstAllele.row))

## Reformat into GENEPOP format

    #Temporary column for group-mutate
      AlleleFrame$ID <- AlleleFrame$Pop

    #mutate 'Pop' column to add sequential unique sample IDs within each population drop temporary group variable
      AlleleFrame <- AlleleFrame%>%group_by(ID)%>%
        mutate(Pop=paste0(Pop,"_",1:length(Pop)))%>%
        ungroup()%>%select(-ID)%>%data.frame()

# Paste together the Loci as one long integer separated for each loci by a space
      Loci <- do.call(paste,c(AlleleFrame[,-1], sep=" "))

#Grab the Population tags that each invididual had following the format ID_,__
      PopVec <- paste0(AlleleFrame[,"Pop"]," ,  ")

#Find population names and insert Genepop "Pop" delimiters
      NameExtract <- substr(PopVec,1,regexpr("_",PopVec)-1)

      if(length(table(NameExtract))>1){PopLengths <- table(NameExtract)[-length(table(NameExtract))]} else {PopLengths=1}

      if(length(table(NameExtract))==2){PopPosition = PopLengths+1}

      if(length(table(NameExtract))>2){
        PopPosition <- c(PopLengths[1]+1,rep(NA,(length(PopLengths)-1)))
        for (i in 2:length(PopLengths)){
          PopPosition[i] <- PopLengths[i]+PopPosition[i-1]
        }
      }

#Paste these to the Loci
      Loci <- paste(PopVec,Loci,sep="")

#Insert the value of "Pop" which partitions the data among populations #only if more than one population
      if(length(table(NameExtract))!=1){Loci <- insert_vals(Vec=Loci,breaks=PopPosition,newVal="Pop")}

#Compile the Genepop format
      Output <- c("Alleleotype-GENEPOP by genepopedit",colnames(AlleleFrame)[-1],
                  "Pop",Loci)

# Save the file
      write.table(Output,path,col.names=FALSE,row.names=FALSE,quote=FALSE)

}
