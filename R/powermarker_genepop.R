# Powermarker Genepop
#' @title Convert Powermarker format to Genepop
#' @description Function converts Powermarker to standard Genepop format
#' @param Powermarker the Powermarker data to be manipulated. The first column should be the 'POP_ID' column which identifies the populations.
#' The second column should be 'Sample_ID' which designates the individual sample IDs. The remaining columns contain the locus alleles in the format of 'A/A' etc.
#' @param missing_data The symbol (typically "-", "?", or "9") which will be replaced with 000.
#' @param path the filepath and filename of output.
#' @param sampleid Whether you want the sampleid in your Genepop to be based off the 'Sample_ID' column (TRUE) or the 'POP_ID' column (FALSE) in the Powermarker file
#' @rdname powermarker_genepop
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @export
#'

powermarker_genepop<-function(Powermarker, missing_data, path,sampleid=TRUE){

  classdef <- utils::read.table(Powermarker, header = TRUE, nrows = 5) # find column classes
  classes <- sapply(classdef, class)
  input <- utils::read.table(Powermarker, header = TRUE, colClasses = classes,stringsAsFactors = FALSE)

   #input<-utils::read.table(Powermarker)
  input.hold<-input[c(1,2)]
  input<-input[,-c(1,2)]
  loci_names<-colnames(input)


input.conv <- as.data.frame(sapply(input, gsub, pattern = missing_data, replacement="000"), stringsAsFactors=FALSE)
input.conv <- as.data.frame(sapply(input.conv, gsub, pattern = "A", replacement="001", fixed = TRUE), stringsAsFactors=FALSE)
input.conv <- as.data.frame(sapply(input.conv, gsub, pattern = "C", replacement="002", fixed = TRUE), stringsAsFactors=FALSE)
input.conv <- as.data.frame(sapply(input.conv, gsub, pattern = "G", replacement="003", fixed = TRUE), stringsAsFactors=FALSE)
input.conv <- as.data.frame(sapply(input.conv, gsub, pattern = "T", replacement="004", fixed = TRUE), stringsAsFactors=FALSE)
input.conv <- as.data.frame(sapply(input.conv, gsub, pattern = "/", replacement="", fixed = TRUE), stringsAsFactors=FALSE)


#the number of individuals for all populations but the last (Pop tagged to the end)
if(length(table(input.hold[,1]))>1){PopLengths <- table(factor(input.hold[,1], levels=unique(input.hold[,1])))[-length(table(input.hold[,1]))]} else {PopLengths=1}

if(length(table(input.hold[,1]))==2){PopPosition = PopLengths+1}

if(length(table(input.hold[,1]))>2){
  PopPosition <- c(PopLengths[1]+1,rep(NA,(length(PopLengths)-1)))
  for (i in 2:length(PopLengths)){
    PopPosition[i] <- PopLengths[i]+PopPosition[i-1]
  }
}

PopColLengths <- table(factor(input.hold[,1], levels=unique(input.hold[,1])))


popvector=rep(1:length(PopColLengths),times=PopColLengths) #vector to differentiate the pops based on the locaton of the "Pop" labels.


# paste together the Loci as one long integer separated for each loci by a space
Loci <- do.call(paste,c(input.conv[,], sep=" "))

#Grab the Population tags that each invididual had following the format ID_,__
NamePops <- as.character(input.hold[,2])
NamePops <- gsub("_","-",NamePops) # to ensure sampleIDs match genepopedit
SampleID <- NamePops

for(i in unique(popvector)){
  commonname <- common_string(SampleID[which(popvector==i)])
  NamePops[which(popvector==i)] <- paste0(commonname,"_")
  SampleID[which(popvector==i)] <- gsub(commonname,paste0(commonname,"_"),SampleID[which(popvector==i)])
}

#Get rid of any extra "-"
SampleID <- gsub("-_","_",SampleID)

if(sampleid){PopVec <- paste0(SampleID," ,  ")}
if(!sampleid){
  SampleNumbers <- NULL
  for(i in unique(input.hold[,1])){
    numclass <- 1:length(which(input.hold[,1]==i))
    SampleNumbers <- c(SampleNumbers,numclass)
  }

  PopVec <- paste0(input.hold[,1],"_",SampleNumbers," ,  ")

  }

#Paste these to the Loci
Loci <- paste(PopVec,Loci,sep="")

#Insert the value of "Pop" which partitions the data among populations #only if more than one population
if(length(table(input.hold[,1]))!=1){Loci <- insert_vals(Vec=Loci,breaks=PopPosition,newVal="Pop")}

#Add the first "Pop" label
Loci <- c("Pop",Loci)

Output<-c("Powermarker to Genepop by genepopedit",loci_names,Loci)

#Save the Genepop file
  utils::write.table(Output,path, col.names=FALSE, row.names=FALSE, quote=FALSE)
}

