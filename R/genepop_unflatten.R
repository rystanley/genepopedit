# GenePop unflatten
#' @title Convert to Genepop format from a flattened dataframe.
#' @description Function returns Genepop file meta-data.
#' @param df dataframe with the first column holding sampleIDs (e.g. BON_01) and the remaining columns holding loci. Column names of loci will be used as loci names in the genepop output.
#' df must be an object in the workspace.
#' @param path the filepath and filename of output.
#' @rdname genepop_unflatten
#' @importFrom utils write.table
#' @export
#'
genepop_unflatten <- function(df,path){

  #Make sure all loci are characters and now factors
  df <- as.data.frame(apply(df,2, as.character),stringsAsFactors = F)

  #check to see if the dataframe was created using genepop_flatten and remove the extra columns
  if("Population" %in% colnames(df)){df <- df[,-grep("Population",colnames(df))]}
  if("SampleNum" %in% colnames(df)){df <- df[,-grep("SampleNum",colnames(df))]}

  NamePops <- df[,1] # Sample names of each
  NameExtract <- substr(NamePops,1,regexpr("_",NamePops)-1)

  #Loci
  temp2 <- df[,2:length(df)]

  #Order the vector of population names and the loci by population. This is done in situations where the populations labels might be
  #missaligned (e.g. Pop1 Pop1 .... Pop2 ... Pop3 ... Pop1Pop4 ....) #Note that this will change the order to alpha numeric. This should
  #not have any influence on the analysis

  reorder <- order(NameExtract) #vector which defines sorting

  NameExtract <- NameExtract[reorder]
  temp2 <- temp2[reorder,]
  NamePops <- NamePops[reorder]

  ## Now add the population tags using npops (number of populations and Pops for the inter differences)

  PopLengths <- table(factor(NameExtract, levels=unique(NameExtract)))[-length(table(NameExtract))]

  if(length(table(NameExtract))==2){PopPosition = PopLengths+1}

  if(length(table(NameExtract))>2){
    PopPosition <- c(PopLengths[1]+1,rep(NA,(length(PopLengths)-1)))
    for (i in 2:length(PopLengths)){
      PopPosition[i] <- PopLengths[i]+PopPosition[i-1]
    }
  }

  #Combine loci together add sampleIDs and 'Pop' labels.
  Loci <- do.call(paste,c(temp2[,], sep=" "))
  Loci <- paste0(NamePops," ,  ",Loci)
  if(length(table(NameExtract))!=1){Loci <- insert_vals(Vec=Loci,breaks=PopPosition,newVal="Pop")}
  Loci <- c("Pop",Loci)

  Output <- c("No STACKS version specified",names(temp2),Loci)

  # Save the file
  utils::write.table(Output,path,col.names=FALSE,row.names=FALSE,quote=FALSE)

}

