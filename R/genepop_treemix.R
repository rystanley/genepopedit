# Genepop_Treemix
#' @title Convert a Genepop to input required for phython processing to TREEMIX.
#' @description Convert Genepop to within-clusters binary PED file for TREEMIX.
#' @param genepop the genepop data to be manipulated. This can be either a file path
#' or a dataframe read in with tab separation, header=FALSE , quote="", and stringsAsFactors=FALSE.
#' This will be the standard genepop format with the first n+1 rows corresponding to the n loci names,
#' or a single comma delimited row of loci names followed by the locus data. Populations are
#' separated by "Pop". Each individual ID is linked to the locus data by " ,  " (space,space space) and is read in as
#' as a single row (character).
#' @param where.pgdspider A file path to the PGDspider installation folder.
#' @param where.plink A file path to the PLINK installation folder.
#' @param allocate.pgd.ram An integer value in GB to specify the maximum amount of RAM to allocate to PGDspider. The default is 1 GB, which should be sufficient for most analyses.
#' @param keep_inter A logical condition statement (default : FALSE) specifying whether to keep the map, ped, and clustering files generated during the conversion.
#' @param path file path to directory where the gzipped Treemix input file will be saved.
#' @rdname genepop_treemix
#' @importFrom R.utils gzip
#' @importFrom utils write.table
#' @importFrom utils read.table
#' @export

genepop_treemix<-function(genepop,where.pgdspider,where.plink,allocate.pgd.ram=1,keep_inter=FALSE,path){

  path.start<-path

  writeLines("Note that this function works with Plink 1.9 for now. If you have Plink 2.0, the function may not perform as intended.")

  ## File path checks
  if(substring(path,first = nchar(path),last=nchar(path))!="/"){
    path <- paste0(path,"/")
  }
  path.start.PGD <- gsub(x = path.start, pattern = " ", replacement = "\\")

  where.pgdspider.PGD <- gsub(x = where.pgdspider, pattern = " ",
                              replacement = "\\ ", fixed = TRUE)

  #Ram warning messages.
  if(allocate.pgd.ram%%1 != 0){
    stop("Please specify an integer GB value to allocate to PGDspider.")
  }

  if (Sys.info()["sysname"] == "Windows" & allocate.pgd.ram >
      1024) {
    allocate.pgd.ram = 1024
    writeLines("Note that currently PGDspider can only utilize ~1 GB of ram on windows based operating systems. Periodically check back to https://github.com/rystanley/genepopedit for any updates to this limitation.\n               ")
  }

  #Ram allocation - gb to mb
  allocate.pgd.ram <- allocate.pgd.ram * 1024

  writeLines("Converting GENEPOP to PED format for PLINK. Note that this step will format data needed to execute the Phyton conversions. Refer to TREEMIX documentation for Phython conversion workflow and scripting")
  writeLines("\n\n             ")
  writeLines("Warning messages are expected as part of conversion process using PGDspider. Note this process can be inhibited by inadequate system resources.\n\n             ")

#Create the spid files required by PGD for conversion
GP_PED_SPID_Top<-"# spid-file generated: Thu May 19 13:29:37 ADT 2016\n\n # GENEPOP Parser questions\n PARSER_FORMAT=GENEPOP
  # Enter the size of the repeated motif (same for all loci: one number; different: comma separated list (e.g.: 2,2,3,2):
  GENEPOP_PARSER_REPEAT_SIZE_QUESTION=
  # Select the type of the data:
    GENEPOP_PARSER_DATA_TYPE_QUESTION=SNP
  # How are Microsat alleles coded?
  GENEPOP_PARSER_MICROSAT_CODING_QUESTION=REPEATS

  # PED Writer questions
  WRITER_FORMAT=PED

  # Save MAP file"

map.loc<-paste0("PED_WRITER_MAP_FILE_QUESTION= ", "PGDtest")
GP_PED_SPID_Bottom<-"# Replacement character for allele encoded as 0 (0 encodes for missing data in PED):
  PED_WRITER_ZERO_CHAR_QUESTION=
  # Specify the locus/locus combination you want to write to the PED file:
    PED_WRITER_LOCUS_COMBINATION_QUESTION=
  # Do you want to save an additional MAP file with loci information?
    PED_WRITER_MAP_QUESTION=true"
  spid.file <- c(GP_PED_SPID_Top, map.loc, GP_PED_SPID_Bottom)
  write(x = spid.file, file = paste0(where.pgdspider.PGD, "/", "GP_PED.spid"))


  file.copy(from = genepop, to = where.pgdspider, overwrite = TRUE)
  GenePop.name <- stringr::str_split(string = genepop, pattern = "/")
  GenePop.name <- unlist(GenePop.name)
  GenePop.name <- GenePop.name[grep(x = GenePop.name, pattern = ".txt")]
  file.rename(from = paste0(where.pgdspider, GenePop.name),
              to = paste0(where.pgdspider, "GPD_for_PED_to_BED.txt"))
  remember.TOPLOC.path <- paste0(where.pgdspider, "GPD_for_PED_to_BED.txt")

  if (Sys.info()["sysname"] != "Windows") {
    input.file.call <- "-inputfile GPD_for_PED_to_BED.txt"
    execute.SPIDER <- paste0("java -Xmx", allocate.pgd.ram,
                             "m -Xms512m -jar PGDSpider2-cli.jar")
    spid.call <- "-spid GP_PED.spid"
    input.format <- "-inputformat GENEPOP"
    output.format <- "-outputformat PED"
    goto.spider <- paste0("cd ", where.pgdspider.PGD, "; ",
                          execute.SPIDER)
    output.file.path <- "-outputfile PGDtest.ped"
    run.PGDspider <- paste0(goto.spider, " ", input.file.call,
                            " ", input.format, " ", output.file.path, " ", output.format,
                            " ", spid.call)
    system(run.PGDspider)
  }

  if (Sys.info()["sysname"] == "Windows") {
    input.file.call <- "-inputfile GPD_for_PED_to_BED.txt"
    execute.SPIDER <- paste0("java -Xmx", allocate.pgd.ram,
                             "m -Xms512m -jar PGDSpider2-cli.jar")
    spid.call <- "-spid GP_PED.spid"
    input.format <- "-inputformat GENEPOP"
    output.format <- "-outputformat PED"
    goto.spider <- paste0("cd ", where.pgdspider.PGD, " && ",
                          execute.SPIDER)
    output.file.path <- "-outputfile PGDtest.ped"
    run.PGDspider <- paste0(goto.spider, " ", input.file.call,
                            " ", input.format, " ", output.file.path, " ", output.format,
                            " ", spid.call)
    shell(run.PGDspider)
  }

  ped.path <- paste0(where.pgdspider, "/", "PGDtest.ped")
  map.path <- paste0(where.pgdspider, "/", "PGDtest.map")
  file.copy(from = ped.path, to = where.plink, overwrite = TRUE)
  file.copy(from = map.path, to = where.plink, overwrite = TRUE)
  plink_ped_path <- paste0(where.plink, "/", "PGDtest.ped")
  plink_map_path <- paste0(where.plink, "/", "PGDtest.map")
  writeLines("\n\n             ")
  callplink<- paste0("cd ", where.plink)

  if (Sys.info()["sysname"] != "Windows") {
    plink.input <- paste0(callplink, "; ", "./plink --noweb --file PGDtest --make-bed --out BinaryPED")
    system(plink.input)
  }

  if (Sys.info()["sysname"] == "Windows") {
    plink.input<-paste0(callplink, " && ", "plink --noweb --file PGDtest --make-bed --out BinaryPED")
    shell(plink.input)
  }
  anal.name=stringr::str_split(string = "PGDtest", pattern = "/")
  anal.name=unlist(anal.name)
  anal.name=anal.name[length(anal.name)]
  map.name=paste0(anal.name,".map")
  ped.name=paste0(anal.name,".ped")
  remember.ped.plink<-paste0(where.plink,ped.name)
  remember.map.plink<-paste0(where.plink,map.name)
  writeLines("\nConverted ped file to binary ped.\n\n            ")


#Step 2
#Convert .fam file to .clust file
#Will add a third column - the clusters column - based on the first 3 characters of your Individual IDs
  famtoconvert <- utils::read.table(
    paste0(where.plink,paste0("BinaryPED",".fam")),
    quote = "",
    sep=" ",
    header=FALSE
    )

  #Extrac population names based on the _ seperation
  NameExtract <- substr(famtoconvert[,2],1,regexpr("_",famtoconvert[,2])-1)

  famtoconvert[,3] <- as.character(NameExtract)

  utils::write.table(x=famtoconvert[,1:3],file=paste0(where.plink,"ClusterFile.clust"),quote =FALSE,col.names = FALSE, row.names = FALSE)
  writeLines("\nCluster file created from .fam file.\n     ")
  writeLines("Creating frequency file using bed and cluster files.\n\n     ")
  remember.fam.plink<-paste0(where.plink,"BinaryPED.fam")
  bed.name=paste0(anal.name,".bed")
  remember.bed.plink<-paste0(where.plink,bed.name)


###STEP 3###
#####Now make the frequency cluster file in Plink
  if (Sys.info()["sysname"] != "Windows") {
    execute.PLINK2 <- paste0(callplink, "; ", "plink --noweb --bfile BinaryPED --freq --within ClusterFile.clust --out TreemixInput")
    system(execute.PLINK2)
  }

if (Sys.info()["sysname"] == "Windows") {
  execute.PLINK2<-paste0(callplink, " && ", "plink --noweb --bfile BinaryPED --freq --within ClusterFile.clust --out TreemixInput")
  shell(execute.PLINK2)
}

#Now gzip the output from Plink, then run this gzipped file through the Python script that comes with Treemix.
R.utils::gzip(filename=paste0(where.plink,"TreemixInput.frq.strat"))

#Clean files and return to path.

#Keep the intermediate conversion files.
    if(keep_inter){
      file.copy(from=paste0(where.pgdspider,"PGDtest.map"),to = paste0(path,"PGDtest.map"))
      file.copy(from=paste0(where.pgdspider,"PGDtest.ped"),to = paste0(path,"PGDtest.ped"))
      file.copy(from=paste0(where.plink,"ClusterFile.clust"),to=paste0(path,"ClusterFile.clust"))
      file.rename(from = paste0(path,"PGDtest.map"),to=paste0(path,"treemix_map.map"))
      file.rename(from = paste0(path,"PGDtest.ped"),to=paste0(path,"treemix_ped.ped"))
    }else{
      file.remove(paste0(where.pgdspider,"PGDtest.map"))
      file.remove(paste0(where.pgdspider,"PGDtest.ped"))
      file.remove(paste0(where.plink,"ClusterFile.clust"))
    }

  writeLines("\nCopying gzipped input file to path and removing unnecessary files\n")
    file.copy(from = paste0(where.plink,"TreemixInput.frq.strat.gz"), to = paste0(path,"TreemixInput.frq.strat.gz"))
    file.remove(paste0(where.pgdspider,"GP_PED.spid"))
    file.remove(paste0(where.pgdspider,"spider.conf.xml"))
    file.remove(paste0(where.pgdspider,"PGDSpider-cli.log"))
    file.remove(paste0(where.pgdspider,"GPD_for_PED_to_BED.txt"))
    file.remove(paste0(where.plink,"PGDtest.ped"))
    file.remove(paste0(where.plink,"PGDtest.map"))
    file.remove(paste0(where.plink,"BinaryPED.nosex"))
    file.remove(paste0(where.plink,"BinaryPED.log"))
    file.remove(paste0(where.plink,"BinaryPED.bed"))
    file.remove(paste0(where.plink,"BinaryPED.bim"))
    file.remove(paste0(where.plink,"BinaryPED.fam"))
    file.remove(paste0(where.plink,"TreemixInput.nosex"))
    file.remove(paste0(where.plink,"TreemixInput.log"))
    file.remove(paste0(where.plink,"TreemixInput.frq.strat.gz"))
    #should have .frq.gz as the extension after this
    writeLines("\nRun your new gzipped file through the Python script that comes with Treemix now. Then you are ready to put it into Treemix!")
return()
}

