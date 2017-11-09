# Filter genepop file for MAF
#' @title Filters genepop file for MAF
#' @description Converts genepop file to PED and filters for MAF using plink and returns genepop file with only SNPs that pass filter
#' @param genepop A file path to the genepop format file you wish to filter.
#' @param where.plink A file path to the PLINK installation folder.
#' @param where.pgdspider A file path to the PGDspider installation folder.
#' @param maf Minor allele frequency cutoff (default = 0.05)
#' @rdname genepop_filter_maf
#' @importFrom data.table fread
#' @export


genepop_filter_maf <- function(genepop, where.plink, where.PGDspider, maf=0.05, path) {


  #Using the new genepop files generated - run through plink
  #Create spid file

  #Using the new genepop files generated - run through plink
  #Create spid file
  spidTop <- "# spid-file generated: Thu Mar 10 09:40:22 AST 2016
  # GENEPOP Parser questions
  PARSER_FORMAT=GENEPOP
  # Enter the size of the repeated motif (same for all loci: one number; different: comma separated list (e.g.: 2,2,3,2):
  GENEPOP_PARSER_REPEAT_SIZE_QUESTION=
  # Select the type of the data:
  GENEPOP_PARSER_DATA_TYPE_QUESTION=SNP
  # How are Microsat alleles coded?
  GENEPOP_PARSER_MICROSAT_CODING_QUESTION=REPEATS
  # PED Writer questions
  WRITER_FORMAT=PED
  # Save MAP file"
  map.loc <- paste0("PED_WRITER_MAP_FILE_QUESTION= ", "PGDtest")
  spidBottom <- "# Replacement character for allele encoded as 0 (0 encodes for missing data in PED):
  PED_WRITER_ZERO_CHAR_QUESTION=
  # Specify the locus/locus combination you want to write to the PED file:
  PED_WRITER_LOCUS_COMBINATION_QUESTION=
  # Do you want to save an additional MAP file with loci information?
  PED_WRITER_MAP_QUESTION=true"
  spid.file <- c(spidTop, map.loc, spidBottom)

  #Console message
  writeLines("Converting GENEPOP to PED format.")
  writeLines("

             ")
  writeLines("Warning messages are expected as part of conversion process using PGDspider.

             ")



  ###write spid file genepop conversion
  write(x = spid.file, file = paste0(where.PGDspider, "/", "hyb.spid"))

  file.copy(from = genepop, to = paste0(where.PGDspider, "/", "genepop_test.txt"), overwrite = T)

  ###convert Genepop to PED using PGD spider
  where.PGDspider.PGD <- gsub(x = where.PGDspider, pattern = " ", replacement = "\\ ", fixed = TRUE)


  #### OS X and LINUX CALL

  if(Sys.info()["sysname"] != "Windows"){

    ### create a string to call PGDspider
    input.file.call <- "-inputfile genepop_test.txt"
    execute.SPIDER <- paste0("java -Xmx1024", "m -Xms512m -jar PGDSpider2-cli.jar")
    spid.call <- "-spid hyb.spid"
    input.format <- "-inputformat GENEPOP"
    output.format <- "-outputformat PED"

    goto.spider <- paste0("cd ", where.PGDspider.PGD, "; ", execute.SPIDER)
    output.file.path <- "-outputfile PGDtest.ped"
    ## string to run
    run.PGDspider <- paste0(goto.spider, " ", input.file.call, " ", input.format, " ", output.file.path, " ", output.format, " ", spid.call)
    ### run PGDspider through system
    system(run.PGDspider)

  } # End MAC LINUX IF

  #### Windows call

  if(Sys.info()["sysname"] == "Windows"){

    ### create a string to call PGDspider
    input.file.call <- "-inputfile genepop_test.txt"
    execute.SPIDER <- paste0("java -Xmx1024", "m -Xms512m -jar PGDSpider2-cli.jar")
    spid.call <- "-spid hyb.spid"
    input.format <- "-inputformat GENEPOP"
    output.format <- "-outputformat PED"
    goto.spider <- paste0("cd ", where.PGDspider.PGD, " && ", execute.SPIDER)
    output.file.path <- "-outputfile PGDtest.ped"
    ## string to run
    run.PGDspider <- paste0(goto.spider, " ", input.file.call, " ", input.format, " ", output.file.path, " ", output.format, " ", spid.call)
    ### run PGDspider through system
    shell(run.PGDspider)

  } # End WINDOWS IF


  remember.PEDpath.PGD <- paste0(where.PGDspider, "PGDtest.ped")
  remember.MAPpath.PGD <- paste0(where.PGDspider, "PGDtest.map")


  ###move to plink folder
  ped.path <- paste0(where.PGDspider, "/", "PGDtest.ped")
  map.path <- paste0(where.PGDspider, "/", "PGDtest.map")
  file.copy(from = ped.path, to = where.plink, overwrite = TRUE)

  file.copy(from = map.path, to = where.plink, overwrite = TRUE)


  plink_ped_path <- paste0(where.plink, "/", "PGDtest.ped")
  plink_map_path <- paste0(where.plink, "/", "PGDtest.map")

  ##Run plink for LD
  ####Plink creates file in plink.ld with LD values. Rename this file and move it to folder to have results later (will continue to over write plink.ld for each time plink is run)

  writeLines("Running plink for maf ")
  writeLines("
             ")

  where.plink.go <- gsub(x = where.plink, pattern = " ", replacement = "\\ ", fixed = TRUE)
  go.to.PLINK <- paste0("cd ", where.plink.go)



  ### OSX LINUX PLINK call
  if(Sys.info()["sysname"] != "Windows"){
    execute.PLINK <- paste0(go.to.PLINK, "; ", "./plink --file PGDtest --maf ", maf, " --make-bed")


    ### run PLINK through system
    system(execute.PLINK)
  }

  ### Windows PLINK CALL
  if(Sys.info()["sysname"] == "Windows"){
    execute.PLINK <- paste0(go.to.PLINK, " && ", "plink --file PGDtest --maf ", maf, " --make-bed")
    ### run PLINK through system
    shell(execute.PLINK)
  }

  #Make list of loci that passed MAF filter
  Loci_list_maf=data.table::fread(paste(where.plink, "plink.bim", sep="/"), header=F)
  Loci_keep=as.character(Loci_list_maf$V2)

  #Write as genepop file with only snps that passed filter
  genepopedit::subset_genepop(genepop, subs = Loci_keep, keep = T, path = paste0(path))


}
