# Get top loci
#' @title Creates a panel of the top n unlinked loci, and exports the list of loci
#' @description Extract the genotypes of individuals at the top n (by Fst) unlinked loci. The default threshold of r2>0.2 is employed for assigning 'linked' loci (default for plink).
#' @param genepop A file path to the genepop format file you wish to create your panel from
#' @param where.plink A file path to the PLINK installation folder.
#' @param where.pgdspider A file path to the PGDspider installation folder.
#' @param r2.threshold The minimum r^2 threshold to consider a pair of loci to be in LD
#' @param fst.threshold The minimum FST threshold required to retain a locus
#' @param ld.window Number of adjacent SNPs to compare each SNP against for LD - default is NULL, which translates to a window size of 99999, which essentially asks to compare each SNP against all others
#' @param ldpop A string which populations (default: "All") will be used to calculate linkage disequilibrium. Names must match names returned by genepop_detective().
#' @param allocate.pgd.ram An integer value in GB to specify the maximum amount of RAM to allocate to PGDspider. The default is 1 GB, which should be sufficient for most analyses.
#' @param return.workspace (default: TRUE) Logical query to return the output to the workspace
#' @param save.output Logical query (default: FALSE) to save the output to the same location as the file being analyzed. Each of the outputs of the function will be saved as a separate file with the file name of the orginal data appended with "_Linkages", "Loci_FST", and "Unlinked_Loci_FST" for the pairwise linked loci along with their r^2, all loci with their global Fst, and only the top unlinked loci with their Fst respectively.
#' @rdname genepop_toploci
#' @export
#' @importFrom stringr str_split str_detect
#' @importFrom plyr rbind.fill
#' @importFrom utils write.table
#' @importFrom utils read.table



genepop_toploci <- function(genepop, where.plink, where.pgdspider, r2.threshold = 0.2, fst.threshold = 0.05,  ld.window = NULL, ldpop = "All", allocate.pgd.ram = 1, return.workspace = TRUE, save.output = FALSE){
  stop("This function is depcrecated because the diveRsity package is no longer maintained")

  writeLines("Note that this function works with Plink 1.9 for now. If you have Plink 2.0, the function may not perform as intended.")

  path.start <- getwd()  ### where to write the files created by genepopedit to

  ## check for retrun options
  if(!save.output & !return.workspace){
    stop("No output option specified. At least one output parameter (save.output or return.workspace) must be TRUE.")
  }


  ## find the populations in the file
  pops.exist <- genepop_detective(genepop) ## see what populations are in the file

  if(allocate.pgd.ram%%1 != 0){
    stop("Please specify an integer GB value to allocate to PGDspider.")
  }

  #Set up ram allocation. Note that only 1024 gb of ram will work with Windows
  allocate.pgd.ram <- allocate.pgd.ram*1024

  if(Sys.info()["sysname"] == "Windows" & allocate.pgd.ram>1024){
    allocate.pgd.ram=1024
    writeLines("Note that currently PGDspider can only utilize ~1 GB of RAM on windows based operating systems. Periodically check back to https://github.com/rystanley/genepopedit for any updates to this limitation.
               ")}

  #Variable checks
  if(r2.threshold < 0 | r2.threshold > 1){
    stop("r^2 threshold must be a value between 0 and 1")
  }

  if(r2.threshold<0.2){
    writeLines("Linkage detection threshold is low (<0.2). Linkage will be classified at a higher frequency than default PLINK selection parameters.
               ")
  }

  if(is.null(ld.window)){
    ld.window = 999999 ### sets the LD window to essentially check every SNP pairwise
  }

  if(ld.window < 0){
    stop("LD window must be positive")
  }

  if(length(which(ldpop %in% c("All",pops.exist)))==0){
    stop(paste0("Parameter 'ldpop' must be a string of population names in the dataset, or ", "'All'. ", "Function stopped."),call. = FALSE)
  }

  if(fst.threshold < 0 | fst.threshold > 1){
    stop("FST threshold must be a value between 0 and 1")
  }



  # modify the path to play nice with spaces

  path.start.PGD <- gsub(x = path.start, pattern = " ", replacement = "\\")
  where.pgdspider.PGD <- gsub(x = where.pgdspider, pattern = " ", replacement = "\\ ", fixed = TRUE)

  #Calculate FST
  writeLines("Calculating locus-specific Fst\n\n")

  # FST.dat <- diveRsity::diffCalc(infile = genepop,outfile = NULL,fst = T,bs_locus = F)
  FST.dat <- FST.dat$std_stats[1:length(FST.dat$std_stats$loci)-1,]
  FSTs <- FST.dat$Fst

  FST.df <- data.frame(genepopedit::genepop_detective(genepop,variable="Loci"), FSTs)
  names(FST.df)[1] <- "loci"
  FST.df <- FST.df[base::order(FST.df$FSTs, decreasing = TRUE),
                   ]
  FST.Filter.Vec <- as.character(FST.df[which(FST.df$FSTs >=
                                                fst.threshold), 1])


  #Console message
  writeLines("Converting GENEPOP to MAP-PED format.")
  writeLines("

             ")
  writeLines("Warning messages are expected as part of conversion process using PGDspider.

             ")


  ### Will LD be calculated for All or specified populations.
  if(ldpop != "All"){

    subPOP <- as.character(ldpop)

    ## subset out the population in which LD is to be calculated - this will make a file, which will be deleted after
    subset_genepop(genepop = genepop, spop = subPOP, subs = FST.Filter.Vec, keep = TRUE, path = paste0(path.start, "/", "subset_for_LD.txt"))
    ## remember path to the file created by subset_genepop
    sub_data_path <- paste0(path.start, "/", "subset_for_LD.txt")
  }

  if(ldpop == "All"){

    popLDsubsetDF <- data.frame(op=pops.exist, paste0("PopA", rep(1:length(pops.exist)))) ## make a both the same

    subset_genepop_aggregate(genepop = genepop, agpopframe = popLDsubsetDF, path = paste0(path.start, "/", "subset_for_LD.txt"))
    sub_data_path <- paste0(path.start, "/", "subset_for_LD.txt")
    subset_genepop(genepop = sub_data_path, subs = FST.Filter.Vec, keep = TRUE, path = paste0(path.start, "/", "subset_for_LD.txt"))
    ## now rename

    subset_genepop_rename(genepop = sub_data_path, path = sub_data_path, nameframe = popLDsubsetDF,renumber=TRUE)

  }




  ### convert file to .ped and .map using PGD spider


  ## create spid file

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

  ## write spid file
  write(x = spid.file, file = paste0(path.start, "/", "hyb.spid"))

  ### move spid file to the PGDspider folder
  file.copy(from = paste0(path.start, "/hyb.spid"), to = where.pgdspider, overwrite = TRUE)
  remember.spidpath <- paste0(path.start, "/", "hyb.spid")
  ## move the input file as well to the same location as PGDspider - this makes this step so much easier
  file.copy(from <- sub_data_path, to = where.pgdspider, overwrite = TRUE)

  #### OS X LINUX call
  if(Sys.info()["sysname"] != "Windows"){
    ### create a string to call PGDspider
    input.file.call <- paste0("-inputfile subset_for_LD.txt")
    execute.SPIDER <- paste0("java -Xmx", allocate.pgd.ram, "m -Xms512m -jar PGDSpider2-cli.jar")
    spid.call <- "-spid hyb.spid"
    input.format <- "-inputformat GENEPOP"
    output.format <- "-outputformat PED"
    goto.spider <- paste0("cd ", where.pgdspider.PGD, "; ", execute.SPIDER)
    output.file.path <- paste0("-outputfile PGDtest.ped")
    ## string to run
    run.PGDspider <- paste0(goto.spider, " ", input.file.call, " ", input.format, " ", output.file.path, " ", output.format, " ", spid.call)

    ### run PGDspider through system
    system(run.PGDspider)
  } # END OSX LINUX if

  #### Windows call
  if(Sys.info()["sysname"] == "Windows"){
    ### create a string to call PGDspider
    input.file.call <- paste0("-inputfile subset_for_LD.txt")
    execute.SPIDER <- paste0("java -Xmx", allocate.pgd.ram, "m -Xms512m -jar PGDSpider2-cli.jar")
    spid.call <- "-spid hyb.spid"
    input.format <- "-inputformat GENEPOP"
    output.format <- "-outputformat PED"
    goto.spider <- paste0("cd ", where.pgdspider.PGD, " && ", execute.SPIDER)
    output.file.path <- paste0("-outputfile PGDtest.ped")
    ## string to run
    run.PGDspider <- paste0(goto.spider, " ", input.file.call, " ", input.format, " ", output.file.path, " ", output.format, " ", spid.call)

    ### run PGDspider through system
    shell(run.PGDspider)
  }

  ## move the created ped and map files to the PLINK folder
  ped.path <- paste0(where.pgdspider, "/", "PGDtest.ped")
  map.path <- paste0(where.pgdspider, "/", "PGDtest.map")

  file.copy(from = ped.path, to = where.plink, overwrite = TRUE)
  file.copy(from = map.path, to = where.plink, overwrite = TRUE)

  plink_ped_path <- paste0(where.plink, "/", "PGDtest.ped")
  plink_map_path <- paste0(where.plink, "/", "PGDtest.map")

  #Console message
  writeLines("

             ")
  writeLines("Calculating Linkage.

             ")

  ### prepare a string to call PLINK
  ### modify PLINK path so it plays nice with system
  where.plink.go <- gsub(x = where.plink, pattern = " ", replacement = "\\ ", fixed = TRUE)
  go.to.PLINK <- paste0("cd ", where.plink.go)

  ### OSX LINUX PLINK call
  if(Sys.info()["sysname"] != "Windows"){
    execute.PLINK <- paste0(go.to.PLINK, "; ", "./plink --file PGDtest --r2 --ld-window-r2 ", r2.threshold, " --ld-window ", ld.window, " --noweb")
    ### run PLINK through system
    system(execute.PLINK)
  }

  ### Windows PLINK CALL
  if(Sys.info()["sysname"] == "Windows"){
    execute.PLINK <- paste0(go.to.PLINK, " && ", "plink --file PGDtest --r2 --ld-window-r2 ", r2.threshold, " --ld-window ", ld.window, " --noweb")
    ### run PLINK through system
    shell(execute.PLINK)
  }

  #Console message
  writeLines("

             ")
  writeLines("Creating Fst optimized unlinked panel.

             ")

  ## copy the LD file created by PLINK to the working directory
  file.copy(from = paste0(where.plink, "plink.ld"), to = path.start)

  ## the format of the LD file is a bit messed up - it is not a regular matrix - have to modify the file a bit to get it to read in properly

  ## rename the file and change it to a txt
  file.rename(from = paste0(path.start, "/", "plink.ld"), to = paste0(path.start, "/", "plink.txt"))

  ## read the file in as characters - there are an uneven number of spaces between data, so have to remove those before can read as a table
  LDs <- readChar(con = "plink.txt", file.info("plink.txt")$size) ## read in as a long character string
  LDs_reform <- gsub(x = LDs, pattern = "\\ +", replacement = " ") ### remove more than 1 space
  ## save as a new .txt file, which is now formatted to be read into R easily
  write(x = LDs_reform, file = "LDsReform.txt")
  ## read back in
  Linked <- utils::read.table("LDsReform.txt", header = TRUE)
  ## keep only the needed columns
  Linked <- Linked[c("SNP_A", "SNP_B", "R2")]
  ### turn both columns into a single vector of unduplicated loci names
  ld.unique <- c(as.character(Linked$SNP_A), as.character(Linked$SNP_B))
  ld.unique <- as.character(ld.unique[which(duplicated(ld.unique)==FALSE)])

  if(length(ld.unique) > 1){

    FST.df2 <- FST.df[which(FST.df$loci %in% ld.unique),]
    FST.ld.ordered <- as.character(FST.df2[order(FST.df2$FSTs,decreasing=T),"loci"])

    holdlist <- list()
    for(i in FST.ld.ordered){
      hold <- c(i,as.character(Linked[which(Linked[,1]%in%i),2]),
                as.character(Linked[which(Linked[,2]%in%i),1]))

      hold2 <- NULL
      for(k in 1:length(hold)){hold2 <- c(hold2,which(FST.ld.ordered == hold[k]))}
      holdlist[[i]] <- hold2
    }

    linked.ranks.df2 <- plyr::rbind.fill(lapply(holdlist,function(y){as.data.frame(t(y),stringsAsFactors=FALSE)}))

    LRD2 <- as.matrix(linked.ranks.df2)

    highest1 <-  LRD2[1,which.min(LRD2[1,])]
    otherdat1 <- LRD2[1,-which.min(LRD2[1,])]
    otherdat1 <- otherdat1[!is.na(otherdat1)]

    to.keep <- FST.ld.ordered[Optimfunc(linked.ranks.df2)]
    to.drop <- setdiff(as.character(FST.df2$loci),to.keep)

    FST.df2$loci <- as.character(FST.df2$loci)
    FST.df$loci <- as.character(FST.df$loci)

    to.keep <- FST.df[-which(FST.df$loci %in% to.drop),]

    #Console message
    writeLines("

               ")
    writeLines("Writing output.

               ")

    Unlinked.panel <- FST.df[which(FST.df$loci %in% to.keep$loci),]


    your.panel <- FST.df
    your.panel_un <- Unlinked.panel[which(Unlinked.panel$FSTs>=fst.threshold),]
    your.panel <- your.panel[order(your.panel$FSTs,decreasing=TRUE),]
    your.panel_unlinked <- your.panel_un[order(your.panel_un$FSTs,decreasing = TRUE),]
    Linked.df <- Linked
  } ### END IF There are loci in LD

  if(length(ld.unique) < 1){ # if not linked loci
    writeLines("Writing output")
    writeLines("No linked loci detected, all loci returned.")
    your.panel <- FST.df
    your.panel <- your.panel[order(your.panel$FSTs,decreasing=TRUE),]
    your.panel_unlinked <- your.panel
    Linked.df <- NA
  }

  #clean up temporary files used in the analysis.
  file.remove(remember.spidpath)
  file.remove(sub_data_path)
  file.remove(ped.path)
  file.remove(map.path)
  file.remove(paste0(where.pgdspider, "/hyb.spid"))
  file.remove(plink_map_path)
  file.remove(plink_ped_path)
  file.remove(paste0(path.start, "/plink.txt"))
  file.remove(paste0(path.start, "/LDsReform.txt"))
  file.remove(paste0(where.pgdspider,"/subset_for_LD.txt"))

  #Console message
  writeLines("

             ")
  writeLines("Process Completed.")

  ## return output
  Output <- list(Linkages=Linked.df,
                 Fst=your.panel,
                 Fst_Unlinked=your.panel_unlinked
  )



  if(save.output == TRUE){

    write_path <- gsub(x = genepop, pattern = ".txt", replacement = "")

    utils::write.table(x = Output$Linkages, file = paste0(write_path, "_linkages.txt"),row.names = FALSE,quote=FALSE)
    utils::write.table(x = Output$Fst, file = paste0(write_path, "_Loci_FST.txt"),row.names=FALSE,quote = FALSE)
    utils::write.table(x = Output$Fst_Unlinked, file = paste0(write_path, "_Unlinked_Loci_FST.txt"),row.names=FALSE,quote = FALSE)

  }

   if(return.workspace == TRUE){

    return(Output)
  }

  }
