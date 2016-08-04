# Genepop -> Colony
#' @title Convert Genepop to Colony format.
#' @description Function to convert Genepop to Colony format.
#' @param genepop the genepop data to be manipulated. This is read in as a complete file path.
#' This will the standard genepop format with a the first n+1 rows corresponding the the n loci names,
#' or a single comma deliminated row of loci names followed by the locus data. Populations are
#' separated by "Pop". Each individual ID is linked to the locus data by " ,  " (space,space space) and is read in as
#' as a single row (character).
#' @param where.plink A file path to the PLINK installation folder.
#' @param where.pgdspider A file path to the PGDspider installation folder.
#' @param denote.missing The value that denotes missing data in your input file
#' @param allocate.pgd.ram An integer value in GB to specify the maximum amount of RAM to allocate to PGDspider. The default is 1 GB, which should be sufficient for most analyses.
#' @param path file path to directory where the Colony files (4) will be saved.
#' @rdname genepop_colony
#' @importFrom stringr str_split
#' @importFrom stringi stri_sub
#' @importFrom utils write.table
#' @importFrom utils read.table
#' @export


##
genepop_colony <- function(genepop, where.plink, where.pgdspider, denote.missing = "000", allocate.pgd.ram = 1,path){

   path.start <- getwd()  ### where to write the files created by genepopedit to

   ## Warning messages
       if(allocate.pgd.ram%%1 != 0){
        stop("Please specify an integer GB value to allocate to PGDspider.")
      }

    #Set up ram allocation. Note that only 1024 gb of ram will work with Windows
      allocate.pgd.ram <- allocate.pgd.ram*1024

      if(Sys.info()["sysname"] == "Windows" & allocate.pgd.ram>1024){
        allocate.pgd.ram=1024
        writeLines("Note that currently PGDspider can only utilize ~1 GB of ram on windows based operating systems. Periodically check back to https://github.com/rystanley/genepopedit for any updates to this limitation.
                   ")}

  ## get the name of the file specified by genepop
      GeneNAME <-  filename_check(X = genepop)

  ## move the file specified by genepop to the working directory
      file.copy(from = genepop, to = path.start)
      remember.GenePop <- paste0(path.start, "/", GeneNAME)

  ## hold version of the genepop in a temporary folder ## needed when getwd() is the same as the locaiton of genepop
      dir.create("temporaryCOLONYcontainer")
      file.copy(from = genepop, to = "temporaryCOLONYcontainer/")

  ## rename the file in the working directory
      file.rename(from = remember.GenePop, paste0(path.start, "/genepop_colony.txt"))
      remember.GenePopColony <- paste0(path.start, "/genepop_colony.txt")

  ## add the original back to the directory
      file.copy(paste("temporaryCOLONYcontainer",GeneNAME,sep="/"),paste(path.start,GeneNAME,sep="/"))

  ## delete temporary container
      unlink("temporaryCOLONYcontainer", recursive = TRUE, force = FALSE)

  ### convert file to .ped and .map using PGD spider

    #Console message
      writeLines("Converting GENEPOP to MAP-PED format.")
      writeLines("

                 ")
      writeLines("Warning messages are expected as part of conversion process using PGDspider.

                 ")

  ## check that the path to PGDspider has been specified correctly
    where.pgdspider <- path_ending(path = where.pgdspider)

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

      map.loc <- paste0("PED_WRITER_MAP_FILE_QUESTION= ", "Colony_Out")

      spidBottom <- "# Replacement character for allele encoded as 0 (0 encodes for missing data in PED):
      PED_WRITER_ZERO_CHAR_QUESTION=
      # Specify the locus/locus combination you want to write to the PED file:
      PED_WRITER_LOCUS_COMBINATION_QUESTION=
      # Do you want to save an additional MAP file with loci information?
      PED_WRITER_MAP_QUESTION=true"

      spid.file <- c(spidTop, map.loc, spidBottom)

  ## write spid file to working directory
      write(x = spid.file, file = paste0(path.start, "/", "GenePopColony.spid"))
      remember.spidpath.WD <- paste0(path.start, "/", "GenePopColony.spid")

  ### move spid file to the PGDspider folder
      file.copy(from = paste0(path.start, "/GenePopColony.spid"), to = where.pgdspider, overwrite = TRUE)
      remember.spidpath.PGD <- paste0(where.pgdspider, "GenePopColony.spid")

  ## Convert Genepop to MAP/PED
      ## move the input file as well to the same location as PGDspider - this makes this step so much easier
      file.copy(from <- remember.GenePopColony, to = where.pgdspider, overwrite = TRUE)
      remember.GenePopColony.PGD <- paste0(where.pgdspider, "/genepop_colony.txt")

      where.pgdspider.PGD <- gsub(x = where.pgdspider, pattern = " ", replacement = "\\ ", fixed = TRUE)

      #### OS X LINUX call
      if(Sys.info()["sysname"] != "Windows"){
        ### create a string to call PGDspider
        input.file.call <- paste0("-inputfile genepop_colony.txt")
        execute.SPIDER <- paste0("java -Xmx", allocate.pgd.ram, "m -Xms512m -jar PGDSpider2-cli.jar")
        spid.call <- "-spid GenePopColony.spid"
        input.format <- "-inputformat GENEPOP"
        output.format <- "-outputformat PED"
        goto.spider <- paste0("cd ", where.pgdspider.PGD, "; ", execute.SPIDER)
        output.file.path <- paste0("-outputfile Colony_Out.ped")
        ## string to run
        run.PGDspider <- paste0(goto.spider, " ", input.file.call, " ", input.format, " ", output.file.path, " ", output.format, " ", spid.call)

        ### run PGDspider through system
        system(run.PGDspider)
      } # END OSX LINUX if

      #### Windows call
      if(Sys.info()["sysname"] == "Windows"){
        ### create a string to call PGDspider
        input.file.call <- paste0("-inputfile genepop_colony.txt")
        execute.SPIDER <- paste0("java -Xmx", allocate.pgd.ram, "m -Xms512m -jar PGDSpider2-cli.jar")
        spid.call <- "-spid GenePopColony.spid"
        input.format <- "-inputformat GENEPOP"
        output.format <- "-outputformat PED"
        goto.spider <- paste0("cd ", where.pgdspider.PGD, " && ", execute.SPIDER)
        output.file.path <- paste0("-outputfile Colony_Out.ped")
        ## string to run
        run.PGDspider <- paste0(goto.spider, " ", input.file.call, " ", input.format, " ", output.file.path, " ", output.format, " ", spid.call)

        ### run PGDspider through system
        shell(run.PGDspider)
      }

  ## move the created ped and map files to the PLINK folder
      ped.path <- paste0(where.pgdspider, "/", "Colony_Out.ped")
      map.path <- paste0(where.pgdspider, "/", "Colony_Out.map")

      file.copy(from = ped.path, to = where.plink, overwrite = TRUE)
      file.copy(from = map.path, to = where.plink, overwrite = TRUE)

      plink_ped_path <- paste0(where.plink, "/", "Colony_Out.ped")
      plink_map_path <- paste0(where.plink, "/", "Colony_Out.map")

## Call plink.
  #Console message
    writeLines("

               ")
    writeLines("Estimating missing data using PLINK.

               ")

  ## check that PLINK folder has been specified correctly
    where.plink <- path_ending(path = where.plink)

  ### prepare a string to call PLINK
  ### modify PLINK path so it plays nice with system
    where.plink.go <- gsub(x = where.plink, pattern = " ", replacement = "\\ ", fixed = TRUE)
    go.to.PLINK <- paste0("cd ", where.plink.go)

  ### OSX LINUX PLINK call
    if(Sys.info()["sysname"] != "Windows"){
      execute.PLINK <- paste0(go.to.PLINK, "; ", "./plink --file Colony_Out --missing --noweb")
      ### run PLINK through system
      system(execute.PLINK)
    }

  ### Windows PLINK CALL
    if(Sys.info()["sysname"] == "Windows"){
      execute.PLINK <- paste0(go.to.PLINK, " && ", "plink --file Colony_Out --missing --noweb")
      ### run PLINK through system
      shell(execute.PLINK)
    }

  ### PLINK is going to make a few files, will want to remove them later just to be nice
    remember.nosex <- paste0(where.plink, "plink.nosex")
    remember.log <- paste0(where.plink, "plink.log")
    remember.imiss <- paste0(where.plink, "plink.imiss")
    remember.lmiss.PLINK <- paste0(where.plink, "plink.lmiss")

  ## copy the LOCI MISSING file created by PLINK to the working directory
    file.copy(from = paste0(where.plink, "plink.lmiss"), to = path.start)
    remeber.lmiss.WD <- paste0(path.start, "/plink.limiss")

  ## rename the file and change it to a txt
    file.rename(from = paste0(path.start, "/", "plink.lmiss"), to = paste0(path.start, "/", "plink.txt"))
    remember.lmiss.WD.renamed <- paste0(path.start, "/", "plink.txt")

  ## read the file in as characters - there are an uneven number of spaces between data, so have to remove those before can read as a table
    Missings <- readChar(con = "plink.txt", file.info("plink.txt")$size) ## read in as a long character string
    Missings_reform <- gsub(x = Missings, pattern = "\\ +", replacement = " ") ### remove more than 1 space
  ## save as a new .txt file, which is now formatted to be read into R easily
    write(x = Missings_reform, file = "MissingsReform.txt")
  ## read back in
    MissingDataReform <- utils::read.table("MissingsReform.txt", header = TRUE)
    remember.missingreform <- paste0(path.start, "/MissingsReform.txt")

  ### only want to keep the LOCI NAMES, and the PROPORTION MISSING, will duplicate loci names so can add a row of 0 to indicate the loci are codominant
    MissingDataReform <- MissingDataReform[, c("SNP", "SNP", "F_MISS", "F_MISS")] #### NOT SURE THIS IS THE BEST VALUE <- BUT, I have seen the error rate just duplicated
  ## transpose to make the correct format
    MissingDataReform <- t(MissingDataReform)
    yourloci <- MissingDataReform[1, ]
    lociconvert <- paste0("Loci_", 1:length(yourloci))

  ## change the names of the loci to Loci plus number because colony can't deal with >20 characters
  ## add a row of 0 to denote codominant
    MissingDataReform[1, ] <- lociconvert
    MissingDataReform[2, ] <- 0
    MissingDataReform <- data.frame(MissingDataReform)

  #Tabularize and format data
    flatdat <- genepop_flatten(genepop = genepop)
    indivs <- data.frame(Ind_Names = flatdat$SampleID, Number = 1:nrow(flatdat))
    colnames(indivs) <- c("Individual Names", "Corresponding Number")
    lociout <- data.frame(A= yourloci, B = lociconvert)
    colnames(lociout) <- c("Loci Names", "Corresponding Number")
    coldat <- flatdat[, -c(1,2,3)]
    blankmat <- matrix(nrow = nrow(coldat), ncol = (length(coldat)*2))
    locidim <- nchar(x = as.character(coldat[1,1]))/2

  #Define apply functions to split loci into alleles
    LHS <- function(x){stringi::stri_sub(as.character(x), 1, locidim)}
    RHS <- function(x){stringi::stri_sub(as.character(x), locidim+1, 2*locidim)}
    leftALLELE <- apply(X = coldat, MARGIN = 2, FUN = LHS)
    rightALLELE <- apply(X = coldat, MARGIN = 2, FUN = RHS)

  #create a blank matrix and populate alleles (Loci now defined in consecutive allele pairs (columns))
    RH_vec <- which(1:ncol(blankmat) %% 2 == 0)
    LH_vec <- which(!1:ncol(blankmat) %% 2 == 0)
    blankmat[,RH_vec] <- rightALLELE
    blankmat[,LH_vec] <- leftALLELE
    blankmat <- cbind(1:nrow(blankmat), blankmat)

  ## Missing data must be coded 0
    blankmat[which(blankmat == denote.missing)] = "0"

  ## if any individual has all missing data, it must be removed or else Colony won't run - find them and remove
    lastcol <- ncol(blankmat)
    colsumBLANKMAT <- function(k){sum(as.numeric(k))}

  ### Logical to test that there ARE individuals with all missing first
  if(length(which(apply(X = blankmat[, 2:lastcol], MARGIN = 1, FUN = colsumBLANKMAT) == 0)) >1){
    remove_Rows <- which(apply(X = blankmat[, 2:lastcol], MARGIN = 1, FUN = colsumBLANKMAT) == 0)
    blankmat <- blankmat[-remove_Rows, ]
  }

  ### Save files
    #Check file path
      path <- path_ending(path)
    ## Loci
      utils::write.table(x = lociout, file = paste0(path, unlist(stringr::str_split(string = GeneNAME, pattern = ".txt"))[1], "_Loci_Conversion.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
    ## Individuals
      utils::write.table(x = indivs, file = paste0(path, unlist(stringr::str_split(string = GeneNAME, pattern = ".txt"))[1], "_Individual_Conversion.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
    ## Marker Type & Error Rate
      utils::write.table(x = MissingDataReform, file = paste0(path, unlist(stringr::str_split(string = GeneNAME, pattern = ".txt"))[1], "_MarkerTypeErrorRT.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
    ## Genotypes
      utils::write.table(x = blankmat, file = paste0(path, unlist(stringr::str_split(string = GeneNAME, pattern = ".txt"))[1], "_GENOTYPES.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")

  ## Clean up intermediate files.
      file.remove(remember.GenePopColony)
      file.remove(remember.spidpath.PGD)
      file.remove(plink_ped_path)
      file.remove(plink_map_path)
      file.remove(ped.path)
      file.remove(map.path)
      file.remove(remember.nosex)
      file.remove(remember.log)
      file.remove(remember.imiss)
      file.remove(remember.lmiss.PLINK)
      file.remove(remember.lmiss.WD.renamed)
      file.remove(remember.spidpath.WD)
      file.remove(remember.GenePopColony.PGD)
      file.remove(remember.missingreform)
  }
