# Get top loci
#' @title Creates a panel of the top n unlinked loci, and exports the list of loci
#' @description Extract the genotypes of individuals at the top n (by Fst) unlinked loci.
#' @param GenePop A file path to the GENEPOP format file you wish to create your panel from
#' @param LDpop A string which denotes which of the two populations you wish to calculate linkage disequilibrium from. The options are "Pop1" or "Pop2", or "Both" (default) if the LD is to be calculated based on both populations. The order of population is based on the order by which they appear in the genepopfile. refer to genepop_detetive() to check which population is considered 'Pop1' and 'Pop2'.
#' @param panel.size An integer number of loci to include in the panel. If not specified all loci will be returned,
#' @param where.PLINK A file path to the PLINK installation folder.
#' @param where.PGDspider A file path to the PGDspider installation folder.
#' @rdname genepop_toploci
#' @export
#' @importFrom hierfstat read.fstat wc
#' @importFrom stringr str_split str_extract
#' @importFrom plyr rbind.fill


genepop_toploci <- function(GenePop, LDpop = "Both", panel.size=NULL, where.PLINK, where.PGDspider){


  #Parameter limits
      nLOCI <- length(genepop_detective(GenePop,"Loci"))
  #Variable checks
      if(length(which(LDpop %in% c("Pop1","Pop2","Both")))==0){
        stop("Parameter 'LDpop' must be a string of value Pop1, Pop2 or Both. Function stopped.",call. = FALSE)
      }

      if(is.null(panel.size)){
        writeLines("Panel size not specified, Fst calculations for all un-linked loci will be returned")
        panel.size <- nLOCI
        }

      if(panel.size>nLOCI){
        writeLines(paste0("Panel size selected is larger than available loci in file: ",GenePop,". panel.size set to the maximum number of loci: ",nLOCI))
        panel.size <- nLOCI
      }

      path.start <- getwd()  ### where to write the files created by genepopedit to

      pops.exist <- genepop_detective(GenePop) ## see what populations are in the file

  ## Function must have exactly two populations to work - fail if not
      if(length(pops.exist) != 2){
        stop("File must contain two populations. See subset_genepop() and-or subset_genepop_rename() for data manipulation options")
      }

  ### Will LD be calculated for both populations or a specific population.
      if(LDpop != "Both"){

        ## makes a dataframe to match whatever name is given to the populations to Pop1 and Pop2
        popLDsubsetDF <- data.frame(op=pops.exist, rename=c("Pop1", "Pop2")) ## make a DF to pick the right pop from the file

        subPOP <- popLDsubsetDF[which(popLDsubsetDF[,2] == LDpop),1] ### get the name of the pop to be subsetted

        subPOP <- as.character(subPOP)

        ## subset out the population in which LD is to be calculated - this will make a file, which will be deleted after
        subset_genepop(GenePop = GenePop, sPop = subPOP, keep = TRUE, path = paste0(path.start, "/", "subset_for_LD.txt"))
        ## remember path to the file created by subset_genepop
        sub_data_path <- paste0(path.start, "/", "subset_for_LD.txt")
      }

        if(LDpop == "Both"){

          popLDsubsetDF <- data.frame(op=pops.exist, rename=c("Pop1", "Pop1")) ## make a both the same

          subset_genepop_aggregate(GenePop = GenePop, agPopFrame = popLDsubsetDF, path = paste0(path.start, "/", "subset_for_LD.txt"))
          sub_data_path <- paste0(path.start, "/", "subset_for_LD.txt")
          ## now rename

          subset_genepop_rename(GenePop = sub_data_path, path = sub_data_path, nameframe = popLDsubsetDF,renumber=TRUE)

        }

      writeLines("Calculating Fst")

  ### convert to FST format using PGDspider

  # modify the path to play nice with spaces

      path.start.PGD <- gsub(x = path.start, pattern = " ", replacement = "\\")
      where.PGDspider.PGD <- gsub(x = where.PGDspider, pattern = " ", replacement = "\\ ", fixed = TRUE)

  GP_FSTAT_SPID <- "# spid-file generated: Fri Apr 08 10:53:23 ADT 2016

  # GENEPOP Parser questions
  PARSER_FORMAT=GENEPOP

  # Enter the size of the repeated motif (same for all loci: one number; different: comma separated list (e.g.: 2,2,3,2):
  GENEPOP_PARSER_REPEAT_SIZE_QUESTION=
  # Select the type of the data:
  GENEPOP_PARSER_DATA_TYPE_QUESTION=SNP
  # How are Microsat alleles coded?
  GENEPOP_PARSER_MICROSAT_CODING_QUESTION=REPEATS

  # FSTAT Writer questions
  WRITER_FORMAT=FSTAT

  # Specify which data type should be included in the FSTAT file  (FSTAT can only analyze one data type per file):
  FSTAT_WRITER_DATA_TYPE_QUESTION=SNP
  # Save label file
  FSTAT_WRITER_LABEL_FILE_QUESTION=
  # Do you want to save an additional file with labels (population names)?
  FSTAT_WRITER_INCLUDE_LABEL_QUESTION=false
  # Specify the locus/locus combination you want to write to the FSTAT file:
  FSTAT_WRITER_LOCUS_COMBINATION_QUESTION=
  "
    #write conversion parameter file (".spid")
        write(x = GP_FSTAT_SPID, file = paste0(path.start, "/", "GP_FSTAT.spid"))

    ### move spid file to the PGDspider folder
        file.copy(from = paste0(path.start, "/GP_FSTAT.spid"), to = where.PGDspider, overwrite = TRUE)
        remember.spidpath <- paste0(path.start, "/", "GP_FSTAT.spid")
    ## move the input file as well to the same location as PGDspider - this makes this step so much easier
        file.copy(from <- GenePop, to = where.PGDspider, overwrite = TRUE)
        GenePop.name <- stringr::str_split(string = GenePop, pattern = "/")
        GenePop.name <- unlist(GenePop.name)
        GenePop.name <- GenePop.name[grep(x = GenePop.name, pattern = ".txt")]
        file.rename(from = paste0(where.PGDspider, GenePop.name), to = paste0(where.PGDspider, "GPD_for_GET_TOP_LOC.txt"))
        remember.TOPLOC.path <- paste0(where.PGDspider, "GPD_for_GET_TOP_LOC.txt")


  #### OS X and LINUX CALL

      if(Sys.info()["sysname"] != "Windows"){

        ### create a string to call PGDspider
        input.file.call <- "-inputfile GPD_for_GET_TOP_LOC.txt"
        execute.SPIDER <- "java -Xmx1024m -Xms512m -jar PGDSpider2-cli.jar"
        spid.call <- "-spid GP_FSTAT.spid"
        input.format <- "-inputformat GENEPOP"
        output.format <- "-outputformat FSTAT"
        goto.spider <- paste0("cd ", where.PGDspider.PGD, "; ", execute.SPIDER)
        output.file.path <- "-outputfile for_FST.txt"
        ## string to run
        run.PGDspider <- paste0(goto.spider, " ", input.file.call, " ", input.format, " ", output.file.path, " ", output.format, " ", spid.call)


        ### run PGDspider through system
        system(run.PGDspider)

      } # End MAC LINUX IF


  #### Windows call

      if(Sys.info()["sysname"] == "Windows"){

        ### create a string to call PGDspider
        input.file.call <- "-inputfile GPD_for_GET_TOP_LOC.txt"
        execute.SPIDER <- "java -Xmx1024m -Xms512m -jar PGDSpider2-cli.jar"
        spid.call <- "-spid GP_FSTAT.spid"
        input.format <- "-inputformat GENEPOP"
        output.format <- "-outputformat FSTAT"
        goto.spider <- paste0("cd ", where.PGDspider.PGD, " && ", execute.SPIDER)
        output.file.path <- "-outputfile for_FST.txt"
        ## string to run
        run.PGDspider <- paste0(goto.spider, " ", input.file.call, " ", input.format, " ", output.file.path, " ", output.format, " ", spid.call)


        ### run PGDspider through system
        shell(run.PGDspider)

      } # End WINDOWS IF



  ### move the FSTAT format file back to the working directory
      file.copy(from = paste0(where.PGDspider, "/for_FST.txt"), to = path.start, overwrite = TRUE)
  ## remember the path of the file created by genepop_fstat
      fst_data_path <- paste0(path.start, "/", "for_FST.txt")


  ### read in the FSTAT formatted file
     for.fst <- hierfstat::read.fstat("for_FST.txt")
  ## calculate Fst
     FST.dat <- hierfstat::wc(for.fst)
  ### get the Fst values
      FSTs <- FST.dat$per.loc$FST
  ### create a dataframe that is the names of the Loci, and their corresponding Fst
      FST.df <- data.frame(colnames(for.fst)[-1], FSTs)
      names(FST.df)[1] <- "loci"
  ## reorder the dataframe from highest to lowest Fst
      FST.df <- FST.df[base::order(FST.df$FSTs, decreasing = TRUE),]

  writeLines("Calculating Linkage")
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
      file.copy(from = paste0(path.start, "/hyb.spid"), to = where.PGDspider, overwrite = TRUE)
      remember.spidpath <- paste0(path.start, "/", "hyb.spid")
  ## move the input file as well to the same location as PGDspider - this makes this step so much easier
      file.copy(from <- sub_data_path, to = where.PGDspider, overwrite = TRUE)


  #### OS X LINUX call
      if(Sys.info()["sysname"] != "Windows"){
        ### create a string to call PGDspider
        input.file.call <- paste0("-inputfile subset_for_LD.txt")
        execute.SPIDER <- "java -Xmx1024m -Xms512m -jar PGDSpider2-cli.jar"
        spid.call <- "-spid hyb.spid"
        input.format <- "-inputformat GENEPOP"
        output.format <- "-outputformat PED"
        goto.spider <- paste0("cd ", where.PGDspider.PGD, "; ", execute.SPIDER)
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
      execute.SPIDER <- "java -Xmx1024m -Xms512m -jar PGDSpider2-cli.jar"
      spid.call <- "-spid hyb.spid"
      input.format <- "-inputformat GENEPOP"
      output.format <- "-outputformat PED"
      goto.spider <- paste0("cd ", where.PGDspider.PGD, " && ", execute.SPIDER)
      output.file.path <- paste0("-outputfile PGDtest.ped")
      ## string to run
      run.PGDspider <- paste0(goto.spider, " ", input.file.call, " ", input.format, " ", output.file.path, " ", output.format, " ", spid.call)

      ### run PGDspider through system
      shell(run.PGDspider)
    }

  ## move the created ped and map files to the PLINK folder
      ped.path <- paste0(where.PGDspider, "/", "PGDtest.ped")
      map.path <- paste0(where.PGDspider, "/", "PGDtest.map")

      file.copy(from = ped.path, to = where.PLINK, overwrite = TRUE)
      file.copy(from = map.path, to = where.PLINK, overwrite = TRUE)

      plink_ped_path <- paste0(where.PLINK, "/", "PGDtest.ped")
      plink_map_path <- paste0(where.PLINK, "/", "PGDtest.map")


  ### prepare a string to call PLINK
  ### modify PLINK path so it plays nice with system
      where.PLINK.go <- gsub(x = where.PLINK, pattern = " ", replacement = "\\ ", fixed = TRUE)
      go.to.PLINK <- paste0("cd ", where.PLINK.go)

  ### OSX LINUX PLINK call
      if(Sys.info()["sysname"] != "Windows"){
        execute.PLINK <- paste0(go.to.PLINK, "; ", "./plink --file PGDtest --r2 --noweb")
        ### run PLINK through system
        system(execute.PLINK)
      }

  ### Windows PLINK CALL
      if(Sys.info()["sysname"] == "Windows"){
        execute.PLINK <- paste0(go.to.PLINK, " && ", "plink --file PGDtest --r2 --noweb")
        ### run PLINK through system
        shell(execute.PLINK)
      }

  ## copy the LD file created by PLINK to the working directory
      file.copy(from = paste0(where.PLINK, "plink.ld"), to = path.start)

  ## the format of the LD file is a bit messed up - it is not a regular matrix - have to modify the file a bit to get it to read in properly

  ## rename the file and change it to a txt
      file.rename(from = paste0(path.start, "/", "plink.ld"), to = paste0(path.start, "/", "plink.txt"))

  ## read the file in as characters - there are an uneven number of spaces between data, so have to remove those before can read as a table
      LDs <- readChar(con = "plink.txt", file.info("plink.txt")$size) ## read in as a long character string
      LDs_reform <- gsub(x = LDs, pattern = "\\ +", replacement = " ") ### remove more than 1 space
  ## save as a new .txt file, which is now formatted to be read into R easily
      write(x = LDs_reform, file = "LDsReform.txt")
  ## read back in
      Linked <- read.table("LDsReform.txt", header = TRUE)
  ## keep only the needed columns
      Linked <- Linked[c("SNP_A", "SNP_B", "R2")]
  ### turn both columns into a single vector of unduplicated loci names
      ld.unique <- c(as.character(Linked$SNP_A), as.character(Linked$SNP_B))
      ld.unique <- as.character(ld.unique[which(duplicated(ld.unique)==FALSE)])
  # head(FST.df)

  ## get a vector which is the names of the loci in order of highest to lowest Fst
      FST.order.vec <- c(as.character(FST.df$loci))
  ## get those loci that appear in both lists
  # which(FST.order.vec %in% ld.unique)

      loci.in.LD <- which(FST.order.vec %in% ld.unique) ### which Loci identified as being in LD

  ## turn the linked loci into a dataframe, then make sure the SNP names are characters so can be searched against the LD and Fst vectors
      Linked.df <- Linked
      Linked.df$SNP_A <- as.character(Linked.df$SNP_A)
      Linked.df$SNP_B <- as.character(Linked.df$SNP_B)

  ### keep only linked loci < the size of the panel you wish to create
      loci.in.LD.vec <- loci.in.LD[which(loci.in.LD < panel.size)] ## which loci in LD < the size of the panel

  # what rows of the ld data frame contain the the linked loci in the top n of Fst values?
      what.positions.get <- list()
      for(k in 1:length(loci.in.LD.vec)){

        in.ld.min <- loci.in.LD.vec[k]

        get.rows <- which(Linked.df$SNP_A == FST.order.vec[in.ld.min] | Linked.df$SNP_B == FST.order.vec[in.ld.min])
        # print(get.rows)
        #what.rows.get <- list(what.rows.get, get.rows)
        what.positions.get[[k]] <- get.rows
      }

  ## Get the names of the SNPs in LD from the PLINK LD file
      SNP.out <- list()
      for(i in 1:length(what.positions.get)){
        to.get <- what.positions.get[[i]]
        SNP_loop <- NULL
        for(k in 1:length(what.positions.get[[i]])){
          SNP.hold.loop <- c(Linked.df[what.positions.get[[i]][k], "SNP_A"], Linked.df[what.positions.get[[i]][k], "SNP_B"]) ## get the values from teh row defined for SNP_A and B
          SNP_loop <- rbind(SNP_loop, SNP.hold.loop) ## there may be more than one row defined, rbind them together
        }
        SNP.out[[i]] <- SNP_loop ## add to new list
      }

  ## compare the SNPs in the LD file to the ranked order by Fst - save the ranks of the SNPs that are linked - but only unique values, i.e. if
  ## one SNP is linked to two other SNPs, only record 3 names, not 4. Ranks = positions in vector of ranked loci

    linked.ranks <- list()
    for(j in 1:length(SNP.out)){

      to.double.check <- SNP.out[[j]] ## gets names of the jth SNPs and what other SNP they are linked to
      dbl.chk <- to.double.check[which(stringr::str_detect(string = noquote(to.double.check), pattern = as.character(FST.order.vec[loci.in.LD.vec[j]]))==FALSE)]
      ### the jth SNPs are saved as a string in to.double.check, this removes the jth SNP in loci.fst from this, so only have non-duplicated

      where.best.link.fst <- loci.in.LD.vec[j] ## the loci with the greatest Fst among the linked ones, will be the jth in the ranked vector
      where.its.linked <- which(FST.order.vec %in% dbl.chk) ## get the loctions on hte ranked vector of the loci it is linked to
      hold.linked.ranks <- c(where.best.link.fst, where.its.linked) ## cbind them together
      linked.ranks[[j]] <- hold.linked.ranks ## output in a list

    }

  ## turn the list into a dataframe by adding NA where there are too few 'columns' in a string
    linked.ranks.df <- plyr::rbind.fill(lapply(linked.ranks,function(y){as.data.frame(t(y),stringsAsFactors=FALSE)}))

  ### now this is where it get a list of numbers to be removed from the ranked vector of loci names by Fst
    h.rows <- which(linked.ranks.df$V1<panel.size) ##
    to.cut.out <- NULL

    for(i in 1:length(h.rows)){
      a <- h.rows[i]

      ## this ensures that where we get more than one of the linked loci in the top n loci, we keep the one with the highest Fst (lowest number in the rank)
      if((linked.ranks.df$V1[a] < linked.ranks.df$V2[a])==TRUE){
        to.cut <- linked.ranks.df$V2[a]
        to.cut.out <- c(to.cut.out, to.cut)
      }
    }

    #clean up temporary files used in the analysis.
    file.remove(remember.spidpath)
    file.remove(sub_data_path)
    file.remove(ped.path)
    file.remove(map.path)
    file.remove(paste0(where.PGDspider, "/hyb.spid"))
    file.remove(paste0(where.PGDspider, "/GP_FSTAT.spid"))
    file.remove(fst_data_path)
    file.remove(plink_map_path)
    file.remove(plink_ped_path)
    file.remove(paste0(path.start, "/plink.txt"))
    file.remove(paste0(path.start, "/LDsReform.txt"))
    file.remove(remember.TOPLOC.path)
    file.remove(paste0(where.PGDspider,"/subset_for_LD.txt"))
    file.remove(paste0(where.PGDspider,"/for_FST.txt"))
    file.remove(paste0(path.start,"/GP_FSTAT.spid"))
    #file.remove(GenePop)

    #wrap up indicator
    writeLines("Process Completed.")

  ## return loci ordered by fst
      your.panel <- data.frame(loci=FST.order.vec[-to.cut.out][1:panel.size],stringsAsFactors = F)
      your.panel <- merge(your.panel,FST.df,by="loci")
      return(your.panel[order(your.panel$FSTs,decreasing = TRUE),])

}
