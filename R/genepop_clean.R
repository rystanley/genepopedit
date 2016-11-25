# genepop_clean
#' @title Fixes common formatting problems using PGDspider.
#' @description Use PGDspider to impliment the proper delimiters in your genepop file.
#' @param genepop A file path to the genepop format file you wish to create your panel from
#' @param where.pgdspider A file path to the PGDspider installation folder.
#' @param allocate.pgd.ram An integer value in GB to specify the maximum amount of RAM to allocate to PGDspider. The default is 1 GB, which should be sufficient for most analyses.
#' @param sample_delim A logical variable (default = TRUE) that defines whether the function needs to fix the sampleID loci delimiter using PGDspider. Here all files require (space,space scpace) to delimit the sampleIDs from the loci
#' @param population_suffix A logical variable (default = FALSE) that defines whether the function needs to add a population specific suffix to each population. This is required for some genepopedit functionality. Between each "pop" label a Pop(n)_ suffix will be added to each sampleID of genepop.
#' @param pop_sample_delim A logical variable (default = FALSE) that defines whether the function needs to add separation ("_") between population name and sample number. Note if population suffix = TRUE this option is redundant and not executed.
#' @param path the filepath and filename of output.#' @rdname genepop_toploci
#' @importFrom data.table fread as.data.table
#' @importFrom stringr str_split str_detect
#' @rdname genepop_clean
#' @export


genepop_clean <- function(genepop,where.pgdspider,allocate.pgd.ram = 1,sample_delim=TRUE,population_suffix=FALSE,pop_sample_delim=FALSE,path){

  path.start <- getwd()  ### where to write the files created by genepopedit to


  if(sample_delim){
        #Set up ram allocation. Note that only 1024 gb of ram will work with Windows
        allocate.pgd.ram <- allocate.pgd.ram*1024

        if(Sys.info()["sysname"] == "Windows" & allocate.pgd.ram>1024){
          allocate.pgd.ram=1024
          writeLines("Note that currently PGDspider can only utilize ~1 GB of ram on windows based operating systems. Periodically check back to https://github.com/rystanley/genepopedit for any updates to this limitation.
                     ")}


        ## Fix the delimiter spacing the sampleIDs and the loci using PGD spider

        #Spid file that defines how PGDspider will convert the file.
        cleanspid <- "# spid-file generated: Fri Nov 25 10:34:45 AST 1901

        # GENEPOP Parser questions
        PARSER_FORMAT=GENEPOP

        # Enter the size of the repeated motif (same for all loci: one number; different: comma separated list (e.g.: 2,2,3,2):
        GENEPOP_PARSER_REPEAT_SIZE_QUESTION=
        # Select the type of the data:
        GENEPOP_PARSER_DATA_TYPE_QUESTION=SNP
        # How are Microsat alleles coded?
        GENEPOP_PARSER_MICROSAT_CODING_QUESTION=REPEATS

        # GENEPOP Writer questions
        WRITER_FORMAT=GENEPOP

        # Specify which data type should be included in the GENEPOP file  (GENEPOP can only analyze one data type per file):
        GENEPOP_WRITER_DATA_TYPE_QUESTION=SNP
        # Specify the locus/locus combination you want to write to the GENEPOP file:
        GENEPOP_WRITER_LOCUS_COMBINATION_QUESTION="

        path.start.PGD <- gsub(x = path.start, pattern = " ", replacement = "\\")
        where.pgdspider.PGD <- gsub(x = where.pgdspider, pattern = " ", replacement = "\\ ", fixed = TRUE)

        #write conversion parameter file (".spid")
        write(x = cleanspid, file = paste0(path.start, "/", "GP_GP.spid"))

        ### move spid file to the PGDspider folder
        file.copy(from = paste0(path.start, "/GP_GP.spid"), to = where.pgdspider, overwrite = TRUE)
        remember.spidpath <- paste0(path.start, "/", "GP_GP.spid")
        ## move the input file as well to the same location as PGDspider - this makes this step so much easier
        file.copy(from <- genepop, to = where.pgdspider, overwrite = TRUE)
        GenePop.name <- stringr::str_split(string = genepop, pattern = "/")
        GenePop.name <- unlist(GenePop.name)
        GenePop.name <- GenePop.name[grep(x = GenePop.name, pattern = ".txt")]
        file.rename(from = paste0(where.pgdspider, GenePop.name), to = paste0(where.pgdspider, "GPD_for_clean.txt"))
        remember.TOPLOC.path <- paste0(where.pgdspider, "GPD_for_clean.txt")


        #Console message
        writeLines("Cleaning GENEPOP sampleID and loci delimiters.")
        writeLines("

                   ")
        writeLines("Warning messages are expected as part of conversion process using PGDspider.

                   ")

        #### OS X and LINUX CALL

        if(Sys.info()["sysname"] != "Windows"){

          ### create a string to call PGDspider
          input.file.call <- "-inputfile GPD_for_clean.txt"
          execute.SPIDER <- paste0("java -Xmx", allocate.pgd.ram, "m -Xms512m -jar PGDSpider2-cli.jar")
          spid.call <- "-spid GP_GP.spid"
          input.format <- "-inputformat GENEPOP"
          output.format <- "-outputformat GENEPOP"
          goto.spider <- paste0("cd ", where.pgdspider.PGD, "; ", execute.SPIDER)
          output.file.path <- "-outputfile for_clean.txt"
          ## string to run
          run.PGDspider <- paste0(goto.spider, " ", input.file.call, " ", input.format, " ", output.file.path, " ", output.format, " ", spid.call)
          ### run PGDspider through system
          system(run.PGDspider)

        } # End MAC LINUX IF

        #### Windows call

        if(Sys.info()["sysname"] == "Windows"){

          ### create a string to call PGDspider
          input.file.call <- "-inputfile GPD_for_clean.txt"
          execute.SPIDER <- paste0("java -Xmx", allocate.pgd.ram, "m -Xms512m -jar PGDSpider2-cli.jar")
          spid.call <- "-spid GP_GP.spid"
          input.format <- "-inputformat GENEPOP"
          output.format <- "-outputformat GENEPOP"
          goto.spider <- paste0("cd ", where.pgdspider.PGD, " && ", execute.SPIDER)
          output.file.path <- "-outputfile for_clean.txt"
          ## string to run
          run.PGDspider <- paste0(goto.spider, " ", input.file.call, " ", input.format, " ", output.file.path, " ", output.format, " ", spid.call)

          ### run PGDspider through system
          shell(run.PGDspider)

        } # End WINDOWS IF

        ### move the cleaned GENEPOP file back to the working directory
        file.copy(from = paste0(where.pgdspider, "/for_clean.txt"), to = path, overwrite = TRUE)
        ## remember the path of the file created by genepop_clean
        fst_data_path <- paste0(path.start, "/", "for_clean.txt")

        #clean up temporary files used in the analysis.
        file.remove(paste0(path.start, "/", "GP_GP.spid"))

  }


  ### if there is no distinction in the sample IDs ###
  if(population_suffix){

        #Console message
        writeLines("Cleaning GENEPOP sampleIDs and adding population labels.")
        writeLines("

                   ")
        writeLines("Note population labels are genetic (Pop(n)_). If specific names are required please refer to subset_genepop_rename().

                   ")

        SampleIDs <- genepop_detective(path,variable="Inds")

        #Check to see if genepop is a file path or dataframe
        if(is.character(path)){
          genepop <- data.table::fread(path,
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
          genepop <- data.table::as.data.table(genepop,stringsAsFactors = FALSE)
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

        #Evaluate the length of each population
        PopLengths <- NULL
        for(i in 2:length(tempPops)){
          PopLengths <- c(PopLengths,tempPops[i]-tempPops[i-1]-1)
        }
        PopLengths <- c(PopLengths,nrow(temp)-sum(PopLengths))#add the length of the last population

        ## Get the population names (prior to the _ in the Sample ID)
        NamePops <- paste0(rep(paste0("Pop",1:length(PopLengths),"_"),times=PopLengths),temp[,1])
        NameExtract <- substr(NamePops,1,regexpr("_",NamePops)-1)

        #Now stitch the data together
        # paste together the Loci as one long integer separated for each loci by a space
        Loci <- do.call(paste,c(temp2[,], sep=" "))

        #Paste these to the Loci
        Loci <- paste(NamePops,Loci,sep=" ,  ")

        PopLengths <- table(factor(NameExtract, levels=unique(NameExtract)))[-length(table(NameExtract))]

        if(length(table(NameExtract))==2){PopPosition = PopLengths+1}

        if(length(table(NameExtract))>2){
          PopPosition <- c(PopLengths[1]+1,rep(NA,(length(PopLengths)-1)))
          for (i in 2:length(PopLengths)){
            PopPosition[i] <- PopLengths[i]+PopPosition[i-1]
          }
        }

        #Insert the value of "Pop" which partitions the data among populations #only if more than one population
        if(length(table(NameExtract))!=1){Loci <- insert_vals(Vec=Loci,breaks=PopPosition,newVal="Pop")}

        #Add the first "Pop" label
        Loci <- c("Pop",Loci)

        Output <- c(stacks.version,names(temp2),Loci)

        # Save the file
        utils::write.table(Output,path,col.names=FALSE,row.names=FALSE,quote=FALSE)

  }



  #if have population specific suffix values but they are not delimited by a _ between the population label and the sample number.

  if(!population_suffix & pop_sample_delim){
        genepop_ID(genepop = path,path = path)
  }


}




