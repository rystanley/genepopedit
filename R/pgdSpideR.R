# pgdSpideR
#' @title Execute pgdSpider conversions in R
#' @description Function to between file types in R.
#' @param input the full path to the input file.
#' @param input_type what are you converting from (e.g. GENEPOP or PED) these must match the dropdown menus in pgdSpider.
#' @param output the full file path to the desired output. (e.g. GENEPOP or PED) these must match the dropdown menus in pgdSpider.
#' @param output_type what are you converting to (e.g. GENEPOP or PED). This must match the dropdown menus in pgdSpider.
#' @param spid the fill file path to the spid file generated in pgdSpider.
#' @param where.PGDspider A file path to the PGDspider installation folder.
#' @rdname pgdSpideR
#' @importFrom stringr str_split
#' @export

pgdSpideR <- function(input,input_format,output,output_format,spid,where.PGDspider){

where.PGDspider.PGD <- gsub(x = where.PGDspider, pattern = " ", replacement = "\\ ", fixed = TRUE)

### Make sure the spid file required for the conversion is within the executable directory
file.copy(from = spid, to = where.PGDspider, overwrite = TRUE)

## move the input file as well to the same location as PGDspider - this makes this step so much easier
file.copy(from <- input, to = where.PGDspider, overwrite = TRUE)

#File names
input.name <- stringr::str_split(string = input, pattern = "/")
input.name <- unlist(input.name)
input.name <- input.name[grep(x = input.name, pattern = ".txt")]

spid.name <- stringr::str_split(string = spid, pattern = "/")
spid.name <- unlist(spid.name)
spid.name <- spid.name[grep(x = spid.name, pattern = ".spid")]

output.name <- stringr::str_split(string = output, pattern = "/")
output.name <- unlist(output.name)
output.name <- output.name[length(output.name)]

### create a string to call PGDspider
input.file.call <- paste("-inputfile", input.name, sep = " ")
execute.SPIDER <- "java -Xmx1024m -Xms512m -jar PGDSpider2-cli.jar"
spid.call <- paste("-spid",spid.name,sep=" ")
input.format <- paste("-inputformat",input_format,sep=" ")
output.format <- paste("-outputformat",output_format,sep=" ")

## check the operating system
if(Sys.info()["sysname"]=="Windows"){
goto.spider <- paste0("cd ", where.PGDspider.PGD, " && ", execute.SPIDER)
output.file.path <- paste("-outputfile", output.name, sep = " ")

## string to run
run.PGDspider <- paste0(goto.spider, " ", input.file.call, " ",
                        input.format, " ", output.file.path, " ",
                        output.format, " ", spid.call)

### run PGDspider through system
shell(run.PGDspider)} else {

    goto.spider <- paste0("cd ", where.PGDspider.PGD, "; ", execute.SPIDER)
    output.file.path <- paste("-outputfile", output.name, sep = " ")

    ## string to run
    run.PGDspider <- paste0(goto.spider, " ", input.file.call, " ",
                            input.format, " ", output.file.path, " ",
                            output.format, " ", spid.call)

    ### run PGDspider through system
    system(run.PGDspider)

}


## move the converted file to the output pathPLINK folder
file.copy(from = paste0(where.PGDspider,output.name), to = output, overwrite = TRUE)

## remove the files you placed in the pgdSpider folder
file.remove(paste0(where.PGDspider,input.name))
file.remove(paste0(where.PGDspider,output.name))
file.remove(paste0(where.PGDspider,spid.name))

writeLines("Process Completed.")

 }
