#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

# parse input variables
option_list = list(
  make_option(c("-l", "--library"), type="character", default=NULL, 
              help="library file name", metavar="character"),
  make_option(c("-d", "--directory"), type="character", default=NULL, 
              help="path to working data directory", metavar="character")
  
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# check input sequence is provided
if (is.null(opt$library)) {
  stop("Input library is required.")
}
# check input sequence is provided
if (is.null(opt$directory)) {
  stop("Output directory is required.")
}

suppressPackageStartupMessages(library(BSgenome))

repeat_library <- suppressMessages(readDNAStringSet(opt$library))
suppressWarnings(repeat_library <- repeat_library[!is.na(as.integer(sub("#.*", "", sub("DR", "", names(repeat_library))))) & startsWith(names(repeat_library), "DR"),])
writeXStringSet(repeat_library, paste0(opt$directory, "/", sub(".*/", "", opt$library)))
