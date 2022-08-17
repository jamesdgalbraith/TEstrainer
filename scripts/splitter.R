#!/usr/bin/env Rscript

library("optparse")

# parse input variables
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-p", "--pieces"), type="integer", default=parallel::detectCores(all.tests = FALSE, logical = TRUE),
              help="number of pieces to break file into (default is number of threads)", metavar="integer"),
  make_option(c("-t", "--type"), type="character", default=NULL, 
              help="sequence type (DNA or AA)", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="split/", 
              help="path to output [default= %default]", metavar="character"),
  make_option(c("-r", "--rename"), action="store_true", default=TRUE,
              help="Name output files based on file number"),
  make_option(c("-n", "--name"), action="store_false", 
              dest="rename", help="Name output files after first sequence")
  
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# check input sequence is provided
if (is.null(opt$file)) {
  stop("Input sequence is required. See script usage (--help)")
}

# check sequence type is provided
if (opt$type %in% c("dna", "DNA")) {
  compiled_seq <- Biostrings::readDNAStringSet(opt$file)
} else if (opt$type %in% c("aa", "AA", "nt", "NT")) {
  compiled_seq <- Biostrings::readAAStringSet(opt$file)
} else {
  stop("sequence type is needed")
}

# if outdir doesn't exist create
if(!dir.exists(opt$out)){dir.create(opt$out)}

# Count number of sequences
no_seq <- length(compiled_seq)

# Check number of output files is less than or equal desired umber of pieces
if(no_seq < opt$pieces){opt$pieces <- no_seq}

# Split and write to file
for( i in 1:opt$pieces){
  
  out_seq <- compiled_seq[(ceiling((i-1)*no_seq/opt$pieces)+1):(ceiling((i)*no_seq/opt$pieces))]
  if(isTRUE(opt$rename)){
    Biostrings::writeXStringSet(x = out_seq, filepath = paste0(opt$out, "/", sub(".*/", "", opt$file),"_seq_", i, ".fasta"))
  } else {
    Biostrings::writeXStringSet(x = out_seq, filepath = paste0(opt$out, "/", sub("#.*", "", names(out_seq)[1]), ".fasta"))
  }
  
}