#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})


# parse input variables
option_list <- list(
  make_option(c("-g", "--genome"), type="character", default=NULL,
              help="genome name", metavar="character"),
  make_option(c("-l", "--library"), type="character", default=NULL,
              help="Path to library", metavar="character"),
  make_option(c("-n", "--iteration"), default=NULL, type = "integer",
              help="Iteration number (required)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$genome)) {
  stop("Genome name is needed")
} else {
  # set genome names
  opt$genome_name <- sub(".*/", "", opt$genome)
}

suppressPackageStartupMessages({
  library(BSgenome)
  library(plyranges)
  library(tidyverse)
})

# determine queries

queries <- read_tsv(paste0("data/run_", opt$iteration, "/self_queries.txt"), col_names = "qseqid", show_col_types = F)

# make placeholders for files to/not to align
to_align <- tibble(query = character())
not_align <- tibble(query = character())

for (i in seq_along(queries$qseqid)) {
  
  in_seq <- Biostrings::readDNAStringSet(filepath = paste0("data/run_", opt$iteration, "/initial_seq/", queries$qseqid[i])) # read in sequence
  
  self_out <- read_tsv(file = paste0("data/run_", opt$iteration, "/self_search/", queries$qseqid[i], ".out"), # read in blast data
                       col_names = c("seqnames", "sseqid", "pident", "length", "qstart", "qend",
                                     "qlen", "sstart", "send", "slen", "evalue", "bitscore"), show_col_types = F) %>%
    filter(sseqid != seqnames) %>% # remove self hits
    filter(!grepl(".*_family-.*#", seqnames)) %>% # remove original consensus (RM)
    filter(sub("#.*", "", seqnames) != queries$qseqid[i]) %>%
    filter(length > 0.5 * BiocGenerics::width(in_seq[1])) # ensure alignment at least 80% of original consensus length

  suppressWarnings(self_ranges <- self_out %>%
    group_by(seqnames) %>%
    filter(qstart >= quantile(qstart, 0.1), # ensure start and end are within 0.1 and 0.9 quantiles
           qend <= quantile(qend, 0.9)) %>%
    dplyr::arrange(-bitscore) %>%
    dplyr::slice(1:5) %>%
    mutate(start = min(qstart), # select extreme start+end
           end = max(qend)) %>%
    dplyr::ungroup() %>%
    dplyr::select(seqnames, start, end) %>%
    base::unique() %>%
    dplyr::arrange(seqnames, start, end) %>%
    plyranges::as_granges()) # convert to granges

  out_seq <- getSeq(in_seq, self_ranges) # get seq
  names(out_seq) <- seqnames(self_ranges) # give names
  
  out_seq <- c(in_seq[1], out_seq) # add original sequence
  
  writeXStringSet(out_seq, paste0("data/run_", opt$iteration, "/to_align/", queries$qseqid[i])) # write to file
  
  to_align <- rbind(to_align, tibble(query = queries$qseqid[i]))

}

# write list of sequences to be aligned to file
write_tsv(to_align, paste0("data/run_", opt$iteration, "/to_align.txt"), col_names = F)

# create file of unextendable repeats
consensus_seq <- Biostrings::readDNAStringSet(filepath = paste0("data/run_", opt$iteration, "/", opt$library))
names(consensus_seq) <- sub(" .*", "", names(consensus_seq))
writeXStringSet(consensus_seq[!sub("#.*", "", names(consensus_seq)) %in% sub(".fasta", "", to_align$query)],
                filepath = paste0("data/run_", opt$iteration, "/rare_", opt$library))
