#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})


# parse input variables
option_list = list(
  make_option(c("-g", "--genome"), type="character", default=NULL,
              help="Path to genome", metavar="character"),
  make_option(c("-l", "--library"), type="character", default=NULL,
              help="Path to library", metavar="character"),
  make_option(c("-n", "--iteration"), default=NULL, type = "integer",
              help="Iteration number (required)"),
  make_option(c("-f", "--flank"), type="integer", default=1500,
              help="Flank length", metavar="integer")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
# opt <- list(genome = "seq/Echis_ycKpl.FINAL.fasta", library = "data/run_1/Echis_RepeatModeler.fasta", iteration = 1, flank = 1500)
if (is.null(opt$genome)) {
  stop("Path to genome is needed")
} 
if (is.null(opt$iteration)){
  stop("Interation number needed")
}

suppressPackageStartupMessages({
  library(BSgenome)
  library(plyranges)
  library(tidyverse)
})

message("Reading blast")
blast_out <- read_tsv(file = paste0("data/run_", opt$iteration, "/",  opt$library,"_initial_blast.out"),
                      col_names = c("qseqid", "seqnames", "pident", "length", "qstart", "qend",
                                    "qlen", "sstart", "send", "slen", "evalue", "bitscore", "qcovs"),
                      show_col_types = F)

blast_out_trimmed <- blast_out %>%
  dplyr::filter(
    pident >= 70,
    qcovs >= 50
    ) %>%
  dplyr::group_by(qseqid) %>%
  dplyr::arrange(-bitscore) %>%
  dplyr::slice(1:20) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(strand = base::ifelse(sstart < send, "+", "-"),
         start = base::ifelse(sstart < send, sstart - opt$flank, send - opt$flank),
         end = base::ifelse(sstart > send, sstart + opt$flank, send + opt$flank),
         start = ifelse(start < 1, 1, start),
         end = ifelse(end > slen, slen, end),
         seqnames = paste0(qseqid, "___", seqnames)) %>%
  dplyr::select(seqnames, qseqid, start, end, strand) %>%
  plyranges::as_granges()

blast_out_merged <- as_tibble(plyranges::reduce_ranges(blast_out_trimmed)) %>%
  mutate(seqnames = as.character(seqnames)) %>%
  separate(seqnames, into = c("qseqid", "seqnames"), sep = "___") %>%
  dplyr::select(-strand, -width) %>%
  plyranges::as_granges()

blast_out_tbl <- tibble(qseqid = names(table(blast_out_merged$qseqid)), n = as.integer(table(blast_out_merged$qseqid))) %>%
  filter(n > 2)

genome_seq <- Biostrings::readDNAStringSet(filepath = opt$genome)
names(genome_seq) <- sub(" .*", "", names(genome_seq))

consensus_seq <- Biostrings::readDNAStringSet(filepath = paste0("data/run_", opt$iteration, "/", opt$library))
names(consensus_seq) <- sub(" .*", "", names(consensus_seq))

if(!dir.exists(paste0("data/run_", opt$iteration, "/initial_seq/"))){
  dir.create(paste0("data/run_", opt$iteration, "/initial_seq/"))
}

# create unaligned multifasta
for(i in seq_along(blast_out_tbl$qseqid)){

  align_ranges <- blast_out_merged[blast_out_merged$qseqid == blast_out_tbl$qseqid[i]]
  align_ranges <- GenomicRanges::reduce(align_ranges, ignore.strand=T)
  align_seq <- BSgenome::getSeq(genome_seq, align_ranges)
  names(align_seq) <- base::paste0(GenomeInfoDb::seqnames(align_ranges), ":", IRanges::ranges(align_ranges))
  
  align_seq <- c(consensus_seq[names(consensus_seq) == blast_out_tbl$qseqid[i]], align_seq)
  
  Biostrings::writeXStringSet(x = align_seq,
                              filepath = paste0("data/run_", opt$iteration, "/initial_seq/",
                                                sub("#.*", "", blast_out_tbl$qseqid[i]), ".fasta"))

}
