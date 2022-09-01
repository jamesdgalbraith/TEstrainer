#!/usr/bin/Rscript

library(optparse)

option_list <- list(
  make_option(c("-i", "--in_seq"), default=NA, type = "character", help="Path to input sequence (required)"),
  make_option(c("-o", "--out"), default=NA, type = "character", help="Path to output folder (required)"),
  make_option(c("-b", "--rps_table"), default=NA, type = "character", help="Path rpstblastn table")
)

opt <- list(in_seq = "TEstrainer_143116_01Sep/dfam_lepidosaurs.fasta", out = "TEstrainer_143116_01Sep/chimeras/", rps_table = "TEstrainer_143116_01Sep/chimeras/dfam_lepidosaurs.fasta.rps.out")

# opt <- parse_args(OptionParser(option_list=option_list))
opt$seq <- sub(".*/", "", opt$in_seq)
# check variables provided
if(is.na(opt$in_seq)){
  stop("Path to in sequence must be supplied")
}
if(is.na(opt$out)){
  stop("Path to in sequence must be supplied")
}
if(is.na(opt$rps_table)){
  stop("Path to rpstblastn table must be supplied")
}

# make empty variable function
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(BSgenome))

# read in fasta
rm_seq_in <- readDNAStringSet(paste0(opt$in_seq))
names(rm_seq_in) <- sub(" .*", "", names(rm_seq_in))
rm_seq_info <- tibble(seqnames = names(rm_seq_in), width = width(rm_seq_in))

compiled_acceptable <- tibble()
truly_chimeric_ranges <- GRanges()
questionable <- tibble()
no_domains_seq <- DNAStringSet()

# read in acceptable domains, rbind additionally discovered
additional_domains <- read_tsv("data/additional_domains.tsv", show_col_types = FALSE) %>%
  dplyr::select(ref, abbrev)
acceptable_domains <- read_tsv("data/acceptable_domains.tsv", show_col_types = FALSE) %>%
  rbind(additional_domains)

as_tibble(as.data.frame(table(additional_domains$ref))) %>% arrange(-Freq)

if(file.size(opt$rps_table)==0){
  writeXStringSet(rm_seq_in, paste0(opt$out, "/clean_", opt$seq))
  quit()
}

# read rps blast out
rps_blast_out <- read_tsv(file = opt$rps_table,
                                           col_names = c("seqnames", "qstart", "qend", "qlen", "sseqid", "sstart", "send", "slen",
                                                         "pident", "length", "mismatch", "gapopen", "evalue", "bitscore", "qcovs", "stitle"),
                                           show_col_types = FALSE) %>%
                   tidyr::separate(stitle, into = c("ref", "abbrev", "full"), sep = ", ")

rps_blast_out <- rps_blast_out %>% dplyr::filter(length >= 0.5*slen)

# select acceptable and non acceptable
contains_acceptable <- rps_blast_out[rps_blast_out$ref %in% acceptable_domains$ref,]
contains_unacceptable <- rps_blast_out[!rps_blast_out$ref %in% acceptable_domains$ref,]

# Filter repeats containing only suitable domains
completely_acceptable <- rps_blast_out %>%
  filter(seqnames %in% contains_acceptable$seqnames,
         !seqnames %in% contains_unacceptable$seqnames) %>%
  arrange(seqnames, qstart)

# Filter repeats containing no suitable domains
questionable <- rps_blast_out %>%
  filter(!seqnames %in% contains_acceptable$seqnames,
         seqnames %in% contains_unacceptable$seqnames) %>%
  arrange(seqnames, qstart)

# Filter repeats containing both suitable and unsuitable domains
chimeric <- rps_blast_out %>%
  filter(seqnames %in% contains_acceptable$seqnames,
         seqnames %in% contains_unacceptable$seqnames) %>%
  arrange(seqnames, qstart) %>%
  mutate(unacceptable = ifelse(ref %in% acceptable_domains$ref, "false", "true"))

# Potentially mask out genic regions of chimeric repeats
chimeric_ranges <- chimeric %>%
  dplyr::mutate(start = ifelse(qstart < qend, qstart, qend),
                end = ifelse(qstart > qend, qstart, qend),
                strand = ifelse(qstart < qend, "+", "-")) %>%
  dplyr::select(seqnames, start, end, strand, abbrev, unacceptable) %>%
  plyranges::as_granges()

chimeric_ranges_unacceptable_false <- chimeric_ranges[chimeric_ranges$unacceptable == "false",]
chimeric_ranges_unacceptable_true <- chimeric_ranges[chimeric_ranges$unacceptable == "true",]

# find regions which don't overlap with acceptable domains
truly_chimeric_ranges <- filter_by_non_overlaps(chimeric_ranges_unacceptable_true, chimeric_ranges_unacceptable_false)

truly_chimeric <- chimeric[chimeric$seqnames %in% seqnames(truly_chimeric_ranges),]
false_positive_chimeric <- chimeric[!chimeric$seqnames %in% seqnames(truly_chimeric_ranges),] %>%
  dplyr::select(-unacceptable)
compiled_acceptable <- rbind(completely_acceptable, false_positive_chimeric)

# identify sequences with no domains
no_domains_seq <- rm_seq_in[!names(rm_seq_in) %in% c(chimeric$seqnames, completely_acceptable$seqnames, questionable$seqnames)]

## add step to combine data
# write to file (check if any filtered, if not write all in to output)
completely_acceptable_seq <- rm_seq_in[names(rm_seq_in) %in% compiled_acceptable$seqnames]
writeXStringSet(c(completely_acceptable_seq, no_domains_seq), paste0(opt$out, "/clean_", opt$seq))
chimeric_seq <- rm_seq_in[names(rm_seq_in) %in% seqnames(truly_chimeric_ranges)]
writeXStringSet(chimeric_seq, paste0(opt$out, "/chimeric_", opt$seq))
questionable_seq <- rm_seq_in[names(rm_seq_in) %in% questionable$seqnames]
writeXStringSet(questionable_seq, paste0(opt$out, "/questionable_", opt$seq))

truly_chimeric


false_positive_chimeric_normal_ranges <- false_positive_chimeric[false_positive_chimeric$ref %in% acceptable_domains$ref,] %>%
  mutate(start = ifelse(qstart < qend, qstart, qend),
         end = ifelse(qstart > qend, qstart, qend),
         strand = ifelse(qstart < qend, "+", "-")) %>%
  dplyr::select(seqnames, start, end, strand, ref, abbrev) %>%
  as_granges() %>%
  reduce_ranges_directed(ref = ref, abbrev = abbrev)

false_positive_chimeric_odd_ranges <- false_positive_chimeric[!false_positive_chimeric$ref %in% acceptable_domains$ref,] %>%
  mutate(start = ifelse(qstart < qend, qstart, qend),
         end = ifelse(qstart > qend, qstart, qend),
         strand = ifelse(qstart < qend, "+", "-")) %>%
  dplyr::select(seqnames, start, end, strand, ref, abbrev) %>%
  as_granges() %>%
  reduce_ranges_directed(ref = ref, abbrev = abbrev)

overlapping <- join_overlap_intersect(false_positive_chimeric_normal_ranges, false_positive_chimeric_odd_ranges) %>%
  as_tibble() %>%
  dplyr::select(-ref.x) %>%
  mutate(ref.y = sub(" ", "", sub("\\)", "", sub("c\\(", "", gsub('\\\"', '', as.character(ref.y))))),
         abbrev.x = sub(" ", "", sub("\\)", "", sub("c\\(", "", gsub('\\\"', '', as.character(abbrev.x)))))) %>%
  arrange(abbrev.y, abbrev.x)

overlapping %>%
  View()
