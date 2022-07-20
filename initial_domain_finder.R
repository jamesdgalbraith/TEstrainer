#!/usr/bin/Rscript

library(ORFik)
library(tidyverse)
library(plyranges)
library(optparse)
library(BSgenome)

opt <- list()

# read in repbase sequences
repbase_dna <- Biostrings::readDNAStringSet("/media/projectDrive_1/databases/repbase/RepBase5May2021.fasta")

# get info about repbase sequences
repbase_seq_info <- tibble(seqnames = names(repbase_dna), width = width(repbase_dna))
repbase_seq_info <- repbase_seq_info %>%
  tidyr::separate(col = seqnames, into = c("seqnames", "class", "species"), sep = "\t") %>%
  dplyr::mutate(name = 1:nrow(repbase_seq_info))

# get orf info to ensure only TEs with ORFs over 300aa included (removes noise from non-autonomous)
# positive strands
pos <- findORFs(repbase_dna, startCodon = "ATG", minimumLength = 300)
# negative strands
neg <- findORFs(reverseComplement(repbase_dna),
                startCodon = "ATG", minimumLength = 300)
pos <- GRanges(pos, strand = "+")
neg <- GRanges(neg, strand = "-")
res <- c(pos, neg) %>%
  as_tibble() %>%
  mutate(name = as.integer(as.character(seqnames))) %>%
  dplyr::select(-seqnames, -width)

repbase_seq_info <- inner_join(repbase_seq_info, res) %>%
  dplyr::select(seqnames, class, species, width) %>%
  base::unique()

# read in repbase rpstblastn
repbase_rps_out <- readr::read_tsv(file = "data/RepBase5May2021.fasta.rps.out",
                                   col_names = c("seqnames", "qstart", "qend", "qlen", "sseqid", "sstart", "send", "slen",
                                                 "pident", "length", "mismatch", "gapopen", "evalue", "bitscore", "qcovs", "stitle"),
                                   show_col_types = F)
repbase_rps_out <- tidyr::separate(data = repbase_rps_out, col = stitle, into = c("ref", "abbrev", "full"), sep = ", ")
repbase_rps_out <- dplyr::inner_join(repbase_rps_out, repbase_seq_info)
domain_info <- repbase_rps_out %>% dplyr::select(ref, abbrev, full) %>% base::unique()

# get info about repbase classes
repbase_class_info <- dplyr::as_tibble(BiocGenerics::as.data.frame(table(repbase_seq_info$class)))
repbase_class_info <- dplyr::mutate(.data = repbase_class_info, Var1 = as.character(Var1))
repbase_class_info <- dplyr::rename(.data = repbase_class_info, class = Var1)
repbase_class_info_min <- repbase_class_info[repbase_class_info$Freq < 10,]
repbase_class_info_max <- repbase_class_info[repbase_class_info$Freq >= 10,]
repbase_class_info_max <- dplyr::arrange(.data = repbase_class_info_max, -Freq)

# loop over to get info
common_domains <- tibble()
for(i in 1:nrow(repbase_class_info_max)){
  holder <- repbase_rps_out[repbase_rps_out$class == repbase_class_info_max$class[i],] %>%
    dplyr::filter(length >= 0.5*slen)
  if(nrow(holder >= 1)){
    holder_tbl <- as_tibble(as.data.frame(table(holder$ref))) %>%
      mutate(Var1 = as.character(Var1),
             perc = Freq/repbase_class_info_max$Freq[i]) %>%
      dplyr::rename(ref = Var1) %>%
      filter(perc >= 0.05, Freq >= 5) %>%
      inner_join(domain_info) %>%
      mutate(class = repbase_class_info_max$class[i]) %>%
      dplyr::arrange(-Freq)
    common_domains <- rbind(common_domains, holder_tbl)
    }
}

# remove multicopy domains, write to file
base::unique(common_domains %>% select(ref, abbrev)) %>%
  dplyr::filter(!startsWith(tolower(abbrev), "ig"),
                !startsWith(tolower(abbrev), "znmc"),
                !startsWith(tolower(abbrev), "7tm")) %>%
  write_tsv("acceptable_domains.tsv")