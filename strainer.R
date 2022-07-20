#!/usr/bin/Rscript

library(optparse)
library(tidyverse)
library(plyranges)
library(BSgenome)

opt <- list(in_seq = "seq/cleaned_Echis_RepeatModeler.fasta", run_trf = FALSE, run_rps = TRUE)
opt$seq <- sub(".*/", "", opt$in_seq)

# read in fasta
rm_seq_in <- readDNAStringSet(paste0("seq/", opt$seq))
rm_seq_info <- tibble(seqnames = sub(" .*", "", names(rm_seq_in)), width = width(rm_seq_in))

if(isTRUE(opt$run_rps)){
  acceptable_domains <- read_tsv("acceptable_domains.tsv")
  
  # read rps blast out
  rps_blast_out <- read_tsv(file = paste0("data/", opt$seq, ".rps.out"),
                            col_names = c("seqnames", "qstart", "qend", "qlen", "sseqid", "sstart", "send", "slen",
                                          "pident", "length", "mismatch", "gapopen", "evalue", "bitscore", "qcovs", "stitle")) %>%
    separate(stitle, into = c("ref", "abbrev", "full"), sep = ", ") %>%
    dplyr::filter(length >= 0.5*slen)
  
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
  no_domains <- rm_seq_in[!names(rm_seq_in) %in% c(chimeric$seqnames, completely_acceptable$seqnames, questionable$seqnames)]
  
}

if(isTRUE(opt$run_trf)){
  # tandem repeat filtering (optional)
  trf_out_ranges <- read_gff3(paste0("data/", opt$seq, ".trf.gff"))
  
  trf_out_tbl <- as_tibble(trf_out_ranges) %>%
    mutate(copies = as.double(copies),
           period = as.double(period),
           consensus_size = as.double(consensus_size)) %>%
    dplyr::filter(copies >=2.1) %>%
    as_granges() %>%
    reduce_ranges() %>%
    as_tibble() %>%
    mutate(seqnames = as.character(seqnames),
           length = end - start + 1) %>%
    dplyr::group_by(seqnames) %>%
    dplyr::mutate(total_length = sum(length)) %>%
    dplyr::select(-width) %>%
    inner_join(rm_seq_info)
  
  trf_out_tbl[trf_out_tbl$total_length/trf_out_tbl$width >= 0.25,]
  trf_out_ranges[seqnames(trf_out_ranges) == "Echis_carinatus_rnd-1_family-119#DNA/TcMar-Tc2"]
  
}

completely_acceptable_seq <- rm_seq_in[names(rm_seq_in) %in% compiled_acceptable$seqnames]
no_domains_seq <- no_domains
writeXStringSet(c(completely_acceptable_seq, no_domains_seq), paste0("out/cleaned_", opt$seq))
chimeric_seq <- rm_seq_in[names(rm_seq_in) %in% seqnames(truly_chimeric_ranges)]
writeXStringSet(chimeric_seq, paste0("out/chimeric_", opt$seq))
questionable_seq <- rm_seq_in[names(rm_seq_in) %in% questionable$seqnames]
writeXStringSet(questionable_seq, paste0("out/questionable_", opt$seq))
