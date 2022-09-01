#!/usr/bin/Rscript

library(optparse)

option_list <- list(
  make_option(c("-i", "--in_seq"), default=NA, type = "character", help="Input sequence (required)"),
  make_option(c("-d", "--directory"), type="character", default=NULL, help="Path to data directory (required)", metavar="character"),
  make_option(c("-s", "--sort"), type="logical", default=FALSE, help="Set to sort final output")
)

opt <- parse_args(OptionParser(option_list=option_list))

# check variables provided
if(is.na(opt$in_seq)){
  stop("Path to in sequence must be supplied")
}
if(is.na(opt$directory)){
  stop("Path to in sequence must be supplied")
}
opt <- list(in_seq = "seq/dfam_lepidosaurs.fasta", directory = "TEstrainer_143116_01Sep/")
opt$out_seq <- sub(".*/", "", opt$in_seq)

# make empty variable function
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(BSgenome))

in_seq <- readDNAStringSet(opt$in_seq)
in_seq_tbl <- tibble(seqnames = names(in_seq), og_width = width(in_seq)) %>%
  dplyr::mutate(draft_seqnames = sub("#.*", "", seqnames))

# read in and rearrange SA-SSR data
sassr <- read_tsv(paste0(opt$directory, "/trf/", opt$out_seq, ".sassr"),
                  col_names = c("seqnames", "ssr", "count", "start"), skip = 1, show_col_types = F) %>%
  dplyr::mutate(period = as.double(width(ssr))) %>%
  dplyr::mutate(ssr_width = count*period, end = start + ssr_width, start = start +1) %>%
  inner_join(in_seq_tbl, by = "seqnames") %>%
  arrange(seqnames) %>%
  filter(count > 2)

sassr_calc <- as_tibble(reduce(as_granges(sassr))) %>%
  dplyr::select(-strand) %>%
  dplyr::mutate(seqnames = as.character(seqnames)) %>%
  group_by(seqnames) %>%
  mutate(total_width = sum(width)) %>%
  ungroup() %>%
  inner_join(in_seq_tbl, by = "seqnames") %>%
  dplyr::mutate(sassr_perc_tr = 100*total_width/og_width) %>%
  dplyr::select(seqnames, sassr_perc_tr) %>%
  base::unique()

# read in and rearrange TRF data
trf <- read_tsv(paste0(opt$directory, "/trf/", opt$out_seq, ".trf"),
                col_names = c("draft_seqnames", "start", "end", "period", "count", "ssr"), show_col_types = F) %>%
  dplyr::mutate(ssr_width = end - start + 1) %>%
  inner_join(in_seq_tbl, by = "draft_seqnames")

trf_select <- trf  %>%
  filter(count > 2) %>%
  dplyr::select(seqnames, start, end)

trf_calc <- as_tibble(reduce(as_granges(trf_select))) %>%
  dplyr::select(-strand) %>%
  dplyr::mutate(seqnames = as.character(seqnames)) %>%
  group_by(seqnames) %>%
  mutate(total_width = sum(width)) %>%
  ungroup() %>%
  inner_join(in_seq_tbl, by = "seqnames") %>%
  dplyr::mutate(trf_perc_tr = 100*total_width/og_width) %>%
  dplyr::select(seqnames, trf_perc_tr) %>%
  base::unique()

# read in and rearrange MREPS data
mreps <- read_tsv(paste0(opt$directory, "/trf/", opt$out_seq, ".mreps"),
                  col_names = c("draft_seqnames", "start", "end", "ssr_width", "period", "count", "error", "sequence"), show_col_types = F) %>%
  mutate(ssr = substr(x = sequence, start = 0, stop = period)) %>%
  inner_join(in_seq_tbl, by = "draft_seqnames")

mreps_select <- mreps %>%
  filter(count > 2) %>%
  dplyr::select(seqnames, start, end)

mreps_calc <- as_tibble(reduce(as_granges(mreps_select))) %>%
  dplyr::select(-strand) %>%
  dplyr::mutate(seqnames = as.character(seqnames)) %>%
  group_by(seqnames) %>%
  mutate(total_width = sum(width)) %>%
  ungroup() %>%
  inner_join(in_seq_tbl, by = "seqnames") %>%
  dplyr::mutate(mreps_perc_tr = 100*total_width/og_width) %>%
  dplyr::select(seqnames, mreps_perc_tr) %>%
  base::unique()

# Compile data
trf_select <- trf %>% select(seqnames, start, end, period, count, og_width, ssr, ssr_width) %>%
  mutate(package = "trf")
mreps_select <- mreps %>% select(seqnames, start, end, period, count, og_width, ssr, ssr_width) %>%
  mutate(package = "mreps")
sassr_select <- sassr %>% select(seqnames, start, end, period, count, og_width, ssr, ssr_width) %>%
  mutate(package = "sassr")
compiled_tr <- rbind(rbind(trf_select, mreps_select), sassr_select)
stats_tr <- full_join(sassr_calc, full_join(mreps_calc, trf_calc, by = "seqnames"), by = "seqnames")

# Identifying satellite/simple repeats
over50tr <- stats_tr %>%
  filter(sassr_perc_tr >= 50 | mreps_perc_tr >= 50 | trf_perc_tr >= 50) %>%
  mutate(sassr_perc_tr = ifelse(is.na(sassr_perc_tr), 0, sassr_perc_tr),
         mreps_perc_tr = ifelse(is.na(mreps_perc_tr), 0, mreps_perc_tr),
         trf_perc_tr = ifelse(is.na(trf_perc_tr), 0, trf_perc_tr))
satellites <- compiled_tr[compiled_tr$seqnames %in% over50tr$seqnames,] %>%
  dplyr::group_by(seqnames) %>%
  dplyr::arrange(-ssr_width, period) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(perc_tr = 100*ssr_width/og_width) %>%
  inner_join(over50tr, by = "seqnames") %>%
  dplyr::mutate(max_perc_tr = pmax(sassr_perc_tr, mreps_perc_tr, trf_perc_tr)) %>%
  dplyr::arrange(max_perc_tr)

# for macrosatellites trim to get single copy
macrosatellites <- satellites %>%
  filter(max_perc_tr > 90) %>%
  filter(period > 200) %>%
  mutate(start = start + period,
         end = start + period) %>%
  as_granges()
macrosatellites_seq <- getSeq(in_seq, macrosatellites)
names(macrosatellites_seq) <- seqnames(macrosatellites)
writeXStringSet(macrosatellites_seq, paste0(opt$directory, "/trf/", opt$out_seq, ".macrosatellites"))

# leave micro to minisatellite repeats intact
other_satellites <- satellites %>%
  dplyr::filter(!seqnames %in% seqnames(macrosatellites)) %>%
  dplyr::mutate(start = 1, end = og_width) %>%
  as_granges()
other_satellites_seq <- getSeq(in_seq, other_satellites)
names(other_satellites_seq) <- seqnames(other_satellites)
writeXStringSet(other_satellites_seq, paste0(opt$directory, "/trf/", opt$out_seq, ".satellites"))

# filtering for trimming
trim_3 <- compiled_tr %>%
  mutate(tr_len = period*count) %>%
  filter(count >=3, !seqnames %in% over50tr$seqnames, period < 20) %>%
  filter(og_width-end <= period) %>%
  dplyr::group_by(seqnames) %>%
  arrange(-end, period) %>%
  dplyr::slice(1)
trim_5 <- compiled_tr %>%
  mutate(tr_len = period*count) %>%
  filter(count >=3, !seqnames %in% over50tr$seqnames, period < 20) %>%
  filter(start <= period) %>%
  dplyr::group_by(seqnames) %>%
  arrange(-end, period) %>%
  dplyr::slice(1)

# Determine thos to be trimmed both ends
trim_both_5 <- trim_5[trim_5$seqnames %in% trim_3$seqnames,] %>%
  dplyr::mutate(start = ifelse(period < 4, end - 8, end - (2*period))) %>%
  dplyr::mutate(start = ifelse(start < 1, 1, start)) %>%
  dplyr::select(seqnames, start)
trim_both_3 <- trim_3[trim_3$seqnames %in% trim_5$seqnames,] %>%
  dplyr::mutate(end = ifelse(period < 4, start + 8, start + (2*period))) %>%
  dplyr::mutate(end = ifelse(end > og_width, og_width, end)) %>%
  dplyr::select(seqnames, end)
trim_both <- inner_join(trim_both_5, trim_both_3, by = "seqnames")

# determine those to be trimmed one end
trim_only_5 <- trim_5[!trim_5$seqnames %in% trim_3$seqnames,] %>%
  dplyr::mutate(start = ifelse(period < 4, end - 8, end - (2*period))) %>%
  dplyr::mutate(start = ifelse(start < 1, 1, start), end = og_width) %>%
  dplyr::select(seqnames, start, end)
trim_only_3 <- trim_3[!trim_3$seqnames %in% trim_5$seqnames,] %>%
  dplyr::mutate(end = ifelse(period < 4, start + 8, start + (2*period))) %>%
  dplyr::mutate(end = ifelse(end > og_width, og_width, end), start = 1) %>%
  dplyr::select(seqnames, start, end)

to_trim <- as_granges(rbind(trim_both, rbind(trim_only_5, trim_only_3)))
trimmed_seq <- getSeq(in_seq, to_trim)
names(trimmed_seq) <- seqnames(to_trim)
writeXStringSet(trimmed_seq, paste0(opt$directory, "/trf/", opt$out_seq, ".trimmed"))

# determine untouched sequences
untouched_seq <- in_seq[!names(in_seq) %in% c(names(trimmed_seq), names(macrosatellites_seq), names(other_satellites_seq))]
writeXStringSet(untouched_seq, paste0(opt$directory, "/trf/", opt$out_seq, ".nonsatellite"))

# put it all together
compiled_fixed_seq <- c(macrosatellites_seq, other_satellites_seq, trimmed_seq, untouched_seq)

if(opt$sort == TRUE){
  
  compiled_fixed_sorted <- tibble(seqnames = names(compiled_fixed_seq), start = 1, end = width(compiled_fixed_seq),
         numbering = sub(".*rnd-", "", sub("#.*", "", names(compiled_fixed_seq)))) %>%
    mutate(round = as.integer(sub("_.*", "", numbering)), family = as.integer(sub(".*-", "", numbering))) %>%
    arrange(round, family, seqnames) %>%
    as_granges()
  compiled_fixed_sorted_seq <- getSeq(compiled_fixed_seq, compiled_fixed_sorted)
  names(compiled_fixed_sorted_seq) <- seqnames(compiled_fixed_sorted)
  
  # write all corrected sequences
  writeXStringSet(compiled_fixed_sorted_seq, paste0(opt$directory, "/trf/trimmed_", opt$out_seq))
  
} else {

  writeXStringSet(compiled_fixed_seq, paste0(opt$directory, "/trf/trimmed_", opt$out_seq))  
  
}

