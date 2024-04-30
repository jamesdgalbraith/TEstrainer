library(tidyverse)
library(plyranges)
library(Biostrings)

dfam_seq_tbl <- tibble(qseqid = names(Biostrings::readDNAStringSet("raw/Dfam_curatedonly.fasta")),
                       qlen = width(Biostrings::readDNAStringSet("raw/Dfam_curatedonly.fasta"))) %>%
  tidyr::separate(col = qseqid, into = c("seqnames", "superfamily"), sep = "#") %>%
  tidyr::separate(col = superfamily, into = c("subclass", "superfamily"), sep = "/")

dfam_rps <- read_tsv("raw/Dfam_curatedonly_rps.out", show_col_types = F) %>%
  dplyr::filter(length/slen > 0.5) %>%
  dplyr::mutate(start = ifelse(qstart < qend, qstart, qend),
                end = ifelse(qend > qstart, qend, qstart),
                strand = ifelse(qstart < qend, "+", "-")) %>%
  dplyr::select(-qstart, -qend) %>%
  tidyr::separate(col = stitle, into = c("ref", "abbrev", "full"), sep = ", ") %>%
  tidyr::separate(col = qseqid, into = c("seqnames", "superfamily"), sep = "#") %>%
  tidyr::separate(col = superfamily, into = c("subclass", "superfamily"), sep = "/") %>%
  mutate(superfamily = ifelse(grepl("ERV", superfamily), "ERV", superfamily)) %>%
  dplyr::mutate(superfamily = case_when(
    grepl("SINE", superfamily) ~ "SINE",
    grepl("/", superfamily) ~ sub(".*/", "", superfamily),
    TRUE ~ superfamily)
  ) %>%
  mutate(
    superfamily = ifelse(superfamily %in% c("TcMar", "TcMar-ISRm11", "TcMar-m44", "TcMar-Pogo", "TcMar-Fot1", "TcMar-Tc2", "TcMar-Mariner", "TcMar-Tc1", "TcMar-Tc4", "TcMar-Tigger"), "Mariner/Tc1", superfamily),
    superfamily = ifelse(superfamily %in% c("hAT-hAT5", "hAT-hAT19", "hAT-hATx", "hAT-hATm", "hAT-Tag1", "hAT-Ac", "hAT-Blackjack", "hAT-Tip100", "hAT-Charlie"), "hAT", superfamily),
    superfamily = ifelse(superfamily %in% c("Pao"), "BEL", superfamily),
    superfamily = ifelse(superfamily %in% c("CMC-Transib"), "Transib", superfamily),
    superfamily = ifelse(superfamily %in% c("PiggyBac"), "piggyBac", superfamily),
    superfamily = ifelse(superfamily %in% c("Naiad", "Chlamys"), "Penelope", superfamily),
    superfamily = ifelse(superfamily %in% c("CMC-Chapaev", "CMC-Chapaev-3", "CMC-EnSpm"), "EnSpm/CACTA", superfamily),
    superfamily = ifelse(superfamily %in% c("Academ-1"), "Academ", superfamily),
    superfamily = ifelse(superfamily %in% c("Ngaro"), "DIRS", superfamily),
    superfamily = ifelse(superfamily %in% c("PIF-ISL2EU"), "ISL2EU", superfamily),
    superfamily = ifelse(superfamily %in% c("R2-Hero"), "Hero", superfamily),
    superfamily = ifelse(superfamily %in% c("CR1-Zenon"), "CR1", superfamily),
    superfamily = ifelse(superfamily %in% c("endogenous retrovirus", "Endogenous Retrovirus", "ERV-1_PM-I", "ERV-2_PM-I"), "ERV", superfamily),
    superfamily = ifelse(superfamily %in% c("LTR retrotransposon", "LTR Retrotransposon"), "LTR", superfamily),
    superfamily = ifelse(superfamily %in% c("RTE-RTE", "RTE-BovB"), "RTE", superfamily)
  )

all_domain_info <- dfam_rps %>%
  dplyr::select(abbrev, ref, full) %>%
  base::unique()

DNA <- dfam_rps %>% filter(subclass == "DNA")
LINE <- dfam_rps %>% filter(subclass == "LINE")
PLE <- dfam_rps %>% filter(subclass == "PLE")
LTR <- dfam_rps %>% filter(subclass == "LTR")
RC <- dfam_rps %>% filter(subclass %in% c("RC"))
other <- dfam_rps %>% filter(!subclass %in% c("LTR", "DNA", "LINE", "PLE", "RC"))

# ID acceptable domains based on having 5 or 10 or more TEs with this domain
LTR_good_domains <- as_tibble(as.data.frame(table(LTR$abbrev))) %>%
  arrange(-Freq) %>%
  filter(Freq >= 10) %>%
  mutate(abbrev = as.character(Var1))
DNA_good_domains <- as_tibble(as.data.frame(table(DNA$abbrev))) %>%
  arrange(-Freq) %>%
  filter(Freq >= 5) %>% # Due to lower number of certain subclasses in Dfam lower threshold
  mutate(abbrev = as.character(Var1))
PLE_good_domains <- as_tibble(as.data.frame(table(PLE$abbrev))) %>%
  arrange(-Freq) %>%
  filter(Freq >= 10) %>%
  mutate(abbrev = as.character(Var1))
RC_good_domains <- as_tibble(as.data.frame(table(RC$abbrev))) %>%
  arrange(-Freq) %>%
  filter(Freq >= 5)  %>% # Due to less Helitrons in Dfam lower threshold
  mutate(abbrev = as.character(Var1))
LINE_good_domains <- as_tibble(as.data.frame(table(LINE$abbrev))) %>%
  arrange(-Freq) %>%
  filter(Freq >= 10) %>%
  mutate(abbrev = as.character(Var1))

# Identify good domains within subclasses, pull out all domains from them for overlap analysis
LINE_good <- LINE[LINE$abbrev %in% LINE_good_domains$abbrev,]
LINE_good <- LINE[LINE$seqnames %in% LINE_good$seqnames,] %>%
  mutate(good_domain = abbrev %in% LINE_good_domains$abbrev)

PLE_good <- PLE[PLE$abbrev %in% PLE_good_domains$abbrev,]
PLE_good <- PLE[PLE$seqnames %in% PLE_good$seqnames,] %>%
  mutate(good_domain = abbrev %in% PLE_good_domains$abbrev)

LTR_good <- LTR[LTR$abbrev %in% LTR_good_domains$abbrev,]
LTR_good <- LTR[LTR$seqnames %in% LTR_good$seqnames,] %>%
  mutate(good_domain = abbrev %in% LTR_good_domains$abbrev)

DNA_good <- DNA[DNA$abbrev %in% DNA_good_domains$abbrev,]
DNA_good <- DNA[DNA$seqnames %in% DNA_good$seqnames,] %>%
  mutate(good_domain = abbrev %in% DNA_good_domains$abbrev)

RC_good <- RC[RC$abbrev %in% RC_good_domains$abbrev,]
RC_good <- RC[RC$seqnames %in% RC_good$seqnames,] %>%
  mutate(good_domain = abbrev %in% RC_good_domains$abbrev)

## Intersect to find good domains overlapping low copy good domains
# LINE
LINE_good_ranges <- LINE_good %>%
  filter(good_domain == TRUE) %>%
  dplyr::select(-subclass, -superfamily, -qlen, -good_domain, -length, -evalue, -bitscore) %>%
  as_granges()

LINE_potential_ranges <- LINE_good %>%
  filter(good_domain == FALSE) %>%
  dplyr::select(-subclass, -superfamily, -qlen, -good_domain, -length, -evalue, -bitscore) %>%
  as_granges()

LINE_likely_joined <- plyranges::join_overlap_inner_directed(LINE_good_ranges, LINE_potential_ranges, minoverlap = 15) %>%
  as_tibble()

LINE_additional_good <- LINE_likely_joined %>%
  dplyr::select(ref.y, abbrev.y, full.y) %>%
  dplyr::rename(ref = ref.y, abbrev = abbrev.y, full = full.y) %>%
  base::unique()

# PLE
PLE_good_ranges <- PLE_good %>%
  filter(good_domain == TRUE) %>%
  dplyr::select(-subclass, -superfamily, -qlen, -good_domain, -length, -evalue, -bitscore) %>%
  as_granges()

PLE_potential_ranges <- PLE_good %>%
  filter(good_domain == FALSE) %>%
  dplyr::select(-subclass, -superfamily, -qlen, -good_domain, -length, -evalue, -bitscore) %>%
  as_granges()

PLE_likely_joined <- plyranges::join_overlap_inner_directed(PLE_good_ranges, PLE_potential_ranges, minoverlap = 15) %>%
  as_tibble()

PLE_additional_good <- PLE_likely_joined %>%
  dplyr::select(ref.y, abbrev.y, full.y) %>%
  dplyr::rename(ref = ref.y, abbrev = abbrev.y, full = full.y) %>%
  base::unique()

# LTR
LTR_good_ranges <- LTR_good %>%
  filter(good_domain == TRUE) %>%
  dplyr::select(-subclass, -superfamily, -qlen, -good_domain, -length, -evalue, -bitscore) %>%
  as_granges()

LTR_potential_ranges <- LTR_good %>%
  filter(good_domain == FALSE) %>%
  dplyr::select(-subclass, -superfamily, -qlen, -good_domain, -length, -evalue, -bitscore) %>%
  as_granges()

LTR_likely_joined <- plyranges::join_overlap_inner_directed(LTR_good_ranges, LTR_potential_ranges, minoverlap = 15) %>%
  as_tibble()

LTR_additional_good <- LTR_likely_joined %>%
  dplyr::select(ref.y, abbrev.y, full.y) %>%
  dplyr::rename(ref = ref.y, abbrev = abbrev.y, full = full.y) %>%
  base::unique()

# DNA
DNA_good_ranges <- DNA_good %>%
  filter(good_domain == TRUE) %>%
  dplyr::select(-subclass, -superfamily, -qlen, -good_domain, -length, -evalue, -bitscore) %>%
  as_granges()

DNA_potential_ranges <- DNA_good %>%
  filter(good_domain == FALSE) %>%
  dplyr::select(-subclass, -superfamily, -qlen, -good_domain, -length, -evalue, -bitscore) %>%
  as_granges()

DNA_likely_joined <- plyranges::join_overlap_inner_directed(DNA_good_ranges, DNA_potential_ranges, minoverlap = 15) %>%
  as_tibble()

DNA_additional_good <- DNA_likely_joined %>%
  dplyr::select(ref.y, abbrev.y, full.y) %>%
  dplyr::rename(ref = ref.y, abbrev = abbrev.y, full = full.y) %>%
  base::unique()

# RC
RC_good_ranges <- RC_good %>%
  filter(good_domain == TRUE) %>%
  dplyr::select(-subclass, -superfamily, -qlen, -good_domain, -length, -evalue, -bitscore) %>%
  as_granges()

RC_potential_ranges <- RC_good %>%
  filter(good_domain == FALSE) %>%
  dplyr::select(-subclass, -superfamily, -qlen, -good_domain, -length, -evalue, -bitscore) %>%
  as_granges()

RC_likely_joined <- plyranges::join_overlap_inner_directed(RC_good_ranges, RC_potential_ranges, minoverlap = 15) %>%
  as_tibble()

RC_additional_good <- RC_likely_joined %>%
  dplyr::select(ref.y, abbrev.y, full.y) %>%
  dplyr::rename(ref = ref.y, abbrev = abbrev.y, full = full.y) %>%
  base::unique()

## Compile data and write to file
compiled_LINE_good <- rbind(all_domain_info[all_domain_info$abbrev %in% LINE_good_domains$abbrev,], LINE_additional_good) %>% mutate(subclass="LINE")
compiled_PLE_good <- rbind(all_domain_info[all_domain_info$abbrev %in% PLE_good_domains$abbrev,], PLE_additional_good) %>% mutate(subclass="PLE")
compiled_LTR_good <- rbind(all_domain_info[all_domain_info$abbrev %in% LTR_good_domains$abbrev,], LTR_additional_good) %>% mutate(subclass="LTR")
compiled_DNA_good <- rbind(all_domain_info[all_domain_info$abbrev %in% DNA_good_domains$abbrev,], DNA_additional_good) %>% mutate(subclass="DNA")
compiled_RC_good <- rbind(all_domain_info[all_domain_info$abbrev %in% RC_good_domains$abbrev,], RC_additional_good) %>% mutate(subclass="RC")
all_compiled_good <- rbind(compiled_LINE_good, compiled_PLE_good, compiled_LTR_good, compiled_DNA_good, compiled_RC_good)
write_tsv(x = all_compiled_good, file = "data/acceptable_domains_2.tsv")

