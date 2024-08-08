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

all_compiled_good <- read_tsv(file = "data/acceptable_domains_2.tsv")

# tm7 <- dfam_rps %>%
#   filter(startsWith(abbrev, prefix = "7tm"))
# tm7 <- dfam_rps %>%
#   filter(seqnames %in% tm7$seqnames)
# View(tm7)
# 
# # dfam_rps %>%
#   filter(grepl("vomeronasal", full, ignore.case=TRUE)) %>%
#   dplyr::select(subclass, abbrev) %>%
#   base::unique() %>%
#   View()
# 
# dfam_rps %>%
#   filter(grepl("globulin", full, ignore.case=TRUE)) %>%
#   dplyr::select(subclass, abbrev) %>%
#   base::unique() %>%
#   View()
# 
# dfam_rps %>%
#   filter(grepl("sugar", full, ignore.case=TRUE)) %>%
#   dplyr::select(subclass, abbrev) %>%
#   base::unique() %>%
#   View()
