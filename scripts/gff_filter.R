library(plyranges)

in_path <- "data/CHYMOMYZA_COSTATA.filteredRepeats.gff"
out_path <- paste(tools::file_path_sans_ext(in_path), "strained", tools::file_ext(in_path), sep = ".")
library_path <- 'data/CHYMOMYZA_COSTATA/CHYMOMYZA_COSTATA-families.fa.strained.dirty'

unacceptable_families <- tolower(sub("#.*", "", names(Biostrings::readDNAStringSet(library_path))))
read_gff(in_path) %>%
  filter(!tolower(ID) %in% unacceptable_families) %>%
  write_gff2(out_path)


                                                