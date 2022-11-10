# TEstrainer
A pipeline to dramatically improve repeat libraries before genome annotation through curation of each repeat, better identification of satellite repeats, removal of multi-copy genes and reclassification of cleaned libraries.


## Required packages and databases

BLAST+ - Altschul, S.F. et al. (1997). https://doi.org/10.1089/10665270050081478.
MAFFT - Katoh, K and Standley, D.M. (2013). https://doi.org/10.1093/molbev/mst010
mreps - Kolpakov, R et al. (2003). https://doi.org/10.1093/nar/gkg617
sa-ssr - Pickett, B.D. (2016). https://doi.org/10.1093/bioinformatics/btw298
GNU Parallel - Tange, O. (2011) The USENIX Magazine, February 2011:42-47.
cd-hit (Optional, used for clustering) - Li, W (2006) https://doi.org/10.1093/bioinformatics/btl158
RepeatModeler (Optional, used for reclassifying) - Flynn, J.M. (2020). https://doi.org/10.1073/pnas.1921046117

### Python3 packages
pyranges - Stovner, E.B and Sætrom, P. (2020). https://doi.org/10.1093/bioinformatics/btz615
pyfaidx - Shirley MD, et al. (2015). https://doi.org/10.7287/peerj.preprints.970v1
biopython - Cock, P.A. et al. (2009). https://doi.org/10.1093/bioinformatics/btp163
pandas - McKinney, W. (2010). https://doi.org/10.25080/Majora-92bf1922-00a 
numpy - Harris, C.R. et al. (2020). https://doi.org/10.1038/s41586-020-2649-2. 

### R packages
optparse
tidyverse -  Wickham, H. et al. (2019). https://doi.org/doi:10.21105/joss.01686. 
BSgenome - Pagès H (2022). https://bioconductor.org/packages/BSgenome.
plyranges - Lee, S. et al. (2019). https://doi.org/10.1186/s13059-018-1597-8.

### Database
NCBI CDD - Lu, S. et al (2020). https://doi.org/10.1093/nar/gkz991 