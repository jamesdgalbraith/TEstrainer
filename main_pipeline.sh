#!/bin/bash

# make directories
mkdir -p data/run_0/raw out

# cluster
cd-hit-est -n 10 -c 0.95 -i ${RM_LIBRARY} -o data/run_0/${RM_LIBRARY} # cluster seq

# split
PIECES="$(grep '>' data/run_0/${RM_LIBRARY})"
Rscript scripts/splitter.R -t nt -f data/run_0/${RM_LIBRARY} -o data/run_0/raw out/ -p $PIECES # split query
ls data/run_0/raw/${RM_LIBRARY}* | sed 's/.*\///' > data/run_0/to_run.txt # list queries

if [ ! -f "${GENOME}".nsq ]; then
makeblastdb -in ${GENOME} -dbtype nucl -out ${GENOME} # makeblastb if needed
fi

# prestrain
mkdir -p data/run_0/rps_out/ data/run_0/trf_out/
parallel --bar --jobs $threads -a data/run_0/to_run.txt 'rpstblastn -query data/run_0/raw/{} -db /media/projectDrive_1/databases/cdd/Cdd -out data/run_0/rps_out/{}.out -outfmt "6 qseqid qstart qend qlen sseqid sstart send slen pident length mismatch gapopen evalue bitscore qcovs stitle" -evalue 0.01 -num_threads 1'
cat data/run_0/rps_out/{}.out > data/run_0/${RM_LIBRARY}.rps.out
# parallel --bar --jobs $threads -a data/run_0/to_run.txt 'trf data/run_0/trf_out/{} 2 7 7 80 10 50 500 -d -h'
# parallel --bar --jobs $threads -a data/run_0/to_run.txt 'python scripts/trf2gff.py -d {}.2.7.7.80.10.50.500.dat -o data/run_0/trf_out/{}.trf.gff'
# cat data/run_0/raw/{}*.trf.gff > data/run_0/${RM_LIBRARY}.trf.gff
Rscript strainer.R --in_seq data/run_0/${RM_LIBRARY} --run_n 0


# runs

  # extend/align
  

  # trim
  
  # strain
  mkdir -p data/run_0/rps_out/ data/run_0/trf_out/
  parallel --bar --jobs $threads -a data/run_0/to_run.txt 'rpstblastn -query data/run_0/raw/{} -db /media/projectDrive_1/databases/cdd/Cdd -out data/run_0/rps_out/{}.out -outfmt "6 qseqid qstart qend qlen sseqid sstart send slen pident length mismatch gapopen evalue bitscore qcovs stitle" -evalue 0.01 -num_threads 1'
  cat data/run_0/rps_out/{}.out > data/run_0/${RM_LIBRARY}.rps.out
  # parallel --bar --jobs $threads -a data/run_0/to_run.txt 'trf data/run_0/trf_out/{} 2 7 7 80 10 50 500 -d -h'
  # parallel --bar --jobs $threads -a data/run_0/to_run.txt 'python scripts/trf2gff.py -d {}.2.7.7.80.10.50.500.dat -o data/run_0/trf_out/{}.trf.gff'
  # cat data/run_0/raw/{}*.trf.gff > data/run_0/${RM_LIBRARY}.trf.gff
  Rscript strainer.R --in_seq data/run_0/${RM_LIBRARY} --run_n 0