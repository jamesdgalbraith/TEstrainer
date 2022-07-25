#!/bin/bash

# parsing
RM_LIBRARY=$(echo $RM_LIBRARY_PATH | sed 's/.*\///')
echo $GENOME $RM_LIBRARY $THREADS $FLANK

# make directories
mkdir -p data/run_0/raw out/

# cluster
cd-hit-est -n 10 -c 0.95 -i ${RM_LIBRARY_PATH} -o data/run_0/${RM_LIBRARY} # cluster seq

# split
PIECES="$(grep -c '>' data/run_0/${RM_LIBRARY})"
Rscript scripts/splitter.R -t nt -f data/run_0/${RM_LIBRARY} -o data/run_0/raw out/ -p $PIECES # split query
ls data/run_0/raw/${RM_LIBRARY}* | sed 's/.*\///' > data/run_0/to_run.txt # list queries

if [ ! -f "${GENOME}".nsq ]; then
makeblastdb -in ${GENOME} -dbtype nucl -out ${GENOME} # makeblastb if needed
fi

# prestrain
mkdir -p data/run_0/rps_out/ data/run_0/trf_out/
parallel --bar --jobs $THREADS -a data/run_0/to_run.txt 'rpstblastn -query data/run_0/raw/{} -db /media/projectDrive_1/databases/cdd/Cdd -out data/run_0/rps_out/{}.out -outfmt "6 qseqid qstart qend qlen sseqid sstart send slen pident length mismatch gapopen evalue bitscore qcovs stitle" -evalue 0.01 -num_threads 1'
cat data/run_0/rps_out/{}.out > data/run_0/${RM_LIBRARY}.rps.out
# parallel --bar --jobs $threads -a data/run_0/to_run.txt 'trf data/run_0/trf_out/{} 2 7 7 80 10 50 500 -d -h'
# parallel --bar --jobs $threads -a data/run_0/to_run.txt 'python scripts/trf2gff.py -d {}.2.7.7.80.10.50.500.dat -o data/run_0/trf_out/{}.trf.gff'
# cat data/run_0/raw/{}*.trf.gff > data/run_0/${RM_LIBRARY}.trf.gff
Rscript strainer.R --in_seq data/run_0/${RM_LIBRARY} --run_n 0


# runs
  # RUN_NO=1; while  [ $RUN_NO -le $RUNS ]; do echo "Welcome $x times" $(( RUN_NO++ )); done
  # split
  mkdir -p data/run_${RUN_NO}/raw
  cp data/run_$(expr $RUN_NO - 1)/clean_${RM_LIBRARY} data/run_${RUN_NO}/${RM_LIBRARY}
  PIECES=$(grep -c '>' data/run_${RUN_NO}/${RM_LIBRARY})
  Rscript scripts/splitter.R -t nt -f data/run_${RUN_NO}/${RM_LIBRARY} -o data/run_${RUN_NO}/raw -p $PIECES
  ls data/run_${RUN_NO}/raw/* | sed 's/.*\///' > data/run_${RUN_NO}/to_run.txt
  
  
  # extend/align
  ## initial blast
  mkdir -p data/run_${RUN_NO}/initial_blast
  parallel --env --bar --jobs ${THREADS} -a data/run_${RUN_NO}/to_run.txt blastn -task dc-megablast -query data/run_${RUN_NO}/raw/{} -db $GENOME -evalue 1e-5 -outfmt \"6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore\" -out data/run_${RUN_NO}/initial_blast/{}.out -num_threads 1 # search genome  
  cat data/run_${RUN_NO}/initial_blast/*.out > data/run_${RUN_NO}/initial_blast.out # compile data
  
  ## extention and check blast
  Rscript scripts/self_blast_setup.R -g ${GENOME} -l ${RM_LIBRARY} -n ${RUN_NO} # extend seqs
  ls data/run_${RUN_NO}/initial_seq/*fasta | sed 's/.*\///' > data/run_${RUN_NO}/self_queries.txt
  mkdir data/run_${RUN_NO}/self_search/
  parallel --bar --jobs ${THREADS} -a data/run_${RUN_NO}/self_queries.txt blastn -task dc-megablast -query data/run_${RUN_NO}/initial_seq/{} -subject data/run_${RUN_NO}/initial_seq/{} -evalue 1e-5 -outfmt \"6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore\" -out data/run_${RUN_NO}/self_search/{}.out -num_threads 1 # self blast
  
  ## clip seqs and align
  mkdir -p data/run_${RUN_NO}/to_align data/run_${RUN_NO}/mafft
  Rscript scripts/mafft_setup.R -g ${GENOME} -n ${RUN_NO} -l ${RM_LIBRARY} # trim seqs pre-mafft
  
  ## align seqs
  if [[ $THREADS -gt 4 ]]; then; MAFFT_THREADS=$(($(($THREADS / 4)))); else; MAFFT_THREADS=1; fi
  parallel --bar --jobs $MAFFT_THREADS -a data/run_${RUN_NO}/to_align.txt "mafft --thread 4 --quiet --localpair --adjustdirectionaccurately data/run_${RUN_NO}/to_align/{} > data/run_${RUN_NO}/mafft/{}"
  ls data/run_${RUN_NO}/mafft/ > data/run_${RUN_NO}/to_trim.txt
  
  # trim
  mkdir -p data/run_${RUN_NO}/TEtrim_con \
           data/run_${RUN_NO}/TEtrim_complete  \
           data/run_${RUN_NO}/TEtrim_unaln \
           data/run_${RUN_NO}/TEtrim_blast \
           data/run_${RUN_NO}/TEtrim_mafft \
           data/run_${RUN_NO}/TEtrim_further \
           data/run_${RUN_NO}/TEtrim_bp

  parallel --bar --jobs 12 -a data/run_${RUN_NO}/to_trim.txt scripts/TEtrim.py --i {} --threads 4 --flank 1500
  
  
  # strain
  
           
  parallel --bar --jobs $threads -a data/run_${RUN_NO}/to_run.txt 'rpstblastn -query data/run_0/raw/{} -db /media/projectDrive_1/databases/cdd/Cdd -out data/run_0/rps_out/{}.out -outfmt "6 qseqid qstart qend qlen sseqid sstart send slen pident length mismatch gapopen evalue bitscore qcovs stitle" -evalue 0.01 -num_threads 1'
  cat data/run_0/rps_out/${RM_LIBRARY}*.out > data/run_0/${RM_LIBRARY}.rps.out
  # parallel --bar --jobs $threads -a data/run_0/to_run.txt 'trf data/run_0/trf_out/{} 2 7 7 80 10 50 500 -d -h'
  # parallel --bar --jobs $threads -a data/run_0/to_run.txt 'python scripts/trf2gff.py -d {}.2.7.7.80.10.50.500.dat -o data/run_0/trf_out/{}.trf.gff'
  # cat data/run_0/raw/{}*.trf.gff > data/run_0/${RM_LIBRARY}.trf.gff
  Rscript strainer.R --in_seq data/run_${RUN_NO}/trimmed_${RM_LIBRARY} --iteration ${RUN_NO}
  