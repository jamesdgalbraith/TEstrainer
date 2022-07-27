#!/bin/bash

# parsing
RM_LIBRARY=$(echo $RM_LIBRARY_PATH | sed 's/.*\///')
# echo $GENOME $RM_LIBRARY $THREADS $FLANK $RUNS
# 
# # make directories
# mkdir -p data/run_0/raw out/
# 
# # cluster
# cd-hit-est -n 10 -c 0.95 -i ${RM_LIBRARY_PATH} -o data/run_0/${RM_LIBRARY} # cluster seq
# 
# # split
# PIECES="$(grep -c '>' data/run_0/${RM_LIBRARY})"
# Rscript scripts/splitter.R -t nt -f data/run_0/${RM_LIBRARY} -o data/run_0/raw out/ -p $PIECES # split query
# ls data/run_0/raw/${RM_LIBRARY}* | sed 's/.*\///' > data/run_0/to_run.txt # list queries
# 
# if [ ! -f "${GENOME}".nsq ]; then
# makeblastdb -in ${GENOME} -dbtype nucl -out ${GENOME} # makeblastb if needed
# fi
# 
# # prestrain
# mkdir -p data/run_0/rps_out/ data/run_0/trf_out/
# parallel --bar --jobs $THREADS -a data/run_0/to_run.txt 'rpstblastn -query data/run_0/raw/{} -db /media/projectDrive_1/databases/cdd/Cdd -out data/run_0/rps_out/{}.out -outfmt "6 qseqid qstart qend qlen sseqid sstart send slen pident length mismatch gapopen evalue bitscore qcovs stitle" -evalue 0.01 -num_threads 1'
# parallel --bar --jobs 1 -a data/run_0/to_run.txt cat data/run_0/rps_out/{}.out >> data/run_0/${RM_LIBRARY}.rps.out
# # parallel --bar --jobs $threads -a data/run_0/to_run.txt 'trf data/run_0/trf_out/{} 2 7 7 80 10 50 500 -d -h'
# # parallel --bar --jobs $threads -a data/run_0/to_run.txt 'python scripts/trf2gff.py -d {}.2.7.7.80.10.50.500.dat -o data/run_0/trf_out/{}.trf.gff'
# # cat data/run_0/raw/{}*.trf.gff > data/run_0/${RM_LIBRARY}.trf.gff
# Rscript strainer.R --in_seq data/run_0/${RM_LIBRARY} --out data/run_0/ --rps_table data/run_0/${RM_LIBRARY}.rps.out
date
# runs
RUN_NO=1
while  [ $RUN_NO -le $RUNS ]
do
  
  # split
  mkdir -p data/run_${RUN_NO}/raw
  cp data/run_$(expr $RUN_NO - 1)/clean_${RM_LIBRARY} data/run_${RUN_NO}/${RM_LIBRARY}
  PIECES=$(grep -c '>' data/run_${RUN_NO}/${RM_LIBRARY})
  echo "Splitting run "${RUN_NO}
  Rscript scripts/splitter.R -t nt -f data/run_${RUN_NO}/${RM_LIBRARY} -o data/run_${RUN_NO}/raw -p $PIECES
  ls data/run_${RUN_NO}/raw/${RM_LIBRARY}_* | sed 's/.*\///' > data/run_${RUN_NO}/to_run.txt

  # extend/align
  ## initial blast
  mkdir -p data/run_${RUN_NO}/initial_blast
  echo "Initial blast "${RUN_NO}
  parallel --env --bar --jobs ${THREADS} -a data/run_${RUN_NO}/to_run.txt blastn -task dc-megablast -query data/run_${RUN_NO}/raw/{} -db $GENOME -evalue 1e-5 -outfmt \"6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore qcovs\" -out data/run_${RUN_NO}/initial_blast/{}.out -num_threads 1 # search genome
  cat data/run_${RUN_NO}/initial_blast/*.out > data/run_${RUN_NO}/${RM_LIBRARY}_initial_blast.out # compile data

  ## extention and check blast
  echo "Extension run "${RUN_NO}
  Rscript scripts/self_blast_setup.R -g ${GENOME} -l ${RM_LIBRARY} -n ${RUN_NO} # extend seqs
  ls data/run_${RUN_NO}/initial_seq/*fasta | sed 's/.*\///' > data/run_${RUN_NO}/self_queries.txt
  mkdir data/run_${RUN_NO}/self_search/
  echo "Extension check run "${RUN_NO}
  parallel --bar --jobs ${THREADS} -a data/run_${RUN_NO}/self_queries.txt blastn -task dc-megablast -query data/run_${RUN_NO}/initial_seq/{} -subject data/run_${RUN_NO}/initial_seq/{} -evalue 1e-5 -outfmt \"6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore\" -out data/run_${RUN_NO}/self_search/{}.out -num_threads 1 # self blast

  ## clip seqs and align
  echo "Clipping run "${RUN_NO}
  mkdir -p data/run_${RUN_NO}/to_align data/run_${RUN_NO}/mafft data/run_${RUN_NO}/TEtrim_complete
  Rscript scripts/mafft_setup.R -g ${GENOME} -n ${RUN_NO} -l ${RM_LIBRARY} # trim seqs pre-mafft

  ## align seqs
  echo "Primary alignment run "${RUN_NO}
  if [[ $THREADS -gt 4 ]]; then MAFFT_THREADS=$(($(($THREADS / 4)))); else MAFFT_THREADS=1; fi
  parallel --bar --jobs $MAFFT_THREADS -a data/run_${RUN_NO}/to_align.txt "mafft --thread 4 --quiet --localpair --adjustdirectionaccurately data/run_${RUN_NO}/to_align/{} > data/run_${RUN_NO}/mafft/{}"
  ls data/run_${RUN_NO}/mafft/ > data/run_${RUN_NO}/${RM_LIBRARY}_to_trim.txt

  # trim
  mkdir -p data/run_${RUN_NO}/TEtrim_con \
           data/run_${RUN_NO}/TEtrim_unaln \
           data/run_${RUN_NO}/TEtrim_blast \
           data/run_${RUN_NO}/TEtrim_mafft \
           data/run_${RUN_NO}/TEtrim_further \
           data/run_${RUN_NO}/TEtrim_bp

  echo "Trimming run "${RUN_NO}
  parallel --bar --jobs $MAFFT_THREADS -a data/run_${RUN_NO}/${RM_LIBRARY}_to_trim.txt python scripts/TEtrim.py -i data/run_${RUN_NO}/mafft/{} -t 4 -f 1500 -n ${RUN_NO}

  # strain further needed to improve # (add flipper to this)
  echo "Straining further needed run "${RUN_NO}
  mkdir -p data/run_${RUN_NO}/strain_further_rps
  ls data/run_${RUN_NO}/TEtrim_further/ > data/run_${RUN_NO}/further_strain.txt
  parallel --bar --jobs $THREADS -a data/run_${RUN_NO}/further_strain.txt rpstblastn -query data/run_${RUN_NO}/TEtrim_further/{} -db /media/projectDrive_1/databases/cdd/Cdd -out data/run_${RUN_NO}/strain_further_rps/{}.out -outfmt \"6 qseqid qstart qend qlen sseqid sstart send slen pident length mismatch gapopen evalue bitscore qcovs stitle\" -evalue 0.01 -num_threads 1
  parallel --bar --jobs 1 -a data/run_${RUN_NO}/further_strain.txt cat data/run_${RUN_NO}/strain_further_rps/{}.out >> data/run_${RUN_NO}/${RM_LIBRARY}.rps.out
  cat data/run_${RUN_NO}/TEtrim_further/* > data/run_${RUN_NO}/further_${RM_LIBRARY}
  Rscript strainer.R --in_seq data/run_${RUN_NO}/further_${RM_LIBRARY} --out data/run_${RUN_NO}/ --rps_table data/run_${RUN_NO}/${RM_LIBRARY}.rps.out

  # strain complete
  echo "Straining complete run "${RUN_NO}
  mkdir -p data/run_${RUN_NO}/strain_complete_rps
  ls data/run_${RUN_NO}/TEtrim_complete/ > data/run_${RUN_NO}/complete_strain.txt
  parallel --bar --jobs $THREADS -a data/run_${RUN_NO}/complete_strain.txt rpstblastn -query data/run_${RUN_NO}/TEtrim_complete/{} -db /media/projectDrive_1/databases/cdd/Cdd -out data/run_${RUN_NO}/strain_complete_rps/{}.out -outfmt \"6 qseqid qstart qend qlen sseqid sstart send slen pident length mismatch gapopen evalue bitscore qcovs stitle\" -evalue 0.01 -num_threads 1
  parallel --bar --jobs 1 -a data/run_${RUN_NO}/complete_strain.txt cat data/run_${RUN_NO}/strain_complete_rps/{}.out >> data/run_${RUN_NO}/complete_${RM_LIBRARY}.rps.out
  cat data/run_${RUN_NO}/TEtrim_complete/* > data/run_${RUN_NO}/complete_${RM_LIBRARY}
  Rscript strainer.R --in_seq data/run_${RUN_NO}/complete_${RM_LIBRARY} --out data/run_${RUN_NO}/ --rps_table data/run_${RUN_NO}/complete_${RM_LIBRARY}.rps.out

  # compile data
  echo "Compile run "${RUN_NO}
  cat data/run_${RUN_NO}/chimeric_further_${RM_LIBRARY} >> data/run_${RUN_NO}/chimeric_complete_${RM_LIBRARY}
  cat data/run_${RUN_NO}/questionable_further_${RM_LIBRARY} >> data/run_${RUN_NO}/questionable_complete_${RM_LIBRARY}
  cp data/run_${RUN_NO}/clean_further_${RM_LIBRARY} data/run_${RUN_NO}/clean_${RM_LIBRARY}

 echo "ready for " $(( RUN_NO++ ))

done
  
date
  
  
  # mkdir data/run_${RUN_NO}/strain_trf
  # parallel --bar --jobs $threads -a data/run_0/to_run.txt 'trf data/run_0/trf_out/{} 2 7 7 80 10 50 500 -d -h'
  # parallel --bar --jobs $threads -a data/run_0/to_run.txt 'python scripts/trf2gff.py -d {}.2.7.7.80.10.50.500.dat -o data/run_0/trf_out/{}.trf.gff'
  # cat data/run_0/raw/{}*.trf.gff > data/run_0/${RM_LIBRARY}.trf.gff
  
  