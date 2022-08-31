#!/bin/bash

usage() { echo "Usage: [-l Repeat library] [-g Genome ] [-t Threads] [-f Flank ] [-r Runs] [-c Set if clustering ] [-h Print this help]" 1>&2; exit 1; }

# default values
FLANK=1500
THREADS=4
RUNS=0
CLUSTER=FALSE

# parsing
while getopts l:g:t:f:r:ch flag; do
  case "${flag}" in
      l) RM_LIBRARY_PATH=${OPTARG};;
      g) GENOME=${OPTARG};;
      t) THREADS=${OPTARG};;
      f) FLANK=${OPTARG};;
      r) RUNS=${OPTARG};;
      c) CLUSTER=TRUE ;;
      h | *)
        print_usage
        exit_script
  esac
done

RM_LIBRARY=$(echo $RM_LIBRARY_PATH | sed 's/.*\///')
DATA_DIR=$(echo "TEstrainer_"$(date +"%H%M%S_%d%b"))
echo $GENOME $RM_LIBRARY $THREADS $FLANK $RUNS $DATA_DIR

# make directories
mkdir -p ${DATA_DIR}/run_0/raw out/ ${DATA_DIR}/curated/
# cluster and split step
if [ "$CLUSTER" == TRUE ]; then
    cd-hit-est -n 10 -c 0.95 -i ${RM_LIBRARY_PATH} -o ${DATA_DIR}/run_0/${RM_LIBRARY} # cluster seq
    PIECES="$(grep -c '>' ${DATA_DIR}/run_0/${RM_LIBRARY})"
else
    cp ${RM_LIBRARY_PATH} ${DATA_DIR}/run_0/${RM_LIBRARY}
    PIECES="$(grep -c '>' ${DATA_DIR}/run_0/${RM_LIBRARY})"
fi

Rscript scripts/splitter.R -t nt -f ${DATA_DIR}/run_0/${RM_LIBRARY} -o ${DATA_DIR}/run_0/raw/ -p $PIECES # split query

# prestrain
mkdir -p ${DATA_DIR}/run_0/rps_out/ ${DATA_DIR}/run_0/trf_out/
echo "Prestrain"
parallel --bar --jobs $THREADS -a ${DATA_DIR}/run_0/raw/${RM_LIBRARY}_split.txt rpstblastn -query ${DATA_DIR}/run_0/raw/{} -db /media/projectDrive_1/databases/cdd/Cdd -out ${DATA_DIR}/run_0/rps_out/{}.out -outfmt \"6 qseqid qstart qend qlen sseqid sstart send slen pident length mismatch gapopen evalue bitscore qcovs stitle\" -evalue 0.01 -num_threads 1
parallel --bar --jobs 1 -a ${DATA_DIR}/run_0/raw/${RM_LIBRARY}_split.txt cat ${DATA_DIR}/run_0/rps_out/{}.out >> ${DATA_DIR}/run_0/${RM_LIBRARY}.rps.out
Rscript strainer.R --in_seq ${DATA_DIR}/run_0/${RM_LIBRARY} --out ${DATA_DIR}/run_0/ --rps_table ${DATA_DIR}/run_0/${RM_LIBRARY}.rps.out
date

# runs
if [[ $RUNS -gt 0 ]]; then
  
  if [ ! -f "${GENOME}".nsq ]; then
    makeblastdb -in ${GENOME} -dbtype nucl -out ${GENOME} # makeblastb if needed
  fi

  RUN_NO=1
  while  [ $RUN_NO -le $RUNS ]
  do
  
    # split
    mkdir -p ${DATA_DIR}/run_${RUN_NO}/raw
    cp ${DATA_DIR}/run_$(expr $RUN_NO - 1)/clean_${RM_LIBRARY} ${DATA_DIR}/run_${RUN_NO}/${RM_LIBRARY}
    PIECES=$(grep -c '>' ${DATA_DIR}/run_${RUN_NO}/${RM_LIBRARY})
    echo "Splitting run "${RUN_NO}
    Rscript scripts/splitter.R -t nt -f ${DATA_DIR}/run_${RUN_NO}/${RM_LIBRARY} -o ${DATA_DIR}/run_${RUN_NO}/raw -p $PIECES
    ls ${DATA_DIR}/run_${RUN_NO}/raw/${RM_LIBRARY}_* | sed 's/.*\///' > ${DATA_DIR}/run_${RUN_NO}/to_run.txt
  
    # extend/align
    ## initial blast
    mkdir -p ${DATA_DIR}/run_${RUN_NO}/initial_blast ${DATA_DIR}/run_${RUN_NO}/initial_blast/TEtrim_complete
    echo "Initial blast "${RUN_NO}
    parallel --bar --jobs ${THREADS} -a ${DATA_DIR}/run_${RUN_NO}/to_run.txt blastn -task dc-megablast -query ${DATA_DIR}/run_${RUN_NO}/raw/{} -db $GENOME -evalue 1e-5 -outfmt \"6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore qcovs\" -out ${DATA_DIR}/run_${RUN_NO}/initial_blast/{}.out -num_threads 1 # search genome
    cat ${DATA_DIR}/run_${RUN_NO}/initial_blast/*.out > ${DATA_DIR}/run_${RUN_NO}/${RM_LIBRARY}_initial_blast.out # compile data
  
    ## extention and check blast
    echo "Extension run "${RUN_NO}
    Rscript scripts/self_blast_setup.R -g ${GENOME} -l ${RM_LIBRARY} -n ${RUN_NO} -d ${DATA_DIR} # extend seqs
    ls ${DATA_DIR}/run_${RUN_NO}/initial_seq/*fasta | sed 's/.*\///' > ${DATA_DIR}/run_${RUN_NO}/self_queries.txt
    mkdir ${DATA_DIR}/run_${RUN_NO}/self_search/
    echo "Extension check run "${RUN_NO}
    parallel --bar --jobs ${THREADS} -a ${DATA_DIR}/run_${RUN_NO}/self_queries.txt blastn -task dc-megablast -query ${DATA_DIR}/run_${RUN_NO}/initial_seq/{} -subject ${DATA_DIR}/run_${RUN_NO}/initial_seq/{} -evalue 1e-5 -outfmt \"6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore\" -out ${DATA_DIR}/run_${RUN_NO}/self_search/{}.out -num_threads 1 # self blast
  
    ## clip seqs and align
    echo "Clipping run "${RUN_NO}
    mkdir -p ${DATA_DIR}/run_${RUN_NO}/to_align ${DATA_DIR}/run_${RUN_NO}/mafft ${DATA_DIR}/run_${RUN_NO}/TEtrim_complete
    Rscript scripts/mafft_setup.R -g ${GENOME} -n ${RUN_NO} -l ${RM_LIBRARY} -d ${DATA_DIR} # trim seqs pre-mafft
    # ls > ${DATA_DIR}/run_${RUN_NO}/to_align.txt
  
    ## align seqs
    echo "Primary alignment run "${RUN_NO}
    if [[ $THREADS -gt 4 ]]; then MAFFT_THREADS=$(($(($THREADS / 4)))); else MAFFT_THREADS=1; fi
    parallel --bar --jobs $MAFFT_THREADS -a ${DATA_DIR}/run_${RUN_NO}/to_align.txt "mafft --thread 4 --quiet --localpair --adjustdirectionaccurately ${DATA_DIR}/run_${RUN_NO}/to_align/{} > ${DATA_DIR}/run_${RUN_NO}/mafft/{}"
    ls ${DATA_DIR}/run_${RUN_NO}/mafft/ > ${DATA_DIR}/run_${RUN_NO}/${RM_LIBRARY}_to_trim.txt
  
    # trim
    mkdir -p ${DATA_DIR}/run_${RUN_NO}/TEtrim_con \
             ${DATA_DIR}/run_${RUN_NO}/TEtrim_unaln \
             ${DATA_DIR}/run_${RUN_NO}/TEtrim_blast \
             ${DATA_DIR}/run_${RUN_NO}/TEtrim_mafft \
             ${DATA_DIR}/run_${RUN_NO}/TEtrim_further \
             ${DATA_DIR}/run_${RUN_NO}/TEtrim_bp
  
    echo "Trimming run "${RUN_NO}
    parallel --bar --jobs $MAFFT_THREADS -a ${DATA_DIR}/run_${RUN_NO}/${RM_LIBRARY}_to_trim.txt python scripts/TEtrim.py -i ${DATA_DIR}/run_${RUN_NO}/mafft/{} -t 4 -f 1500 -n ${RUN_NO} -d ${DATA_DIR}
  
    # strain further needed to improve # (add flipper to this)
    echo "Straining further needed run "${RUN_NO}
    mkdir -p ${DATA_DIR}/run_${RUN_NO}/strain_further_rps
    ls ${DATA_DIR}/run_${RUN_NO}/TEtrim_further/ > ${DATA_DIR}/run_${RUN_NO}/further_strain.txt
    parallel --bar --jobs $THREADS -a ${DATA_DIR}/run_${RUN_NO}/further_strain.txt rpstblastn -query ${DATA_DIR}/run_${RUN_NO}/TEtrim_further/{} -db /media/projectDrive_1/databases/cdd/Cdd -out ${DATA_DIR}/run_${RUN_NO}/strain_further_rps/{}.out -outfmt \"6 qseqid qstart qend qlen sseqid sstart send slen pident length mismatch gapopen evalue bitscore qcovs stitle\" -evalue 0.01 -num_threads 1
    parallel --bar --jobs 1 -a ${DATA_DIR}/run_${RUN_NO}/further_strain.txt cat ${DATA_DIR}/run_${RUN_NO}/strain_further_rps/{}.out >> ${DATA_DIR}/run_${RUN_NO}/further_${RM_LIBRARY}.rps.out
    cat ${DATA_DIR}/run_${RUN_NO}/TEtrim_further/* > ${DATA_DIR}/run_${RUN_NO}/further_${RM_LIBRARY}
    Rscript strainer.R --in_seq ${DATA_DIR}/run_${RUN_NO}/further_${RM_LIBRARY} --out ${DATA_DIR}/run_${RUN_NO}/ --rps_table ${DATA_DIR}/run_${RUN_NO}/further_${RM_LIBRARY}.rps.out
  
    # strain complete
    echo "Straining complete run "${RUN_NO}
    mkdir -p ${DATA_DIR}/run_${RUN_NO}/strain_complete_rps
    ls ${DATA_DIR}/run_${RUN_NO}/TEtrim_complete/ > ${DATA_DIR}/run_${RUN_NO}/complete_strain.txt
    parallel --bar --jobs $THREADS -a ${DATA_DIR}/run_${RUN_NO}/complete_strain.txt rpstblastn -query ${DATA_DIR}/run_${RUN_NO}/TEtrim_complete/{} -db /media/projectDrive_1/databases/cdd/Cdd -out ${DATA_DIR}/run_${RUN_NO}/strain_complete_rps/{}.out -outfmt \"6 qseqid qstart qend qlen sseqid sstart send slen pident length mismatch gapopen evalue bitscore qcovs stitle\" -evalue 0.01 -num_threads 1
    parallel --bar --jobs 1 -a ${DATA_DIR}/run_${RUN_NO}/complete_strain.txt cat ${DATA_DIR}/run_${RUN_NO}/strain_complete_rps/{}.out >> ${DATA_DIR}/run_${RUN_NO}/complete_${RM_LIBRARY}.rps.out
    cat ${DATA_DIR}/run_${RUN_NO}/TEtrim_complete/* > ${DATA_DIR}/run_${RUN_NO}/complete_${RM_LIBRARY}
    Rscript strainer.R --in_seq ${DATA_DIR}/run_${RUN_NO}/complete_${RM_LIBRARY} --out ${DATA_DIR}/run_${RUN_NO}/ --rps_table ${DATA_DIR}/run_${RUN_NO}/complete_${RM_LIBRARY}.rps.out
  
    # compile data
    echo "Compile run "${RUN_NO}
    cat ${DATA_DIR}/run_${RUN_NO}/chimeric_further_${RM_LIBRARY} ${DATA_DIR}/run_${RUN_NO}/chimeric_complete_${RM_LIBRARY} > ${DATA_DIR}/run_${RUN_NO}/chimeric_${RM_LIBRARY}
    rm ${DATA_DIR}/run_${RUN_NO}/chimeric_further_${RM_LIBRARY} ${DATA_DIR}/run_${RUN_NO}/chimeric_complete_${RM_LIBRARY}
    cat ${DATA_DIR}/run_${RUN_NO}/questionable_further_${RM_LIBRARY} ${DATA_DIR}/run_${RUN_NO}/questionable_complete_${RM_LIBRARY} > ${DATA_DIR}/run_${RUN_NO}/questionable_${RM_LIBRARY}
    rm ${DATA_DIR}/run_${RUN_NO}/questionable_further_${RM_LIBRARY} ${DATA_DIR}/run_${RUN_NO}/questionable_complete_${RM_LIBRARY}
    cp ${DATA_DIR}/run_${RUN_NO}/clean_further_${RM_LIBRARY} ${DATA_DIR}/run_${RUN_NO}/clean_${RM_LIBRARY}
  
    
  
    # exit loop if no more extension needed
    if [ ! -s ${DATA_DIR}/run_${RUN_NO}/clean_further_${RM_LIBRARY} ]; then
      echo Finished
    else
      echo "ready for " $(( RUN_NO++ ))
    fi
  
  done
  
  # compile extended and filter sequences
  cat ${DATA_DIR}/run_*/questionable_${RM_LIBRARY} > ${DATA_DIR}/curated/questionable_${RM_LIBRARY}
  cat ${DATA_DIR}/run_*/chimeric_${RM_LIBRARY} > ${DATA_DIR}/curated/chimeric_${RM_LIBRARY}
  cat ${DATA_DIR}/run_*/clean_complete_${RM_LIBRARY} > ${DATA_DIR}/curated/cleaned_${RM_LIBRARY}

  # add any those could have been extended further to clean
  if [ -s ${DATA_DIR}/run_${RUNS}/clean_further_${RM_LIBRARY} ]
  then
    cat ${DATA_DIR}/run_${RUNS}/clean_further_${RM_LIBRARY} >> ${DATA_DIR}/curated/cleaned_${RM_LIBRARY}
  fi

else

  # if no curation performed use initial sweep
  cp ${DATA_DIR}/run_0/*_${RM_LIBRARY} ${DATA_DIR}/curated/
  mv ${DATA_DIR}/curated/clean_${RM_LIBRARY} ${DATA_DIR}/curated/cleaned_${RM_LIBRARY}

fi

# Identify simple repeats and satellites, trim ends of LINEs/SINEs 
mkdir -p ${DATA_DIR}/trf/split
TRF_PIECES="$(grep -c '>' ${DATA_DIR}/curated/cleaned_${RM_LIBRARY})"
Rscript scripts/splitter.R -f ${DATA_DIR}/curated/cleaned_${RM_LIBRARY} -p ${TRF_PIECES} -t DNA -o ${DATA_DIR}/trf/split -n

# Running TRF
parallel --bar --jobs 128 -a ${DATA_DIR}/trf/split/cleaned_${RM_LIBRARY}_split.txt trf ${DATA_DIR}/trf/split/{} 2 7 7 80 10 50 500 -d -h -ngs ">" ${DATA_DIR}/trf/split/{}.trf
while read a; do
  awk '{OFS="\t"}{if ($4 >= 2) print FILENAME,$1,$2,$3,$4,$14}' ${DATA_DIR}/trf/split/${a}.trf | sed 's/.*\///;s/.fasta.trf//' >> ${DATA_DIR}/trf/cleaned_${RM_LIBRARY}.trf;
done<${DATA_DIR}/trf/split/cleaned_${RM_LIBRARY}_split.txt

# Running SA-SSR
sa-ssr -e -l 20 -L 50000 -m 1 -M 5000 -t 128 ${DATA_DIR}/curated/cleaned_${RM_LIBRARY} ${DATA_DIR}/trf/cleaned_${RM_LIBRARY}.sassr

# Running mreps and parse
parallel --bar --jobs 128 -a ${DATA_DIR}/trf/split/cleaned_${RM_LIBRARY}_split.txt bash scripts/mreps_parser.sh -i ${DATA_DIR}/trf/split/{}
while read a; do
  awk '{OFS="\t"}{print FILENAME,$1,$2,$3,$4,$5,$6,$7}' ${DATA_DIR}/trf/split/${a}.mreps.txt | sed 's/.*\///;s/.fasta.mreps.txt//' >> ${DATA_DIR}/trf/cleaned_${RM_LIBRARY}.mreps;
done<${DATA_DIR}/trf/split/cleaned_${RM_LIBRARY}_split.txt

# Interpret mreps, TRF and SA-SSR
Rscript scripts/simple_repeat_filter_trim.R -i ${DATA_DIR}/curated/cleaned_${RM_LIBRARY} -d ${DATA_DIR}

# Classify improved consensi using RepeatModeler's RepeatClassifier
cd ${DATA_DIR}/trf/
RepeatClassifier -pa $THREADS -consensi cleaned_${RM_LIBRARY}
cd ../../

# Trim chimeric elements
rpstblastn -query ${DATA_DIR}/curated/chimeric_${RM_LIBRARY} -db /media/projectDrive_1/databases/cdd/Cdd -out ${DATA_DIR}/curated/chimeric_${RM_LIBRARY}.out -outfmt "6 qseqid qstart qend qlen sseqid sstart send slen pident length mismatch gapopen evalue bitscore qcovs stitle" -evalue 0.01 -num_threads $THREADS
