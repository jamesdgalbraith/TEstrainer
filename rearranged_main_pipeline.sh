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

# determine further variables
RM_LIBRARY=$(echo $RM_LIBRARY_PATH | sed 's/.*\///')
DATA_DIR=$(echo "TEstrainer_"$(date +"%H%M%S_%d%b"))
if [[ $THREADS -gt 4 ]]; then MAFFT_THREADS=$(($(($THREADS / 4)))); else MAFFT_THREADS=1; fi
echo $GENOME $RM_LIBRARY $THREADS $FLANK $RUNS $DATA_DIR $MAFFT_THREADS

# make directories
mkdir -p ${DATA_DIR}/run_0/ out/ ${DATA_DIR}/curated/

# cluster and split step
if [ "$CLUSTER" == TRUE ]; then
  cd-hit-est -n 10 -c 0.95 -i ${RM_LIBRARY_PATH} -o ${DATA_DIR}/run_0/further_${RM_LIBRARY} # cluster seq
  PIECES="$(grep -c '>' ${DATA_DIR}/run_0/${RM_LIBRARY})"
else
  cp ${RM_LIBRARY_PATH} ${DATA_DIR}/run_0/further_${RM_LIBRARY}
  PIECES="$(grep -c '>' ${DATA_DIR}/run_0/further_${RM_LIBRARY})"
fi

# runs
if [[ $RUNS -gt 0 ]]; then
  
  if [ ! -f "${GENOME}".nsq ]; then
    makeblastdb -in ${GENOME} -dbtype nucl -out ${GENOME} # makeblastb if needed
  fi

  RUN_NO=1
  while  [ $RUN_NO -le $RUNS ]
  do
  
    # split
    mkdir -p ${DATA_DIR}/run_${RUN_NO}/raw ${DATA_DIR}/run_${RUN_NO}/TEtrim_complete
    cp ${DATA_DIR}/run_$(expr $RUN_NO - 1)/further_${RM_LIBRARY} ${DATA_DIR}/run_${RUN_NO}/${RM_LIBRARY}
    PIECES=$(grep -c '>' ${DATA_DIR}/run_${RUN_NO}/${RM_LIBRARY})
    echo "Splitting run "${RUN_NO}
    Rscript scripts/splitter.R -t nt -f ${DATA_DIR}/run_${RUN_NO}/${RM_LIBRARY} -o ${DATA_DIR}/run_${RUN_NO}/raw -p $PIECES
    
    # extend/align
    ## initial blast
    mkdir -p ${DATA_DIR}/run_${RUN_NO}/initial_blast ${DATA_DIR}/run_${RUN_NO}/initial_blast/TEtrim_complete
    echo "Initial blast "${RUN_NO}
    parallel --bar --jobs ${THREADS} -a ${DATA_DIR}/run_${RUN_NO}/raw/${RM_LIBRARY}_split.txt blastn -task dc-megablast -query ${DATA_DIR}/run_${RUN_NO}/raw/{} -db $GENOME -evalue 1e-5 -outfmt \"6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore qcovs\" -out ${DATA_DIR}/run_${RUN_NO}/initial_blast/{}.out -num_threads 1 # search genome
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
    mkdir -p ${DATA_DIR}/run_${RUN_NO}/to_align ${DATA_DIR}/run_${RUN_NO}/mafft
    Rscript scripts/mafft_setup.R -g ${GENOME} -n ${RUN_NO} -l ${RM_LIBRARY} -d ${DATA_DIR} # trim seqs pre-mafft
    
    ## align seqs
    echo "Primary alignment run "${RUN_NO}
    parallel --bar --jobs $MAFFT_THREADS -a ${DATA_DIR}/run_${RUN_NO}/to_align.txt mafft --thread 4 --quiet --localpair --adjustdirectionaccurately ${DATA_DIR}/run_${RUN_NO}/to_align/{} ">" ${DATA_DIR}/run_${RUN_NO}/mafft/{}
    
    # trim
    mkdir -p ${DATA_DIR}/run_${RUN_NO}/TEtrim_con \
             ${DATA_DIR}/run_${RUN_NO}/TEtrim_unaln \
             ${DATA_DIR}/run_${RUN_NO}/TEtrim_blast \
             ${DATA_DIR}/run_${RUN_NO}/TEtrim_mafft \
             ${DATA_DIR}/run_${RUN_NO}/TEtrim_further \
             ${DATA_DIR}/run_${RUN_NO}/TEtrim_bp
  
    echo "Trimming run "${RUN_NO}
    parallel --bar --jobs $MAFFT_THREADS -a ${DATA_DIR}/run_${RUN_NO}/to_align.txt python scripts/TEtrim.py -i ${DATA_DIR}/run_${RUN_NO}/mafft/{} -t 4 -f 1500 -n ${RUN_NO} -d ${DATA_DIR}
    cat ${DATA_DIR}/run_${RUN_NO}/TEtrim_complete/*fasta > ${DATA_DIR}/run_${RUN_NO}/complete_${RM_LIBRARY}
    cat ${DATA_DIR}/run_${RUN_NO}/TEtrim_further/*fasta > ${DATA_DIR}/run_${RUN_NO}/further_${RM_LIBRARY}
  
    # exit loop if no more extension needed
    if [ -s ${DATA_DIR}/run_${RUN_NO}/further_${RM_LIBRARY} ]; then
      echo "ready for " $(( RUN_NO++ ))
    else
      echo Finished
    fi
  
  done

  # Compile all completed
  cat ${DATA_DIR}/run_*/complete_${RM_LIBRARY} > ${DATA_DIR}/${RM_LIBRARY}
  # if more could have been extended append to compilation
  if [ -s ${DATA_DIR}/run_${RUNS}/further_${RM_LIBRARY} ]; then
    cat ${DATA_DIR}/run_${RUNS}/further_${RM_LIBRARY} >> ${DATA_DIR}/${RM_LIBRARY}
  fi
  
else

  # for if no curation performed
  cp ${DATA_DIR}/run_0/further_${RM_LIBRARY} ${DATA_DIR}/${RM_LIBRARY}

fi

# Identify simple repeats and satellites, trim ends of LINEs/SINEs
echo "Splitting for simple/satellite packages"
mkdir -p ${DATA_DIR}/trf/split
PIECES="$(grep -c '>' ${DATA_DIR}/${RM_LIBRARY})"
Rscript scripts/splitter.R -f ${DATA_DIR}/${RM_LIBRARY} -p ${PIECES} -t DNA -o ${DATA_DIR}/trf/split -n
# Running TRF
echo "Running TRF"
parallel --bar --jobs 128 -a ${DATA_DIR}/trf/split/${RM_LIBRARY}_split.txt trf ${DATA_DIR}/trf/split/{} 2 7 7 80 10 50 500 -d -h -ngs ">" ${DATA_DIR}/trf/split/{}.trf
parallel --bar --jobs 128 -a ${DATA_DIR}/trf/split/${RM_LIBRARY}_split.txt python3 scripts/trf_parser.py --trf ${DATA_DIR}/trf/split/{}.trf --out ${DATA_DIR}/trf/split/{}.trf.tsv
parallel --bar --jobs 128 -a ${DATA_DIR}/trf/split/${RM_LIBRARY}_split.txt cat ${DATA_DIR}/trf/split/{}.trf.tsv ">>" ${DATA_DIR}/trf/${RM_LIBRARY}.trf
# Running SA-SSR
echo "Running SA-SSR"
sa-ssr -e -l 20 -L 50000 -m 1 -M 5000 -t 128 ${DATA_DIR}/${RM_LIBRARY} ${DATA_DIR}/trf/${RM_LIBRARY}.sassr
# Running mreps and parse
echo "Running mreps"
parallel --bar --jobs 128 -a ${DATA_DIR}/trf/split/${RM_LIBRARY}_split.txt bash scripts/mreps_parser.sh -i ${DATA_DIR}/trf/split/{}
cat ${DATA_DIR}/trf/split/*.mreps.txt | sed 's/.*\///;s/.fasta.mreps.txt//' > ${DATA_DIR}/trf/${RM_LIBRARY}.mreps;
# Interpret mreps, TRF and SA-SSR
echo "Trimming and sorting based on mreps, TRF, SA-SSR"
Rscript scripts/simple_repeat_filter_trim.R -i ${DATA_DIR}/${RM_LIBRARY} -d ${DATA_DIR}
cp ${DATA_DIR}/trf/trimmed_${RM_LIBRARY} ${DATA_DIR}/${RM_LIBRARY}


# Identify and trim chimeric elements, remove proteins
mkdir -p ${DATA_DIR}/chimeras/split/
PIECES="$(grep -c '>' ${DATA_DIR}/${RM_LIBRARY})"
Rscript scripts/splitter.R -t nt -f ${DATA_DIR}/${RM_LIBRARY} -o ${DATA_DIR}/chimeras/split/ -p $PIECES
parallel --bar --jobs $THREADS -a ${DATA_DIR}/chimeras/split/${RM_LIBRARY}_split.txt rpstblastn -query ${DATA_DIR}/chimeras/split/{} -db /media/projectDrive_1/databases/cdd/Cdd -out ${DATA_DIR}/chimeras/split/{}.out -outfmt \"6 qseqid qstart qend qlen sseqid sstart send slen pident length mismatch gapopen evalue bitscore qcovs stitle\" -evalue 0.01 -num_threads 1
parallel --bar --jobs $THREADS -a ${DATA_DIR}/chimeras/split/${RM_LIBRARY}_split.txt cat ${DATA_DIR}/chimeras/split/{}.out > ${DATA_DIR}/chimeras/${RM_LIBRARY}.rps.out
# Rscript scripts/strainer.R --in_seq ${DATA_DIR}/${RM_LIBRARY} --out ${DATA_DIR}/chimeras/${RM_LIBRARY} --rps_table ${DATA_DIR}/chimeras/${RM_LIBRARY}.rps.out


# # Classify improved consensi using RepeatModeler's RepeatClassifier
# echo "Reclassifying repeats"
# mkdir -p ${DATA_DIR}/classify/
# cp ${DATA_DIR}/${RM_LIBRARY} ${DATA_DIR}/classify/
# cd ${DATA_DIR}/classify/
# RepeatClassifier -pa {$THREADS} -consensi ${RM_LIBRARY}
# # Classify improved consensi using a method based on RepeatModeler's RepeatClassifier
# echo "Reclassifying repeats"
# mkdir -p ${DATA_DIR}/classify/split
# PIECES="$(grep -c '>' ${DATA_DIR}/${RM_LIBRARY})"
# Rscript scripts/splitter.R -t nt -f ${DATA_DIR}/${RM_LIBRARY} -o ${DATA_DIR}/classify/split/ -p $PIECES
# echo "BLASTN"
# parallel --bar --jobs $THREADS -a ${DATA_DIR}/classify/split/${RM_LIBRARY}_split.txt blastn -db data/RepeatMasker.lib -xdrop_ungap 500 -xdrop_gap_final 1000 -xdrop_gap 125 -min_raw_gapped_score 250 -dust no -gapopen 25 -gapextend 5 -word_size 7 -outfmt \"6 std qlen slen stitle\" -query ${DATA_DIR}/classify/split/{} -out ${DATA_DIR}/classify/split/{}.nt.out
# cat ${DATA_DIR}/classify/split/${RM_LIBRARY}*.nt.out > ${DATA_DIR}/classify/${RM_LIBRARY}.nt.out
# echo "BLASTX"
# parallel --bar --jobs $THREADS -a ${DATA_DIR}/classify/split/${RM_LIBRARY}_split.txt blastx -db data/RepeatPeps.lib -word_size 2 -outfmt \"6 std qlen slen stitle\" -query ${DATA_DIR}/classify/split/{} -out ${DATA_DIR}/classify/split/{}.aa.out
# cat ${DATA_DIR}/classify/split/${RM_LIBRARY}*.aa.out > ${DATA_DIR}/classify/${RM_LIBRARY}.aa.out
# echo "TRF"
# parallel --bar --jobs $THREADS -a ${DATA_DIR}/classify/split/${RM_LIBRARY}_split.txt trf ${DATA_DIR}/classify/split/{} 2 7 7 80 10 50 500 -d -h -ngs ">" ${DATA_DIR}/classify/split/{}.trf
# parallel --bar --jobs $THREADS -a ${DATA_DIR}/classify/split/${RM_LIBRARY}_split.txt python3 scripts/trf_parser.py --trf ${DATA_DIR}/classify/split/{}.trf --out ${DATA_DIR}/classify/split/{}.trf.tsv
# parallel --bar --jobs $THREADS -a ${DATA_DIR}/classify/split/${RM_LIBRARY}_split.txt cat ${DATA_DIR}/classify/split/{}.trf.tsv ">>" ${DATA_DIR}/classify/${RM_LIBRARY}.trf
