#!/bin/bash

usage() { echo "Usage: [-l Repeat library] [-g Genome ] [-t Threads] [-f Flank ] [-r Runs] [-c Set if clustering ] [-D Set if fixing Dfam data ] [-h Print this help]" 1>&2; exit 1; }

# default values
FLANK=1500
THREADS=4
RUNS=0
CLUSTER=FALSE
DFAM=FALSE
CLASSIFY=FALSE
# for potential folder name
TIME=$(date +"%s")
TIME=${TIME: -4}

# parsing
while getopts l:g:t:f:r:d:cCDh flag; do
  case "${flag}" in
    l) RM_LIBRARY_PATH=${OPTARG};;
    g) GENOME=${OPTARG};;
    t) THREADS=${OPTARG};;
    f) FLANK=${OPTARG};;
    r) RUNS=${OPTARG};;
    d) DATA_DIR=${OPTARG};;
    c) CLUSTER=TRUE ;;
    C) CLASSIFY=TRUE;;
    D) DFAM=TRUE ;;
    h | *)
      print_usage
      exit_script
  esac
done

# determine further variables and check files exist
if [ ! -f ${RM_LIBRARY_PATH} ]; then echo "Library not found"; usage; fi
if [ -z ${RM_LIBRARY_PATH} ]; then echo "Library must be supplied"; usage; else RM_LIBRARY=$(echo $RM_LIBRARY_PATH | sed 's/.*\///'); fi
if [[ $RUNS -gt 0 ]]; then 
  if [ -z ${GENOME} ]; then echo "If refining genome must be supplied"; usage; fi
  if [ ! -f ${GENOME} ]; then echo "Refining genome not found"; usage; fi
fi
if [ -z "$DATA_DIR" ]; then DATA_DIR=$(echo "TS_"${RM_LIBRARY}"_"${TIME}); fi

# create data dir if missing
if [ ! -d ${DATA_DIR} ] 
then
    mkdir ${DATA_DIR}
fi


if [[ $THREADS -gt 4 ]]; then MAFFT_THREADS=$(($(($THREADS / 4)))); else MAFFT_THREADS=1; fi
echo $RM_LIBRARY_PATH $GENOME $RM_LIBRARY $THREADS $FLANK $RUNS $DATA_DIR $MAFFT_THREADS "Cluster is" ${CLUSTER} "Dfam is" ${DFAM}

# make directories
mkdir -p ${DATA_DIR}/run_0/ out/ ${DATA_DIR}/curated/

# initial copy
if [ "$DFAM" == TRUE ]; then Rscript scripts/Dfam_extractor.R -l ${RM_LIBRARY_PATH} -d ${DATA_DIR} ; else cp ${RM_LIBRARY_PATH} ${DATA_DIR}/${RM_LIBRARY};fi

# cluster and split step
if [ "$CLUSTER" == TRUE ]; then
  cd-hit-est -n 10 -c 0.95 -i ${RM_LIBRARY_PATH} -o ${DATA_DIR}/run_0/further_${RM_LIBRARY} # cluster seq
  PIECES="$(grep -c '>' ${DATA_DIR}/run_0/further_${RM_LIBRARY})"
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
      break
    fi

  done

  # Compile all completed
  cat ${DATA_DIR}/run_*/complete_${RM_LIBRARY} > ${DATA_DIR}/${RM_LIBRARY}
  # if more could have been extended append to compilation
  if [ -s ${DATA_DIR}/run_${RUNS}/further_${RM_LIBRARY} ]; then
    cat ${DATA_DIR}/run_${RUNS}/further_${RM_LIBRARY} >> ${DATA_DIR}/${RM_LIBRARY}
  fi
  sed -i 's/ .*//' ${DATA_DIR}/${RM_LIBRARY}

else

echo "No curation"

fi

# Identify simple repeats and satellites, trim ends of LINEs/SINEs
echo "Splitting for simple/satellite packages"
mkdir -p ${DATA_DIR}/trf/split
PIECES="$(grep -c '>' ${DATA_DIR}/${RM_LIBRARY})"
Rscript scripts/splitter.R -f ${DATA_DIR}/${RM_LIBRARY} -p ${PIECES} -t DNA -o ${DATA_DIR}/trf/split -n
# Running TRF
echo "Running TRF"
parallel --bar --jobs ${THREADS} -a ${DATA_DIR}/trf/split/${RM_LIBRARY}_split.txt trf ${DATA_DIR}/trf/split/{} 2 7 7 80 10 50 500 -d -h -ngs ">" ${DATA_DIR}/trf/split/{}.trf
parallel --bar --jobs ${THREADS} -a ${DATA_DIR}/trf/split/${RM_LIBRARY}_split.txt python3 scripts/trf_parser.py --trf ${DATA_DIR}/trf/split/{}.trf --out ${DATA_DIR}/trf/split/{}.trf.tsv
find ./${DATA_DIR}/trf/split/ -type f -name "*trf.tsv" -exec cat {} + | cat > ${DATA_DIR}/trf/${RM_LIBRARY}.trf
# Running SA-SSR
echo "Running SA-SSR"
sa-ssr -e -l 20 -L 50000 -m 1 -M 5000 -t ${THREADS} ${DATA_DIR}/${RM_LIBRARY} ${DATA_DIR}/trf/${RM_LIBRARY}.sassr
# Running mreps and parse
echo "Running mreps"
parallel --bar --jobs ${THREADS} -a ${DATA_DIR}/trf/split/${RM_LIBRARY}_split.txt bash scripts/mreps_parser.sh -i ${DATA_DIR}/trf/split/{}
find ./${DATA_DIR}/trf/split/ -type f -name "*mreps" -exec cat {} + | cat > ${DATA_DIR}/trf/${RM_LIBRARY}.mreps
# Interpret mreps, TRF and SA-SSR
echo "Trimming and sorting based on mreps, TRF, SA-SSR"
Rscript scripts/simple_repeat_filter_trim.R -i ${DATA_DIR}/${RM_LIBRARY} -d ${DATA_DIR}
cp ${DATA_DIR}/trf/trimmed_${RM_LIBRARY} ${DATA_DIR}/${RM_LIBRARY}

# Identify and trim chimeric elements, remove proteins
mkdir -p ${DATA_DIR}/chimeras/split/
PIECES="$(grep -c '>' ${DATA_DIR}/${RM_LIBRARY})"
cp ${DATA_DIR}/${RM_LIBRARY} ${DATA_DIR}/chimeras/prestrain_${RM_LIBRARY}
Rscript scripts/splitter.R -t nt -f ${DATA_DIR}/${RM_LIBRARY} -o ${DATA_DIR}/chimeras/split/ -p $PIECES
parallel --bar --jobs $THREADS -a ${DATA_DIR}/chimeras/split/${RM_LIBRARY}_split.txt rpstblastn -query ${DATA_DIR}/chimeras/split/{} -db /media/projectDrive_1/databases/cdd/Cdd -out ${DATA_DIR}/chimeras/split/{}.out -outfmt \"6 qseqid qstart qend qlen sseqid sstart send slen pident length mismatch gapopen evalue bitscore qcovs stitle\" -evalue 0.01 -num_threads 1
find ./${DATA_DIR}/chimeras/split/ -type f -name "*.out" -exec cat {} + | cat > ${DATA_DIR}/chimeras/${RM_LIBRARY}.rps.out
Rscript scripts/strainer.R --in_seq ${DATA_DIR}/${RM_LIBRARY} --directory ${DATA_DIR}
cat ${DATA_DIR}/chimeras/clean_${RM_LIBRARY} ${DATA_DIR}/chimeras/chimeric_${RM_LIBRARY} > ${DATA_DIR}/${RM_LIBRARY}

# Delete temp files
rm -r ${DATA_DIR}/*/split/

if [ "$CLASSIFY" == TRUE ]; then
  # Classify improved consensi using RepeatModeler's RepeatClassifier
  echo "Reclassifying repeats"
  mkdir -p ${DATA_DIR}/classify/
  cp ${DATA_DIR}/${RM_LIBRARY} ${DATA_DIR}/classify/
  cd ${DATA_DIR}/classify/
  RepeatClassifier -debug -pa ${THREADS} -consensi ${RM_LIBRARY}
  cd -
  cp ${DATA_DIR}/classify/${RM_LIBRARY}.classified ${DATA_DIR}/${RM_LIBRARY}
fi