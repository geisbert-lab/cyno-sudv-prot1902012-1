#!/bin/bash -l

# qsub options
#$ -l h_rt=24:00:00
#$ -pe omp 10
#$ -j y
#$ -o log-$JOB_NAME.qlog

## setup -----------------------------------------------------------------------
# functions
mesg () { echo -e "[MSG] $@"; }
err () { echo -e "[ERR] $@"; exit 1; }
checkcmd () {
  if [ $? -eq 0 ]
  then
    mesg "$@ succeeded"
  else
    err "$@ failed"
  fi
}

# help message
# pre-set variables
LOFREQ="pipeline/lofreq/lofreq"
THRESHOLD="100"
HELP="usage: qsub -P PROJECT -N JOBNAME $(basename "$0") OPTIONS

options (REQUIRED):
  -n NHP/host genome bowtie2 index
  -v viral genome bowtie2 index
  -f viral genome FASTA
  -o output directory
  -s sample ID
  -x forward FASTQ file
  -y reverse FASTQ file (optional)
  -t [OPTIONAL] minimum aligned read depth (default: $THRESHOLD)
  -h show this message and exit
"
# parsing arguments
while getopts ":hn:v:f:o:s:x:y:t:" opt
do
  case ${opt} in
    n ) HOST="${OPTARG}"
      ;;
    v ) VIRAL="${OPTARG}"
      ;;
    f ) REFSEQ="${OPTARG}"
      ;;
    o ) ODIR="${OPTARG}"
      ;;
    s ) SAMPLE="${OPTARG}"
      ;;
    x ) R1="${OPTARG}"
      ;;
    y ) R2="${OPTARG}"
      ;;
    t ) THRESHOLD="${OPTARG}"
      ;;
    h ) echo "$HELP" && exit 0
      ;;
    \? ) err "Invalid option ${opt}\n${HELP}"
      ;;
  esac
done
shift $((OPTIND -1))

## job info --------------------------------------------------------------------
echo "=========================================================="
echo "Start date: $(date)"
echo "Running on node: $(hostname)"
echo "Current directory: $(pwd)"
echo "Job name: $JOB_NAME"
echo "Job ID: $JOB_ID"
echo "=========================================================="
echo ""

## check inputs ----------------------------------------------------------------
mesg "STEP 0: CHECKING INPUTS"

# load bowtie2
module load bowtie2/2.4.2
checkcmd "Loading bowtie2"

# load samtools
module load htslib/1.18
module load samtools/1.18
checkcmd "Loading samtools"

# load R
module load R/4.0.2
checkcmd "Loading R"

# double-check that lofreq exists
if [ -z "$($LOFREQ version 2> /dev/null)" ]
then
  err "LoFreq error: $LOFREQ"
fi

# bowtie2 indicies should have 6 files
if [ -z "$HOST" ]
then
  err "No host index provided"
elif [ "$(ls -1 ${HOST}.* 2> /dev/null | wc -l)" -eq 6 ]
then
  mesg "Host index: $HOST"
else
  err "Invalid host index: $HOST"
fi
if [ -z "$VIRAL" ]
then
  err "No host index provided"
elif [ "$(ls -1 ${VIRAL}.* 2> /dev/null | wc -l)" -eq 6 ]
then
  mesg "Host index: $VIRAL"
else
  err "Invalid host index: $VIRAL"
fi

# reference sequence
if [ -z "$REFSEQ" ]
then
  err "No reference FASTA provided"
elif [ -f "$REFSEQ" ]
then
  mesg "Reference FASTA: $REFSEQ"
else
  err "Invalid reference FASTA: $REFSEQ"
fi

# R1/R0 FASTQ file
if [ -z "$R1" ]
then
  err "No FASTQ file provided (-x)"
elif [ -f "$R1" ]
then
  mesg "First FASTQ file: $R1"
else
  err "Invalid first FASTQ file: $R1"
fi

# check for R2 FASTQ file
if [ -z "$R2" ]
then 
  mesg "Only 1 FASTQ file detected. Running as unpaired."
elif [ -f "$R2" ]
then
  mesg "Second FASTQ file: $R2"
else
  err "Invalid second FASTQ file: $R2"
fi

# output directory
if [ -z "$ODIR" ]
then
  err "No output directory provided"
elif [ -d "$ODIR" ]
then
  mesg "Valid output directory: $ODIR"
else
  mesg "Creating output directory: $ODIR"
  mkdir -p "$ODIR"
fi

# sample ID
if [ -z "$SAMPLE" ]
then
  err "No sample ID provided"
else
  mesg "Working with sample: $SAMPLE"
  ODIR="$ODIR/$SAMPLE"
  # make output directory
  mkdir -p "$ODIR"
fi

# done checking inputs!
mesg "Done checking inputs!"
echo ""

## alignment to host ---------------------------------------------------------
mesg "STEP 1: DISCARD HOST READS"

# build command based on whether has paired reads
if [ -z "$R2" ] # unpaired
then
  CMD="bowtie2 --threads 8 --quiet -x '$HOST' --very-sensitive --un-conc-gz '$ODIR/filtered.fq.gz' -U '$R1' > /dev/null"
else
  CMD="bowtie2 --threads 8 --quiet -x '$HOST' --very-sensitive --un-conc-gz '$ODIR/filtered-r%.fq.gz' -1 '$R1' -2 '$R2' > /dev/null"
fi

# run alignment filter
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Host alignment" 

# update R1 and R2
if [ -z "$R2" ] # unpaired
then
  R1="$ODIR/filtered.fq.gz"
else
  R1="$ODIR/filtered-r1.fq.gz"
  R2="$ODIR/filtered-r2.fq.gz"
fi

## alignment to virus --------------------------------------------------------
mesg "STEP 2: ALIGN TO VIRUS"

# build command based on whether has paired reads
if [ -z "$R2" ] # unpaired
then
  CMD="bowtie2 --threads 8 -x '$VIRAL' -U '$R1' > '$ODIR/alignment.sam'"
else # paired
  CMD="bowtie2 --threads 8 -x '$VIRAL' -1 '$R1' -2 '$R2' > '$ODIR/alignment.sam'"
fi

# run alignment
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Viral alignment"

# compress SAM to BAM
CMD="samtools view --threads 8 -b -h '$ODIR/alignment.sam' > '$ODIR/alignment-raw.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Compression"
rm "$ODIR/alignment.sam"
echo ""

## process BAM -----------------------------------------------------------------
mesg "STEP 3: PROCESS ALIGNMENT"

# sort BAM
CMD="samtools sort --threads 8 '$ODIR/alignment-raw.bam' > '$ODIR/alignment-sorted.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Sorting"
rm "$ODIR/alignment-raw.bam"

# score indels to get final BAM
CMD="$LOFREQ indelqual --dindel --ref '$REFSEQ' '$ODIR/alignment-sorted.bam' > '$ODIR/alignment.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Indelqual"
rm "$ODIR/alignment-sorted.bam"

# index final BAM
CMD="samtools index '$ODIR/alignment.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Indexing"
echo ""

## calculate coverage ----------------------------------------------------------
mesg "STEP 4: CALCULATE COVERAGE"

# coverage with samtools depth
CMD="samtools depth --threads 8 -a -H '$ODIR/alignment.bam' > '$ODIR/coverage.tsv'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Coverage"
echo ""

## assemble consensus ----------------------------------------------------------
mesg "STEP 5: ASSEMBLE CONSENSUS"

# consensus with samtools
CMD="samtools consensus --threads 8 -a --use-qual --min-depth $THRESHOLD --call-fract 0.5 --mode simple  --output '$ODIR/consensus-tmp.fa' '$ODIR/alignment.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Consensus"

# update consensus header
mesg "Updating consensus header"
echo ">$SAMPLE" > "$ODIR/consensus.fa"
cat "$ODIR/consensus-tmp.fa" | grep "^[^>]" >> "$ODIR/consensus.fa"
rm "$ODIR/consensus-tmp.fa"
echo ""

## quantify SNVs ---------------------------------------------------------------
mesg "STEP 6: QUANTIFY SNVs"

# run lofreq
# keeping mapping quality parameters same between samtools and LoFreq
CMD="$LOFREQ call-parallel --pp-threads 8 --call-indels --min-cov $THRESHOLD --ref '$REFSEQ' '$ODIR/alignment.bam' > '$ODIR/snvs.vcf'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "LoFreq"
echo ""

## package version -------------------------------------------------------------
mesg "Pipeline complete! Printing package versions..."
module list
echo "LoFreq"
$LOFREQ version
echo ""
