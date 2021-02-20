#!/bin/bash

## -- Conda Environment
## conda create -n gostripes -c conda-forge -c bioconda \
## fastqc multiqc parallel seqtk star samtools cutadapt csvtk

## -- GNU Parallel Note
## System PERL interefered with GNU parallel.
## module unload PERL

## -- I used S288C STRIPE-seq as example data.
## module load sra-toolkit
## fasterq-dump -S -O seqs SRR10759419

set -e

## Command line arguments.

while [[ $# -gt 0 ]]; do
  key=$1
  case $key in
    -s|--sample-sheet)
      SAMPLES="$2"
      shift 2
      ;;
    -a|--assembly)
      ASSEMBLY="$2"
      shift 2
      ;;
    -g|--annotation)
      GTF="$2"
      shift 2
      ;;
    -o|--out-dir)
      WORKDIR="$2"
      shift 2
      ;;
    -t|--threads)
      CORES="$2"
      shift 2
      ;;
    --paired)
      PAIRED=true
      shift
      ;;
    --genomeSAindexNbases)
      GNB="$2"
      shift 2
      ;;
    -m|--max-mem)
      MAXMEM="$2"
      shift 2
      ;;
  esac
done

if [[ -z ${WORKDIR+x} ]]; then WORKDIR=$(pwd); fi
if [[ -z ${CORES+x} ]]; then CORES=1; fi
if [[ -z ${GNB+x} ]]; then GNB=14; fi
if [[ -z ${MAXMEM+x} ]]; then MAXMEM="768M"; fi

## Setup.

export CORES
export MAXMEM

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
for func in ${SCRIPT_DIR}/bin/*sh; do source $func; done

cd $WORKDIR

## Preliminary fastq quality control.

if [[ ! -d seqs_qc ]]; then mkdir seqs_qc; fi

printf "%-35s%s%s\n" "[$(date)]..." "Running FastQC on input FASTQs"
if [[ $PAIRED = true ]]; then
  cat  \
    <(csvtk cut -f fastq_1 $SAMPLES | csvtk del-header) \
    <(csvtk cut -f fastq_2 $SAMPLES | csvtk del-header) | \
    xargs fastqc -o seqs_qc -t $CORES
else
  csvtk cut -f fastq_1 $SAMPLES | csvtk del-header | \
    xargs fastqc -o seqs_qc -t $CORES
fi
printf "%-35s%s%s\n" "[$(date)]..." "Finished running FastQC"

## Keep only R1 reeds with N8TATAG3, and remove UMI, spacer, and poly-G.
## Tolerate 1 mismatched base.

if [[ ! -d fastqs/trimmed ]]; then mkdir -p fastqs/trimmed; fi

printf "%-35s%s%s\n" "[$(date)]..." "Filtering and trimming ^N{8}TATAGGG from FASTQs with Cutadapt"
if [[ $PAIRED = true ]]; then
  export -f CUTADAPT_PAIRED
  parallel --link \
    -a <(csvtk cut -f fastq_1 $SAMPLES | csvtk del-header) \
    -a <(csvtk cut -f fastq_2 $SAMPLES | csvtk del-header) \
    -a <(csvtk cut -f name $SAMPLES | csvtk del-header) \
    "CUTADAPT_PAIRED"
else
  export -f CUTADAPT_SINGLE
  parallel --link \
    -a <(csvtk cut -f fastq_1 $SAMPLES | csvtk del-header) \
    -a <(csvtk cut -f name $SAMPLES | csvtk del-header) \
    "CUTADAPT_SINGLE"
fi
printf "%-35s%s%s\n" "[$(date)]..." "Finished filtering and trimming FASTQs"

## Processed fastq quality control.

if [[ ! -d fastqs/trimmed_qc ]]; then mkdir fastqs/trimmed_qc; fi

printf "%-35s%s%s\n" "[$(date)]..." "Running FastQC on filtered and trimmed FASTQs"
find fastqs/trimmed -name "*fastq" | xargs fastqc -o fastqs/trimmed_qc -t $CORES
printf "%-35s%s%s\n" "[$(date)]..." "Finished running FastQC on filtered and trimmed FASTQs"

## Create STAR genome index.

if [[ ! -d index ]]; then mkdir index; fi

printf "%-35s%s%s\n" "[$(date)]..." "Generating the STAR genome index"
STAR \
  --runThreadN $CORES \
  --runMode genomeGenerate \
  --genomeDir index \
  --genomeFastaFiles $ASSEMBLY \
  --sjdbGTFfile $GTF \
  --genomeSAindexNbases $GNB
printf "%-35s%s%s\n" "[$(date)]..." "Finished generating the STAR genome index"

## Align with STAR.

if [[ ! -d bams/aligned ]]; then mkdir -p bams/aligned; fi

printf "%-35s%s%s\n" "[$(date)]..." "Aligning reads to genome using STAR"
if [[ $PAIRED = true ]]; then
  export -f STAR_PAIRED
  parallel --link \
    -a <(csvtk cut -f fastq_1 $SAMPLES | csvtk del-header) \
    -a <(csvtk cut -f fastq_2 $SAMPLES | csvtk del-header) \
    -a <(csvtk cut -f name $SAMPLES | csvtk del-header) \
    "STAR_PAIRED"
else
  export -f STAR_SINGLE
  parallel --link \
    -a <(csvtk cut -f fastq_1 $SAMPLES | csvtk del-header) \
    -a <(csvtk cut -f name $SAMPLES | csvtk del-header) \
    "STAR_SINGLE"
fi
printf "%-35s%s%s\n" "[$(date)]..." "Finished aligning reads"

## Remove PCR duplicates and other poor reads.

if [[ ! -d bams/cleaned ]]; then mkdir bams/cleaned; fi

printf "%-35s%s%s\n" "[$(date)]..." "Removing poor alignments and for paired-end data, PCR duplicates"
if [[ $PAIRED = true ]]; then
  export -f SAMTOOLS_PAIRED
  parallel \
    -a <(csvtk cut -f name $SAMPLES | csvtk del-header) \
    "SAMTOOLS_PAIRED"
else
  export -f SAMTOOLS_SINGLE
  parallel \
    -a <(csvtk cut -f name $SAMPLES | csvtk del-header) \
    "SAMTOOLS_SINGLE"
fi
printf "%-35s%s%s\n" "[$(date)]..." "Finished removing poor alignments, and for paired-end data, PCR duplicates"

## Index the BAMs.

printf "%-35s%s%s\n" "[$(date)]..." "Indexing BAMs"
parallel -j$CORES samtools index {} ::: bams/cleaned/*bam
printf "%-35s%s%s\n" "[$(date)]..." "Finished indexing BAMs"

## Multiqc report.

if [[ ! -d multiqc ]]; then mkdir multiqc; fi

printf "%-35s%s%s\n" "[$(date)]..." "Generating MultiQC report"
multiqc -o multiqc -d .
printf "%-35s%s%s\n" "[$(date)]..." "Finished making MultiQC report"
