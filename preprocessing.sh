#!/bin/bash

## -- Conda Environment
## conda create -n gostripes -c conda-forge -c bioconda \
## fastqc multiqc parallel seqkit star samtools cutadapt csvtk

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
    -p|--paired)
      PAIRED=true
      shift
      ;;
  esac
done

if [[ -z ${WORKDIR+x} ]]; then WORKDIR="."; fi
if [[ -z ${CORES+x} ]]; then CORES=1; fi

## Setup.

find bin -name "*sh" | xargs -n1 source

cd $WORKDIR
export CORES

## Preliminary fastq quality control.

if [[ ! -d seqs_qc ]]; then mkdir seqs_qc; fi

if [[ $PAIRED = true ]]; then
  cat  \
    <(csvtk cut -f fastq_1 $SAMPLES | csvtk del-header) \
    <(csvtk cut -f fastq_2 $SAMPLES | csvtk del-header) | \
    xargs fastqc -o seqs_qc -t $CORES
else
  csvtk cut -f fastq_1 $SAMPLES | csvtk del-header | \
    xargs fastqc -o seqs_qc -t $CORES
fi

## Keep only R1 reeds with N8TATAG3, and remove UMI, spacer, and poly-G.
## Tolerate 1 mismatched base.

if [[ ! -d fastqs/trimmed ]]; then mkdir -p fastqs/trimmed; fi

if [[ $PAIRED = true ]]; then
  parallel --link \
    -a <(csvtk cut -f fastq_1 $SAMPLES | csvtk del-header) \
    -a <(csvtk cut -f fastq_2 $SAMPLES | csvtk del-header) \
    -a <(csvtk cut -f name $SAMPLES | csvtk del-header) \
    'cutadapt \
      -g ^NNNNNNNNTATAGGG \
      -j $CORES \
      -e 1 \
      --discard-untrimmed \
      -o fastqs/trimmed/{1/} \
      -p fastqs/trimmed/{2/} \
      {1} {2} \
      1> fastqs/trimmed/{3}_log.txt'
else
  parallel --link \
    -a <(csvtk cut -f fastq_1 $SAMPLES | csvtk del-header) \
    -a <(csvtk cut -f name $SAMPLES | csvtk del-header) \
    'cutadapt \
      -g ^NNNNNNNNTATAGGG \
      -j $CORES \
      -e 1 \
      --discard-untrimmed \
      -o fastqs/trimmed/{1/} \
      {1} \
      1> fastqs/trimmed/{2}_log.txt'
fi

## Processed fastq quality control.

if [[ ! -d fastqs/trimmed_qc ]]; then mkdir fastqs/trimmed_qc; fi

find fastqs/trimmed -name "*fastq" | xargs fastqc -o fastqs/trimmed_qc -t $CORES

## Create STAR genome index.

if [[ ! -d index ]]; then mkdir index; fi

STAR \
  --runThreadN $CORES \
  --runMode genomeGenerate \
  --genomeDir index \
  --genomeFastaFiles $ASSEMBLY \
  --sjdbGTFfile $GTF \
  --genomeSAindexNbases 10

## Align with STAR.

if [[ ! -d bams/aligned ]]; then mkdir -p bams/aligned; fi

if [[ $PAIRED = true ]]; then
  parallel --link \
    -a <(csvtk cut -f fastq_1 $SAMPLES | csvtk del-header) \
    -a <(csvtk cut -f fastq_2 $SAMPLES | csvtk del-header) \
    -a <(csvtk cut -f name $SAMPLES | csvtk del-header) \
    'STAR \
      --runThreadN $CORES \
      --genomeDir index \
      --readFilesIn fastqs/trimmed/{1/} fastqs/trimmed/{2/} \
      --outFileNamePrefix bams/aligned/{3}_'
else
  parallel --link \
    -a <(csvtk cut -f fastq_1 $SAMPLES | csvtk del-header) \
    -a <(csvtk cut -f name $SAMPLES | csvtk del-header) \
    'STAR \
      --runThreadN $CORES \
      --genomeDir index \
      --readFilesIn fastqs/trimmed/{1/} \
      --outFileNamePrefix bams/aligned/{2}_'
fi

## Remove PCR duplicates and other poor reads.

if [[ ! -d bams/cleaned ]]; then mkdir bams/cleaned; fi

if [[ $PAIRED = true ]]; then
  parallel \
    -a <(csvtk cut -f name $SAMPLES | csvtk del-header) \
    'samtools sort -n -@ $CORES bams/aligned/{}_Aligned.out.sam | \
      samtools fixmate -m - - | \
      samtools sort -@ $CORES - | \
      samtools markdup - - | \
      samtools view -F 3852 -f 3 -O BAM -@ $CORES \
      -o bams/cleaned/{}.bam'
else
  parallel \
    -a <(csvtk cut -f name $SAMPLES | csvtk del-header) \
    'samtools sort -@ $CORES bams/aligned/{}_Aligned.out.sam | \
      samtools view -F 2820 -O BAM -@ $CORES \
      -o bams/cleaned/{}.bam'
fi

## Index the BAMs.

parallel -j$CORES 'samtools index "{}"' ::: bams/cleaned/*bam

## Multiqc report.

if [[ ! -d multiqc ]]; then mkdir multiqc; fi

multiqc -o multiqc -d .
