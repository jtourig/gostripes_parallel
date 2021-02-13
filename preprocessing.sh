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

## Command line arguments.

while [[ $# -gt 0 ]]; do
  key=$1
  case $key in
    -s|--sample-sheet)
      SAMPLES="$2"
      shift 2
      ;;
    -as|--assembly)
      ASSEMBLY="$2"
      shift 2
      ;;
    -an|--annotation)
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
  esac
done

## Setup.

cd $WORKDIR
export CORES

## Preliminary fastq quality control.

if [ ! -d seqs_qc ]; then mkdir seq_qc; fi

cat  \
  <(csvtk cut -f fastq_1 $SAMPLES | csvtk del-header) \
  <(csvtk cut -f fastq_2 $SAMPLES | csvtk del-header) | \
  xargs fastqc -o seqs_qc -t $CORES

## Keep only R1 reeds with N8TATAG3, and remove UMI, spacer, and poly-G.
## Tolerate 1 mismatched base.

if [ ! -d processed ]; then mkdir processed; fi

parallel --link \
  -a <(csvtk cut -f fastq_1 $SAMPLES | csvtk del-header) \
  -a <(csvtk cut -f fastq_2 $SAMPLES | csvtk del-header) \
  -a <(csvtk cut -f name $SAMPLES | csvtk del-header) \
  'cutadapt \
    -g ^NNNNNNNNTATAGGG \
    -j $CORES \
    -e 1 \
    --discard-untrimmed \
    -o processed/{1/} \
    -p processed/{2/} \
    {1} {2} \
    1> processed/{3}_log.txt'

## Processed fastq quality control.

if [ ! -d processed_qc ]; then mkdir processed_qc; fi

find processed -name "*fastq" | xargs fastqc -o processed_qc -t $CORES

## Create STAR genome index.

if [ ! -d index ]; then mkdir index; fi

STAR \
  --runThreadN $CORES \
  --runMode genomeGenerate \
  --genomeDir index \
  --genomeFastaFiles genome/$(basename $ASSEMBLY) \
  --sjdbGTFfile genome/$(basename $GTF) \
  --genomeSAindexNbases 10

## Align with STAR.

if [ ! -d aligned ]; then mkdir aligned; fi

parallel --link \
  -a <(csvtk cut -f fastq_1 $SAMPLES | csvtk del-header) \
  -a <(csvtk cut -f fastq_2 $SAMPLES | csvtk del-header) \
  -a <(csvtk cut -f name $SAMPLES | csvtk del-header) \
  'STAR \
    --runThreadN $CORES \
    --genomeDir index \
    --readFilesIn {1} {2} \
    --outFileNamePrefix aligned/{3}_'

## Remove PCR duplicates and other poor reads.

if [ ! -d cleaned ]; then mkdir cleaned; fi

parallel \
  -a <(csvtk cut -f name $SAMPLES | csvtk del-header) \
  'samtools sort -n -@ $CORES aligned/{}_Aligned.out.sam | \
    samtools fixmate -m - - | \
    samtools sort -@ $CORES - | \
    samtools markdup - - | \
    samtools view -F 3852 -f 3 -O BAM -@ $CORES \
    -o cleaned/{}.bam'

## Index the BAMs.

parallel -j$CORES 'samtools index "{}"' ::: cleaned/*bam

## Multiqc report.

if [ ! -d multiqc ]; then mkdir multiqc; fi

multiqc -o multiqc -d .
