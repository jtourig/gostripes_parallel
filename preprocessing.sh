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

## Variables.

CORES=4

WORKDIR='.'
SAMPLES='sample_sheet.csv'

GTF='ftp://ftp.ensembl.org/pub/release-102/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.102.gtf.gz'
ASSEMBLY='ftp://ftp.ensembl.org/pub/release-102/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz'

## Setup.

cd $WORKDIR
export CORES
export SAMPLES

## Preliminary fastq quality control.

if [ ! -d seqs_qc ]; then mkdir seq_qc; fi

find seqs -name "*fastq" | xargs fastqc -o seqs_qc -t $CORES

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

## Download genome.

if [ ! -d genome ]; then mkdir genome; fi

parallel -j$CORES 'curl {} | gunzip > genome/{/.}' ::: $GTF $ASSEMBLY

## Create STAR genome index.

if [ ! -d genome/index ]; then mkdir genome/index; fi

STAR \
  --runThreadN $CORES \
  --runMode genomeGenerate \
  --genomeDir genome/index \
  --genomeFastaFiles genome/$(basename $ASSEMBLY .gz) \
  --sjdbGTFfile genome/$(basename $GTF .gz) \
  --genomeSAindexNbases 10

## Align with STAR.

if [ ! -d aligned ]; then mkdir aligned; fi

parallel --link \
  -a <(csvtk cut -f fastq_1 $SAMPLES | csvtk del-header) \
  -a <(csvtk cut -f fastq_2 $SAMPLES | csvtk del-header) \
  -a <(csvtk cut -f name $SAMPLES | csvtk del-header) \
  'STAR \
    --runThreadN $CORES \
    --genomeDir genome/index \
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
