#!/bin/bash

## Download reads.

if [ ! -d seqs ]; then mkdir seqs; fi

curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR107/019/SRR10759419/SRR10759419_1.fastq.gz | \
  gunzip > seqs/SRR10759419_1.fastq

curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR107/019/SRR10759419/SRR10759419_2.fastq.gz | \
  gunzip > seqs/SRR10759419_2.fastq

## Create sample sheet.

cat \
  <(echo name,fastq_1,fastq_2) \
  <(echo S288C,seqs/SRR10759419_1.fastq,seqs/SRR10759419_2.fastq) \
  > sample_sheet.csv

## Download genome FASTA and GTF.

if [ ! -d genome ]; then mkdir genome; fi

curl ftp://ftp.ensembl.org/pub/release-102/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.102.gtf.gz | \
  gunzip > genome/Saccharomyces_cerevisiae.R64-1-1.102.gtf

curl ftp://ftp.ensembl.org/pub/release-102/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz | \
  gunzip > genome/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa

## Run GOSTRIPEs.

bash ../preprocessing.sh \
  -s sample_sheet.csv \
  -a genome/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa \
  -g genome/Saccharomyces_cerevisiae.R64-1-1.102.gtf \
  -t 4 \
  -o $(pwd) \
  --paired

## Remove test files.

rm -r genome index fastqs seqs multiqc bams seqs_qc sample_sheet.csv
