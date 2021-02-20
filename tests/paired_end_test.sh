#!/bin/bash

## Download reads.

if [[ ! -d seqs ]]; then mkdir seqs; fi

curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR107/019/SRR10759419/SRR10759419_1.fastq.gz | \
  gunzip > seqs/SRR10759419_1.fastq

if [[ ! -f seqs/SRR10759419_1.fastq ]]; then
  >&2 echo "ERROR: seqs/SRR10759419_1.fastq doesn't exist."
  exit 1
fi

curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR107/019/SRR10759419/SRR10759419_2.fastq.gz | \
  gunzip > seqs/SRR10759419_2.fastq

if [[ ! -f seqs/SRR10759419_2.fastq ]]; then
  >&2 echo "ERROR: seqs/SRR10759419_2.fastq doesn't exist."
  exit 1
fi

## Create sample sheet.

cat \
  <(echo name,fastq_1,fastq_2) \
  <(echo S288C,seqs/SRR10759419_1.fastq,seqs/SRR10759419_2.fastq) \
  > sample_sheet.csv

if [[ ! -f sample_sheet.csv ]]; then
  >&2 echo "ERROR: sample_sheet.csv doesn't exist."
  exit 1
fi

## Download genome FASTA and GTF.

if [[ ! -d genome ]]; then mkdir genome; fi

curl ftp://ftp.ensembl.org/pub/release-102/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.102.gtf.gz | \
  gunzip > genome/Saccharomyces_cerevisiae.R64-1-1.102.gtf

if [[ ! -f genome/Saccharomyces_cerevisiae.R64-1-1.102.gtf ]]; then
  >&2 echo "ERROR: genome/Saccharomyces_cerevisiae.R64-1-1.102.gtf doesn't exist."
  exit 1
fi

curl ftp://ftp.ensembl.org/pub/release-102/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz | \
  gunzip > genome/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa

if [[ ! -f genome/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa ]]; then
  >&2 echo "ERROR: genome/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa doesn't exist."
  exit 1
fi

## Run GOSTRIPEs.

bash ../gostripes.sh \
  -s sample_sheet.csv \
  -a genome/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa \
  -g genome/Saccharomyces_cerevisiae.R64-1-1.102.gtf \
  -t 8 \
  -o $(pwd) \
  --paired

## Check outputs exists.

if [[ ! -f seqs_qc/SRR10759419_1_fastqc.html || ! -f seqs_qc/SRR10759419_2_fastqc.html ]]; then
  >&2 echo "ERROR: FastQC files missing."
fi

if [[ ! -f fastqs/trimmed/SRR10759419_1.fastq || ! -f fastqs/trimmed/SRR10759419_2.fastq || ]]; then
  >&2 echo "ERROR: Trimmed FASTQ files and/or trimming log missing."
fi

if [[ ! -f fastqs/trimmed_qc/SRR10759419_1_fastqc.html || ! -f fastqs/trimmed_qc/SRR10759419_2_fastqc.html ]]; then
  >&2 echo "ERROR: FastQC of trimmed FASTQ files missing."
fi

if [[ ! -f index/SA ]]; then
  >&2 echo "ERROR: STAR genome index is missing."
fi

if [[ ! -f bams/aligned/S288C_Aligned.out.sam ]]; then
  >&2 echo "ERROR: Aligned SAM files missing."
fi

if [[ ! -f bams/cleaned/S288C.bam || ! -f bams/cleaned/S288C.bam.bai ]]; then
  >&2 echo "ERROR: Sorted and cleaned BAM file and/or BAM file index missing."
fi

if [[ ! -f multiqc/multiqc_report.html ]]; then
  >&2 echo "ERROR: MultiQC report missing."
fi

## Remove test files.

#rm -r genome index fastqs seqs multiqc bams seqs_qc sample_sheet.csv
