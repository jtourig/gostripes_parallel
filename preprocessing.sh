#!/bin/bash

## -- Conda Environment
## conda install -n gostripes -c conda-forge -c bioconda \
## fastqc multiqc parallel seqkit star samtools

## -- GNU Parallel Note
## System PERL interefered with GNU parallel.
## module unload PERL

## -- I used S288C STRIPE-seq as example data.
## module load sra-toolkit
## fasterq-dump -S -O seqs SRR10759419

## Variables.

WORKDIR='.'
GTF='ftp://ftp.ensembl.org/pub/release-102/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.102.gtf.gz'
ASSEMBLY='ftp://ftp.ensembl.org/pub/release-102/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz'
SEQDIR=seqs

## Preliminary fastq quality control.
cd $WORKDIR

if [ ! -d seqs_qc ]; then
  mkdir seqs_qc
fi

find seqs -name "*fastq" | xargs fastqc -o seqs_qc -t 4

## Keep only R1 reeds with N8TATAG3, and remove UMI, spacer, and poly-G.

if [ ! -d processed ]; then
  mkdir processed
fi

parallel -j1 -k -I ,, 'seqkit grep -j4 -srP -p"^[ATGCN]{8}TATAG{3}" ",," | seqkit subseq -r 16:-1 > processed/"{/}"' ::: ${SEQDIR}/*_1*

## Properly pair R1 and R2 reads.

parallel -j1 -k --link 'seqkit grep -j4 -f <(seqkit seq -ni "{1}") "{2}" > processed/"{2/}"' ::: processed/*_1* ::: seqs/*_2*

## Processed fastq quality control.

if [ ! -d processed_qc ]; then
  mkdir processed_qc
fi

find processed "*fastq" | xargs fastqc -o processed_qc -t 4

## Download genome.

if [ ! -d genome ]; then
  mkdir genome
fi

parallel -j2 'curl {} | gunzip > genome/{/.}' ::: $GTF $ASSEMBLY

## Create STAR genome index.

if [ ! -d genome/index ]; then
  mkdir genome/index
fi

STAR \
  --runThreadN 4 \
  --runMode genomeGenerate \
  --genomeDir genome/index \
  --genomeFastaFiles genome/$(basename $ASSEMBLY .gz) \
  --sjdbGTFfile genome/$(basename $GTF .gz) \
  --genomeSAindexNbases 10

## Align with STAR.

if [ ! -d aligned ]; then
  mkdir aligned
fi

parallel -k --link -j1 \
  'STAR \
    --runThreadN 4 \
    --genomeDir genome/index \
    --readFilesIn "{1}" "{2}" \
    --outFileNamePrefix aligned/"$(echo {1/.} | cut -d_ -f1)"_' \
  ::: processed/*_1* \
  ::: processed/*_2*

## Remove PCR duplicates and other poor reads.

if [ ! -d cleaned ]; then
  mkdir cleaned
fi

parallel -k -j1 \
  'samtools sort -n "{}" | \
    samtools fixmate -m - - | \
    samtools sort - | \
    samtools markdup - - | \
    samtools view -F 3852 -f 3 -O BAM \
    -o cleaned/$(echo "{/}" | cut -f1 -d_).bam' \
  ::: aligned/*Aligned*

## Multiqc report.
if [ ! -d multiqc ]; then
  mkdir multiqc
fi

multiqc -o multiqc -d .
