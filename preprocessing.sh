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

CORES=4

WORKDIR='.'
SEQDIR=''

R1_IDENTIFIER='_1'
R2_IDENTIFIER='_2'

GTF='ftp://ftp.ensembl.org/pub/release-102/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.102.gtf.gz'
ASSEMBLY='ftp://ftp.ensembl.org/pub/release-102/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz'

## Setup.

cd $WORKDIR
export CORES

## Preliminary fastq quality control.

if [ ! -d seqs_qc ]; then mkdir seq_qc; fi

find seqs -name "*fastq" | xargs fastqc -o seqs_qc -t $CORES

## Keep only R1 reeds with N8TATAG3, and remove UMI, spacer, and poly-G.

if [ ! -d processed ]; then mkdir processed; fi

parallel -I,, \
  'seqkit grep -j$CORES -srP -p"^[ATGCN]{8}TATAG{3}" ",," | seqkit subseq -r 16:-1 > processed/"{/}"' \
  ::: ${SEQDIR}/*${R1_IDENTIFIER}*

## Properly pair R1 and R2 reads.

parallel --link \
  'seqkit grep -j$CORES -f <(seqkit seq -ni "{1}") "{2}" > processed/"{2/}"' \
  ::: processed/*${R1_IDENTIFIER}* \
  ::: seqs/*${R2_IDENTIFIER}*

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

parallel -k --link -j1 \
  'STAR \
    --runThreadN $CORES \
    --genomeDir genome/index \
    --readFilesIn "{1}" "{2}" \
    --outFileNamePrefix aligned/{=1s/_.*//=}_' \
  ::: processed/*${R1_IDENTIFIER}* \
  ::: processed/*${R2_IDENTIFIER}*

## Remove PCR duplicates and other poor reads.

if [ ! -d cleaned ]; then mkdir cleaned; fi

parallel -k -j1 \
  'samtools sort -n -@ $CORES "{}" | \
    samtools fixmate -m - - | \
    samtools sort -@ $CORES - | \
    samtools markdup - - | \
    samtools view -F 3852 -f 3 -O BAM -@ $CORES \
    -o cleaned/{=s/_.*//=}.bam' \
  ::: aligned/*Aligned*

## Index the BAMs.

parallel -j$CORES -k 'samtools index "{}"' ::: cleaned/*bam

## Multiqc report.

if [ ! -d multiqc ]; then mkdir multiqc; fi

multiqc -o multiqc -d .
