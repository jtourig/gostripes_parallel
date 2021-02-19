#!/bin/bash

SAMTOOLS_PAIRED() {
  printf "%-5s%s" "..." "Processing ${1}"
  samtools sort -n -@ $CORES bams/aligned/${1}_Aligned.out.sam | \
    samtools fixmate -m - - | \
    samtools sort -@ $CORES - | \
    samtools markdup - - | \
    samtools view -F 3852 -f 3 -O BAM -@ $CORES \
    -o bams/cleaned/${1}.bam
}
export SAMTOOLS_PAIRED

SAMTOOLS_SINGLE() {
  printf "%-5s%s" "..." "Processing ${1}"
  samtools sort -@ $CORES bams/aligned/${1}_Aligned.out.sam | \
    samtools view -F 2820 -O BAM -@ $CORES \
    -o bams/cleaned/${1}.bam
}
export SAMTOOLS_SINGLE
