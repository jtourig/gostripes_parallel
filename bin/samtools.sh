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

SAMTOOLS_SINGLE() {
  printf "%-5s%s" "..." "Processing ${1}"
  samtools sort -@ $CORES bams/aligned/${1}_Aligned.out.sam | \
    samtools view -F 2820 -O BAM -@ $CORES | \
    umi_tools dedup \
      --output-stats=bams/cleaned/${1}_dedup_stats.txt \
      --log=bams/cleaned/${1}_dedup_log.txt \
      -S bams/cleaned/${1}.bam
}
