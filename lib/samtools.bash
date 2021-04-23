#!/bin/bash

SAMTOOLS_PAIRED() {
  printf "%-5s%s" "..." "Processing ${1}"
  samtools sort -n -@ $CORES -m $MAXMEM bams/aligned/${1}_Aligned.out.sam | \
    samtools fixmate -m - - | \
    samtools sort -@ $CORES -m $MAXMEM - | \
    samtools markdup - - | \
    samtools view -F 3852 -f 3 -O BAM -@ $CORES \
    -o bams/cleaned/${1}.bam
}

SAMTOOLS_SINGLE() {

  if [[ ! -d tmpdir ]]; then mkdir tmpdir; fi

  printf "%-5s%s" "..." "Processing ${1}"
  samtools sort -@ $CORES -m $MAXMEM bams/aligned/${1}_Aligned.out.sam | \
    samtools view -F 2820 -O BAM -@ $CORES | \
    umi_tools dedup \
      --output-stats=bams/cleaned/${1}_dedup_stats.txt \
      --log=bams/cleaned/${1}_dedup_log.txt \
      --temp-dir=tmpdir \
      -S bams/cleaned/${1}.bam
}
