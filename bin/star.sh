#!/bin/bash

STAR_PAIRED() {
  printf "%-5s%s" "..." "Aligning ${1} and ${2}"
  STAR \
    --runThreadN $CORES \
    --genomeDir index \
    --readFilesIn fastqs/trimmed/$(basename ${1}) fastqs/trimmed/$(basename ${2}) \
    --outFileNamePrefix bams/aligned/${3}_
}
export STAR_PAIRED

STAR_SINGLE() {
  printf "%-5s%s" "..." "Aligning ${1}"
  STAR \
    --runThreadN $CORES \
    --genomeDir index \
    --readFilesIn fastqs/trimmed/$(basename ${1}) \
    --outFileNamePrefix bams/aligned/${2}_
}
export STAR_SINGLE
