#!/bin/bash

STAR_PAIRED() {
  STAR \
    --runThreadN $CORES \
    --genomeDir index \
    --readFilesIn fastqs/trimmed/$(basename ${1}) fastqs/trimmed/$(basename ${2}) \
    --outFileNamePrefix bams/aligned/${3}_
}
export STAR_PAIRED

STAR_SINGLE() {
  STAR \
    --runThreadN $CORES \
    --genomeDir index \
    --readFilesIn fastqs/trimmed/$(basename ${1}) \
    --outFileNamePrefix bams/aligned/${2}_
}
export STAR_SINGLE
