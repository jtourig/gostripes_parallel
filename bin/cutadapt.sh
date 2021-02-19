#!/bin/bash

CUTADAPT_PAIRED() {
  cutadapt \
    -g ^NNNNNNNNTATAGGG \
    -j $CORES \
    -e 1 \
    --discard-untrimmed \
    -o fastqs/trimmed/$(basename ${1}) \
    -p fastqs/trimmed/$(basename ${2}) \
    ${1} ${2} \
    1> fastqs/trimmed/${3}_log.txt
}
export CUTADAPT_PAIRED

CUTADAPT_SINGLE() {
  cutadapt \
    -g ^NNNNNNNNTATAGGG \
    -j $CORES \
    -e 1 \
    --discard-untrimmed \
    -o fastqs/trimmed/$(basename ${1}) \
    ${1} \
    1> fastqs/trimmed/${2}_log.txt
}
export CUTADAPT_SINGLE
