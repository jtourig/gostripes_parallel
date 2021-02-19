#!/bin/bash

CUTADAPT_PAIRED() {
  printf "%-5s%s" "..." "Filtering and trimming ${1} and ${2}"
  seqtk mergepe $1 $2 | \
  cutadapt \
    -g "^N{8}TATAG{3}" \
    -j $CORES \
    -O 14 \
    -e 1 \
    --discard-untrimmed \
    --interleaved \
    - | \
  cutadapt \
    -g "^CCTACACGACGCTCTTCCGATCTN{8}TATAG{3}" \
    -j $CORES \
    -O 21 \
    -e 2 \
    -m 15 \
    --interleaved \
    -o fastqs/trimmed/$(basename ${1}) \
    -p fastqs/trimmed/$(basename ${2}) \
    -
}

CUTADAPT_SINGLE() {
  printf "%s-5%s" "..." "Filtering and trimming ${1}"
  cutadapt \
    -g "^N{8}TATAG{3}" \
    -j $CORES \
    -e 1 \
    --discard-untrimmed \
    ${1} | \
  cutadapt \
    -g "^CCTACACGACGCTCTTCCGATCTN{8}TATAG{3}" \
    -j $CORES \
    -O 21 \
    -e 2 \
    -m 15 \
    -o fastqs/trimmed/$(basename ${1}) \
    -
}
