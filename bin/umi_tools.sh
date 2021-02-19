#!/bin/bash

UMI_TOOLS_EXTRACT() {
  umi_tools extract \
    -p "NNNNNNNN" \
    -I $1 \
    -L fastqs/trimmed/${2}_umi_log.txt
}
