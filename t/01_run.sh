#!/bin/bash
set -o errexit
set -o pipefail

################################################################
## check that the pipeline runs and files end up in results
################################################################

# Single read

T_DIR="test_01"
mkdir -p ${T_DIR}

nextflow run ../main.nf \
    -profile "test_local" \
    --files "test_files/bams/sr/tiny/*.bam" \
    -resume 
