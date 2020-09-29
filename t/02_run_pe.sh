#!/bin/bash
set -o errexit
set -o pipefail

################################################################
## check that the pipeline runs and files end up in results
################################################################

# Single read

T_DIR="test_02"
mkdir -p ${T_DIR}

nextflow run ../main.nf \
    -profile "test_local" \
    --seqmode "PE" \
    --files "test_files/bams/pe/*.bam" \
    --outdir "${T_DIR}/results/" \
    -w ${T_DIR}/work \
    -resume


    # check structure and content of output

    # check structure and content of output

    if [ ! -d $T_DIR/results ]; then
        echo "missing results files"
        exit 1
    fi

    if [ ! $(ls -1q $T_DIR/results | wc -l) -eq 4 ] ; then ## double check
        echo "wrong number of results subdirs"
        exit 1
    fi

    if [ ! -d $T_DIR/results/Aligned_BAMS ]; then
        echo "missing results/Aligned_BAMS files"
        exit 1
    fi

    if [ ! -d $T_DIR/results/Aligned_BAMS/Dedup ]; then
        echo "missing results/Aligned_BAMS/Dedup files"
        exit 1
    fi

    if [ ! -d $T_DIR/results/FilterFastq ]; then
        echo "missing results/FilterFastq files"
        exit 1
    fi

    if [ ! -d $T_DIR/results/QC_plots ]; then
        echo "missing results/QC_plots files"
        exit 1
    fi

    if [ ! -d $T_DIR/results/QUALfiltered ]; then
        echo "missing results/QUALfiltered files"
        exit 1
    fi


    if [ ! $(ls -1q $T_DIR/results/Aligned_BAMS | wc -l) -eq 7 ] ; then ## 3*2 + Dedup
        echo "wrong number of results/Aligned_BAMS"
        exit 1
    fi

    if [ ! $(ls -1q $T_DIR/results/Aligned_BAMS/Dedup | wc -l) -eq 3 ] ; then
        echo "wrong number of results/Aligned_BAMS/Dedup"
        exit 1
    fi

    if [ ! $(ls -1q $T_DIR/results/FilterFastq | wc -l) -eq 24 ] ; then ## 3*4*2
        echo "wrong number of results/FilterFastq"
        exit 1
    fi

    if [ ! $(ls -1q $T_DIR/results/QC_plots | wc -l) -eq 5 ] ; then 
        echo "wrong number of results/QC_plots"
        exit 1
    fi

    if [ ! $(ls -1q $T_DIR/results/QUALfiltered  | wc -l) -eq 18 ] ; then ## 6*3
        echo "wrong number of results/QUALfiltered"
        exit 1
    fi


## And more on file level
