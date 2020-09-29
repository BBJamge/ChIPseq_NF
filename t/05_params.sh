#!/bin/bash
set -o errexit
set -o pipefail

################################################################
## check in parameters
################################################################


T_DIR="test_05"
mkdir -p ${T_DIR}


echo "seqmode"

nextflow run ../main.nf \
    -profile "test_local" \
    --seqmode "NE" \
    --files "test_files/fastq/pe/*_{1,2}.fq" \
    --outdir "${T_DIR}/results/" \
    -w ${T_DIR}/work \
    -resume \
     > ${T_DIR}/tmp && exit 1


if [ ! $(grep "Invalid sequencing mode: use --seqmode to choose either PE or SE"   ${T_DIR}/tmp | wc -l ) -gt 0 ]; then
    echo "wrong error"
    exit 1
fi


rm ${T_DIR}/tmp

echo "reads given, exists"


nextflow run ../main.nf \
    -profile "test_local" \
    --seqmode "SR" \
    --files "empty/*.bam" \
    --outdir "${T_DIR}/results/" \
    -w ${T_DIR}/work \
    -resume \
     > ${T_DIR}/tmp && exit 1

if [ ! $(grep "Invalid path: point to path containing BAM/fq files using --reads flag"   ${T_DIR}/tmp | wc -l ) -gt 0 ]; then
      echo "wrong error"
      exit 1
fi

rm ${T_DIR}/tmp


nextflow run ../main.nf \
    -profile "test_local" \
    --seqmode "PE" \
    --files "empty/*_{1,2}.fq" \
    --outdir "${T_DIR}/results/" \
    -w ${T_DIR}/work \
    -resume \
     > ${T_DIR}/tmp && exit 1


if [ ! $(grep "Invalid path: point to path containing fq files using --reads flag"   ${T_DIR}/tmp | wc -l ) -gt 0 ]; then
      echo "wrong error"
      exit 1
fi

rm ${T_DIR}/tmp

# later genome etc

echo "match reads and fastq"

nextflow run ../main.nf \
    -profile "test_local" \
    --seqmode "PE" \
    --files "test_files/fastq/pe/*.fq" \
    --outdir "${T_DIR}/results/" \
    -w ${T_DIR}/work \
    -resume \
     > ${T_DIR}/tmp && exit 1


if [ ! $(grep "Please end reads path with \*\_{1,2}.fq when using PE fastq files"  ${T_DIR}/tmp | wc -l ) -gt 0 ]; then
    echo "wrong error"
    exit 1
fi

rm ${T_DIR}/tmp
