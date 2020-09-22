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
    --run "test" \
    --runtype "test" \
    --reads "test_files/bams/sr/tiny/*.bam" \
    --output "${T_DIR}/results" \
    --seqmode "SE"\
    --macsconfig "test_files/test.config_sr" \
    --genomeIndex "tair10" \
    --threads 8 \
    -w ${T_DIR}/work \
    -resume


# check structure and content of output

if [ ! -d $T_DIR/results ]; then
    echo "missing results files"
    exit 1
fi

if [ ! $(ls -1q $T_DIR/results | wc -l) -eq 1 ] ; then
    echo "wrong number of results subdirs"
    exit 1
fi

if [ ! -d $T_DIR/results/test ]; then
    echo "missing results/test files"
    exit 1
fi

if [ ! $(ls -1q $T_DIR/results/test | wc -l) -eq 1 ] ; then
    echo "wrong number of results/test subdirs"
    exit 1
fi

if [ ! -d $T_DIR/results/test/main ]; then
    echo "missing results/test/main files"
    exit 1
fi

if [ ! $(ls -1q $T_DIR/results/test/main | wc -l) -eq 6 ] ; then
    echo "wrong number of results/test subdirs"
    exit 1
fi

if [ ! -d $T_DIR/results/test/main/aligned_files ] || \
    [ ! -d $T_DIR/results/test/main/bigwigs ] || \
    [ ! -d $T_DIR/results/test/main/deeptools ] || \
    [ ! -d ${T_DIR}'/results/test/main/tiny2_C5N3AANXX_7#22452_CGATGT' ] || \
    [ ! -d ${T_DIR}'/results/test/main/tiny2_C5N3AANXX_7#22453_TGACCA' ] || \
    [ ! -d ${T_DIR}'/results/test/main/tiny2_C5N3AANXX_7#22457_CTTGTA' ] ; then
    echo "missing results/test/main/.. files"
    exit 1
fi

if [ ! $(ls -1q $T_DIR/results/test/main/aligned_files | wc -l) -eq 9 ] ; then
    echo "wrong number of results/test/main/aligned_files files"
    exit 1
fi
if [ ! $(ls -1q $T_DIR/results/test/main/bigwigs | wc -l) -eq 4 ] ; then
    echo "wrong number of results/test/main/bigwigs files"
    exit 1
fi
if [ ! $(ls -1q $T_DIR/results/test/main/deeptools | wc -l) -eq 8 ] ; then
    echo "wrong number of results/test/main/deeptools file"
    exit 1
fi
if [ ! $(ls -1q $T_DIR/results/test/main/'tiny2_C5N3AANXX_7#22452_CGATGT' | wc -l) -eq 3 ] ; then
    echo "wrong number of results/test/main/'tiny2_C5N3AANXX_7#22452_CGATGT' files"
    exit 1
fi
if [ ! $(ls -1q $T_DIR/results/test/main/'tiny2_C5N3AANXX_7#22453_TGACCA' | wc -l) -eq 3 ] ; then
    echo "wrong number of results/test/main/'tiny2_C5N3AANXX_7#22453_TGACCA' files"
    exit 1
fi
if [ ! $(ls -1q $T_DIR/results/test/main/'tiny2_C5N3AANXX_7#22457_CTTGTA' | wc -l) -eq 3 ] ; then
    echo "wrong number of results/test/main/'tiny2_C5N3AANXX_7#22457_CTTGTA' files"
    exit 1
fi

if [ ! $(ls -1q $T_DIR/results/test/main/'tiny2_C5N3AANXX_7#22457_CTTGTA/statistics' | wc -l) -eq 2 ] ; then
    echo "wrong number of results/test/main/'tiny2_C5N3AANXX_7#22457_CTTGTA/statistics' files"
    exit 1
fi

if [ ! $(ls -1q $T_DIR/results/test/main/'tiny2_C5N3AANXX_7#22457_CTTGTA/trim_galore' | wc -l) -eq 3 ] ; then
    echo "wrong number of results/test/main/'tiny2_C5N3AANXX_7#22457_CTTGTA/trim_galore' files"
    exit 1
fi

if [ ! -f $T_DIR/results/test/main/'tiny2_C5N3AANXX_7#22457_CTTGTA/tiny2_C5N3AANXX_7#22457_CTTGTA.bowtie2.log' ]; then
    echo "bowtie2 log"
    exit 1
fi


## check that the files are correct, could add more

if [ ! -f $T_DIR/results/test/main/aligned_files/'tiny2_C5N3AANXX_7#22452_CGATGT.filtered_uniq.bam' ] || \
  [ $(md5sum $T_DIR/results/test/main/aligned_files/'tiny2_C5N3AANXX_7#22452_CGATGT.filtered_uniq.bam' | \
  awk '{print $1}') != "6478d0d369c33a59370b7e1d689b583c" ]; then
    echo "wrong or missing $T_DIR/results/test/main/aligned_files/'tiny2_C5N3AANXX_7#22452_CGATGT.filtered_uniq.bam'"
    exit 1
fi

if [ ! -f $T_DIR/results/test/main/bigwigs/testbam2.log2r.bw ] || \
  [ $(md5sum $T_DIR/results/test/main/bigwigs/testbam2.log2r.bw | \
  awk '{print $1}') != "f26d116bad4f0ddc1605c895dbde35c2" ]; then
    echo "wrong or missing $T_DIR/results/test/main/bigwigs/testbam2.log2r.bw"
    exit 1
fi


if [ ! -f $T_DIR/results/test/main/deeptools/matrix.chroms.gz ] ; then
    echo "missing $T_DIR/results/test/main/deeptools/matrix.chroms.gz"
    exit 1
fi
