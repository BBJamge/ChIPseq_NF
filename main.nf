#!/usr/bin/env nextflow

/**************
* Parameters
**************/


params.ref_seq = "tair10"
params.files  =   "bams/*.bam"
params.seqmode = "SR" // "PE"
params.extendReads = 200
params.bin = "10"
params.genomeSize = 119146348
params.outdir = "results"
log.info """\
outdir	: ${params.outdir}
"""

if( !(params.seqmode in ['PE','SR'] )) { exit 1, "Invalid sequencing mode: use --seqmode to choose either PE or SE" }


/**************
* Start
**************/

/**********************************
* read files of correct type
***********************************/

// Infer file type

if(params.seqmode == "PE"){
  if (params.files.split('/').takeRight(1)==['*.fastq'] | params.files.split('/').takeRight(1)==['*.fq']){
    exit 1. "Please end reads path with *_{1,2}.fq when using PE fastq files"
  }
  if (params.files.split('/').takeRight(1)==['*_{1,2}.fastq'] | params.files.split('/').takeRight(1)==['*_{1,2}.fq']){
    params.type = "fastq"
  }
}
if(params.files.split('/').takeRight(1)==['*.bam']){ //PE or SE
  params.type = "bam"
}
if(params.files.split('/').takeRight(1)==['*.fq'] | params.files.split('/').takeRight(1)==['*.fastq'] ){
  params.type = "fastq"
}
if( !params.type ){
  exit 1. "Type params not set"
}

/******************
* Create channels *
*******************/

if( params.seqmode == "PE" & params.type=="fastq"){
  reads = Channel
                  .fromFilePairs( params.files, size: -1 )
                  .ifEmpty { error "Invalid path: point to path containing fq files using --reads flag" }
}
else {
reads = Channel
                .fromPath( params.files )
                .ifEmpty { error "Invalid path: point to path containing BAM/fq files using --reads flag" }
                .map { file -> tuple(file.baseName, file) }
}
if( params.type == "fastq"){
  reads.into{reads;fq_files} // put into read_pairs later
}








//Set parameters
if ( params.ref_seq == "tair10"){
  params.index="library://elin.axelsson/index/index_bowtie2_tair10:v2.4.1-release-47"
}


// index
process get_bowtie2_index {

  input:
  file params.index

  output:
  file params.ref_seq into bw2index

  script:
  """
  singularity run ${params.index}
  """
}


/*************
* processes for BAMtoFASTQ
*************/

process Bam2Fastq {
    tag "$id"

    input:
    set id, file(bam) from reads

    output:
    set val(id), file("${id}*fastq") into bam_fastq

    when:
    params.type=="bam"

    script:
    if (params.seqmode == 'SR')
      """
      bedtools bamtofastq -i ${bam} -fq ${id}.fastq
      """
    else if (params.seqmode == 'PE')
      """
      bedtools bamtofastq -i ${bam} -fq ${id}_1.fastq -fq2 ${id}_2.fastq
      """
    }

/***********************
* COPY FASTQ FILES TO CHANNEL(if fastq files provided)
***********************/
  if(params.type=="fastq"){
    read_pairs = read_pairs.mix(fq_files)
  }

/**************
* Split bam_fastq into 4 for 4 different process
**************/

bam_fastq.into{bam_fastq; bam_qc; bam_for_stats; bam_for_trim}

/**************
* QC
**************/
process FASTQC {
    publishDir "$params.outdir/QUALfiltered",mode:'copy'

    input:
    set id, file(input) from bam_qc

    output:
    file "*.{zip,html,txt}" into fastqc_results

    script:
    """
    fastqc ${input}
    """
}

process FilterFastq {
    publishDir "$params.outdir/FilterFastq",mode:'copy'

    input:
    set datasetID, file(reads) from bam_fastq


    output:
    set datasetID, file(reads) into filtered_fastq
    file("${datasetID}*trimming_report.txt")
    file("${datasetID}*fastqc.{html,zip}")

    script:
    if (params.seqmode == 'PE')
        """
        trim_galore --dont_gzip --stringency 1 --fastqc --length 5 --paired ${reads}
        """
    else if (params.seqmode == 'SR')
        """
        trim_galore --dont_gzip --stringency 1 --fastqc --length 5 ${reads}
        """
    else
    error "Invalid sequencing mode: choose either PE or SR"
}




process BamStats {
    publishDir "$params.outdir/QUALfiltered",mode:'copy'
    label 'env_qc_small'

    input:
    set id, file(x) from bam_for_stats

    output:
    set id, file("${id}_qualstats.txt") into bam_stats

    script:
    """
    fastx_quality_stats -i ${x} -o ${id}_qualstats.txt  -Q33
    """
}

process TrimFastq {
    publishDir "$params.outdir/QUALfiltered",mode:'copy'
    label 'env_qc_small'

    input:
    set id, file(x) from bam_for_trim

    output:
    set id, file("${id}_qualtrimmed.fastq") into bam_trim

    script:
    """
	  seqtk trimfq ${x} > ${id}_qualtrimmed.fastq
    """
}


process AlignBowtie2 {
    //publishDir "$params.outdir/Aligned_BAMS",mode:'copy'
    label 'env_align_medium'

    input:
    set id, file(x) from filtered_fastq
    file(genomeIndex) from bw2index

    output:
    set id, file("${id}.sam") into outsam

    script:
    if (params.seqmode == 'PE')
        """
        bowtie2 -x ${genomeIndex}/${genomeIndex} -1 ${x[0]} -2 ${x[0]} -S ${id}.sam
        """
    else if (params.seqmode == 'SR')
        """
        bowtie2 -x ${genomeIndex}/${genomeIndex} -U ${x} -S ${id}.sam
        """
}

process Sam2Bam {
    label 'env_sam_medium'
    publishDir "$params.outdir/Aligned_BAMS",mode:'copy'

    input:
    set id, file(x) from outsam

    output:
    set id, file("${id}.sorted.bam"), file("${id}.sorted.bam.bai") into alignedbam
//    set id2, file("${id2}.sorted.bam.bai") into alignedbambai

    script:
    """
    samtools view -S -b ${x} > ${id}.bam
    samtools sort ${id}.bam > ${id}.sorted.bam
    samtools index ${id}.sorted.bam
    """
}


process rmdup {
    label 'env_picard_medium'
    publishDir "$params.outdir/Aligned_BAMS/Dedup",mode:'copy'

    input:
    set id, file(x), file(bai) from alignedbam

    output:
    set id, file("${id}.dedup.bam") into dedupbam


    script:
    """
    java -jar /usr/picard/picard.jar MarkDuplicates \
     I=${x} \
     O=${id}.dedup.bam \
     M=${id}.dedup.metrics.txt \
     ASSUME_SORTED=true \
     REMOVE_DUPLICATES=true
    """
}

dedupbam.into{dedupbam;dedupbam_index}

process index_rmdup {

  input:
  set id, file(bam) from dedupbam_index

  output:
  file("${id}.dedup.sorted.bam") into dedupsortbam
  file("${id}.dedup.sorted.bam.bai") into dedupbambai

  script:
  """
  samtools sort ${id}.dedup.bam > ${id}.dedup.sorted.bam
  samtools index ${id}.dedup.sorted.bam
  """
  }

process corPlot {
    publishDir "$params.outdir/QC_plots",mode:'copy'

    input:
    file (allbams) from dedupsortbam.toSortedList()
    file (allbais) from dedupbambai.toSortedList()

    output:
    file("QCall*") into cor_npz

    script:
    """
    multiBamSummary bins --bamfiles ${allbams} --ignoreDuplicates -o QCall.npz
    """
}

//dedupsortbam.collect().toSortedList().println()



workflow.onComplete {
	println ( workflow.success ? "Successfull!" : "Messed Up something" )
}
