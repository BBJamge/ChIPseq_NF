#!/usr/bin/env nextflow

/**************
* Parameters
**************/


params.ref_seq = "tair10"
params.files  =   "bams/*.bam"
params.extendReads = 200
params.bin = "10"
params.genomeSize = 119146348
params.outdir = "results"
log.info """\
outdir	: ${params.outdir}
"""



/**************
* Start
**************/

//Set parameters
if ( params.ref_seq == "tair10"){
  params.index="library://elin.axelsson/index/index_bowtie2_tair10:v2.4.1-release-47"
}


// first put all bam files into channel "bamfiles"

bamfiles = Channel
  .fromPath(params.files)
  .map { file -> [ file.baseName, file] }
// each item in the channel has a file name and an actual file


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
    set id, file(bam) from bamfiles

    output:
    set val(id), file('*.fastq') into bam_fastq

    script:
    """
    bedtools bamtofastq -i ${bam} -fq ${id}.fastq
    """
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
    set id, file(input) from bam_fastq

    output:
    set id, file("${id}_filtered.fastq") into filtered_fastq
    file("${id}_filtered.fastq")

    script:
    """
    fastq_quality_filter -q 30 -p 90 -Q33 -i ${input} -v -o ${id}_filtered.fastq
    """
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
