#!/usr/bin/env nextflow

/**************
* Parameters
**************/

params.ref_seq = "/groups/berger/lab/cluster_files/Bhagyshree/TAIR10/Sequence/Bowtie2Index/genome"
//params.ref_id = "genome"
params.files  =   "bams/*.bam"
//params.tmpdir = "$TMPDIR"
params.seqtkdir = "/groups/berger/lab/cluster_files/Bhagyshree/seqtk"
params.extendReads = 200
params.bin = "10"
params.genomeSize = 119146348
params.outdir = "results"
log.info """\
outdir	: ${params.outdir}
"""

//bowtie2_index = Channel.fromPath( params.ref_fol + "/" + params.ref_id + "*" )

/**************
* Start
**************/
// first put all bam files into channel "bamfiles"

bamfiles = Channel
  .fromPath(params.files)
  .map { file -> [ file.baseName, file] }
// each item in the channel has a file name and an actual file


/*************
* processes for BAMtoFASTQ
*************/

//bamfiles.combine(Channel.from(params.bin.split(","))).set {to_fastq}

process Bam2Fastq {
    label 'env_bed_small'
    tag "$id"
    //publishDir "$params.outdir/",mode:'copy'

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
    label 'env_qc_small'

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
    label 'env_qc_small'

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
	${params.seqtkdir}/seqtk trimfq ${x} > ${id}_qualtrimmed.fastq
    """
}

process AlignBowtie2 {
    //publishDir "$params.outdir/Aligned_BAMS",mode:'copy'
    label 'env_align_medium'

    input:
    set id, file(x) from filtered_fastq
    //file(indices) from bowtie2_index

    output:
    set id, file("${id}.sam") into outsam

    script:
    """
    bowtie2 -x "$params.ref_seq" -U ${x} -S ${id}.sam
    """
}

process Sam2Bam {
    label 'env_sam_medium'
    publishDir "$params.outdir/Aligned_BAMS",mode:'copy'

    input:
    set id, file(x) from outsam

    output:
    set id, file("${id}.sorted.bam*") into alignedbam

    script:
    """
    samtools view -S -b ${x} > ${id}.bam
    samtools sort ${id}.bam > ${id}.sorted.bam
    samtools index ${id}.sorted.bam
    """
}

workflow.onComplete {
	println ( workflow.success ? "Successfull!" : "Messed Up something" )
}
