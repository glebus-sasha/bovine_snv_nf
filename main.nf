#!/usr/bin/env nextflow
// Include processes
include { QCONTROL }            from './processes/qcontrol.nf'
include { TRIM }                from './processes/trim.nf'
include { CUTADAPT }            from './processes/cutadapt.nf'
include { ALIGN }               from './processes/align.nf'
include { FLAGSTAT }            from './processes/flagstat.nf'
include { BAMINDEX }            from './processes/bamindex.nf'
include { VARCALL }             from './processes/varcall.nf'
include { WHATSHAP }            from './processes/whatshap.nf'
include { REPORT }              from './processes/report.nf'
include { VARCALL_MPILEUP }     from './processes/varcall_mpileup.nf'

// Logging pipeline information
log.info """\
\033[0;36m  ==========================================  \033[0m
\033[0;34m         B O V I N E   S N V   N F            \033[0m
\033[0;36m  ==========================================  \033[0m

    reference:  ${params.reference}
    reads:      ${params.reads}
    outdir:     ${params.outdir}
    """
    .stripIndent(true)

// Make the results directory if it needs
def result_dir = new File("${params.outdir}")
result_dir.mkdirs()

// Define the input channel for reference file
reference = params.reference ? Channel.fromPath("${params.reference}").collect(): null

// Define the input channel for FASTQ files, if provided
input_fastqs = params.reads ? Channel.fromFilePairs("${params.reads}/*[rR]{1,2}*.*{fastq,fq}*", checkIfExists: true) : null

// Define the input channel for bwa index files, if provided
bwaidx = params.bwaidx ? Channel.fromPath("${params.bwaidx}/*", checkIfExists: true).collect() : null

// Define the input channel for fai index files, if provided
faidx = params.bwaidx ? Channel.fromPath("${params.faidx}/*.fai", checkIfExists: true).collect() : null

// Define the bed_file channel
// If params.regions is provided, create a channel from the specified path and collect it into a list
// Otherwise, create a channel from the path "assets/dummy.bed" and collect it into a list
// https://github.com/nextflow-io/nextflow/issues/1694
bed_file = params.regions ? Channel.fromPath("${params.regions}").collect() : Channel.fromPath("assets/dummy.bed").collect()

// Define the workflow
workflow one{ 
    QCONTROL(input_fastqs)
    TRIM(input_fastqs)
    CUTADAPT(TRIM.out.trimmed_reads, "GCAG")
    ALIGN(CUTADAPT.out.cutadapted_reads, reference, bwaidx, bed_file)
    FLAGSTAT(ALIGN.out.bam)
    BAMINDEX(ALIGN.out.bam)
//    VARCALL(reference, BAMINDEX.out.bai, faidx, bed_file)
    VARCALL_MPILEUP(reference, BAMINDEX.out.bai, faidx, bed_file)
    WHATSHAP(reference, faidx, BAMINDEX.out.bai, VARCALL_MPILEUP.out.vcf)
    REPORT(TRIM.out.json.collect(), QCONTROL.out.zip.collect(), FLAGSTAT.out.flagstat.collect(), WHATSHAP.out.stats_tsv.collect())

    // Make the pipeline reports directory if it needs
    if ( params.reports ) {
        def pipeline_report_dir = new File("${params.outdir}/pipeline_info/")
        pipeline_report_dir.mkdirs()
    }
}

workflow two { 
    QCONTROL(input_fastqs)
    TRIM(input_fastqs)
    ALIGN(TRIM.out.trimmed_reads, reference, bwaidx, bed_file)
    FLAGSTAT(ALIGN.out.bam)
    BAMINDEX(ALIGN.out.bam)
    VARCALL_MPILEUP(reference, BAMINDEX.out.bai, faidx, bed_file)
    WHATSHAP(reference, faidx, BAMINDEX.out.bai, VARCALL_MPILEUP.out.vcf)
    REPORT(TRIM.out.json.collect(), QCONTROL.out.zip.collect(), FLAGSTAT.out.flagstat.collect(), WHATSHAP.out.stats_tsv.collect())

    // Make the pipeline reports directory if it needs
    if ( params.reports ) {
        def pipeline_report_dir = new File("${params.outdir}/pipeline_info/")
        pipeline_report_dir.mkdirs()
    }
}

workflow three {
    input_fastqs.view()

}


workflow {
    two()
}



