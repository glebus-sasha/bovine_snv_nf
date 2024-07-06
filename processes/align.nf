// Define the `ALIGN` process that aligns reads to the reference genome
process ALIGN {
    container = 'glebusasha/bwa_samtools'
    tag "$reference ${sid} $bedfile"
    cpus params.cpus
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/ALIGN"
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(fq_1_trimmed), path(fq_2_trimmed), path(bedfile)
    path reference
    path idx
    
    output:
    tuple val("${sid}_${bedfile.baseName}"), path("*.sorted.bam"), emit: bam
    
    script:
    def bed_option = bedfile.getBaseName() == 'dummy' ? "" : "-L ${bedfile}"    // If the base name of bedfile is 'dummy', set bed_option to an empty string
    """
        bwa mem \
            -t ${task.cpus} ${reference} ${fq_1_trimmed} ${fq_2_trimmed} | \
        samtools view -bh ${bed_option} | \
        samtools sort -o ${sid}_${bedfile.baseName}.sorted.bam

    """
}