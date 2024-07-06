// Define the `ANNOTATE` process that performs annotation
process ANNOTATE {
    container = 'ensemblorg/ensembl-vep:latest'
    containerOptions "-B ${params.vepcache}:/opt/vep/.vep"
    tag "$vcf"
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/ANNOTATE"
    cpus params.cpus
    debug true
    cache "lenient"
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(vcf)

    output:
    path '*.vep', emit: vep
    path '*.vep.html', emit: html

    script:
    """
    vep \
    -i $vcf \
    -o ${sid}.vep \
    --stats_file ${sid}.vep.html \
    --fork ${task.cpus} \
    --database \
    --everything \
    --species bos_taurus 
    """
}
