// Define the `WHATSHAP` process that performs variant calling
process WHATSHAP {
    container = 'hangsuunc/whatshap:v1'
    tag "$reference $bamFile $bedfile"
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/WHATSHAP"

//    cache "lenient" 
//    debug true
    errorStrategy 'ignore'
	
    input:
    path reference
    path faidx
    tuple val(sid), path(bai), path(bamFile)
    tuple val(sid), path(vcf)

    output:
    tuple val("${sid}"), path("${sid}_phased.vcf"),      emit: phased_vcf
    tuple val("${sid}"), path("${sid}_phased.tsv"),      emit: stats_tsv
    
    script:
    """
    whatshap phase -o ${sid}_phased.vcf --reference=${reference} ${vcf} ${bamFile}
    whatshap stats --tsv ${sid}_phased.tsv ${sid}_phased.vcf
    """
}