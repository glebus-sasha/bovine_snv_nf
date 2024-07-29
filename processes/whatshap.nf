// Define the `WHATSHAP` process that performs variant calling
process WHATSHAP {
    container = 'angsuunc/whatshap'
    tag "$reference $bamFile $bedfile"
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/VARCALL"

//    cache "lenient" 
//    debug true
    errorStrategy 'ignore'
	
    input:
    path reference
    tuple val(sid), path(bai), path(bamFile)
    path val(sid), path(vcf)

    output:
    tuple val("${sid}"), path("${sid}_phased.vcf"),      emit: phased_vcf
    
    script:
    """
    whatshap phase -o ${sid}_phased.vcf --reference=${reference} ${vcf} ${bamFile}
    # whatshap genotype --reference ${reference} -o ${sid}_genotyped.vcf ${sid}_phased.vcf ${bamFile}
    # whatshap stats --tsv ${sid}_genotyped.vcf
    """
}