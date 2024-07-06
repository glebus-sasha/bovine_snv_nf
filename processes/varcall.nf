// Define the `VARCALL` process that performs variant calling
process VARCALL {
    container = 'staphb/bcftools:latest'
    tag "$reference $bamFile"
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/VARCALL"
//	debug true
//  errorStrategy 'ignore'
	
    input:
    path reference
    tuple val(sid), path(bai), path(bamFile)
    path fai

    output:
    tuple val(sid), path("${sid}.vcf"),         emit: tup_vcf
    path("${sid}.vcf"),                         emit: vcf

    script:
    """    
    bcftools mpileup -Ou -f $reference $bamFile | bcftools call -mv > ${sid}.vcf
    """
}
