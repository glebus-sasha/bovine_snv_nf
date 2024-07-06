// Define the `BSFSTATS` process that performs variant calling
process BSFSTATS {
    container = 'staphb/bcftools:latest'
    tag "$vcf"
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/BSFSTATS"
//	debug true
//  errorStrategy 'ignore'
	
    input:
    tuple val(sid), path(vcf)
    
    output:
    path("${sid}.bcfstats"),      emit: bcfstats

    script:
    """    
    bcftools stats $vcf > ${sid}.bcfstats
    """
}
