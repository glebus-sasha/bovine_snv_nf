// Define the `VARCALL_MPILEUP` process that performs variant calling
process VARCALL_MPILEUP {
    container = 'staphb/bcftools:latest'
    tag "$reference $bamFile"
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/VARCALL_MPILEUP"
//	debug true
  errorStrategy 'ignore'
	
    input:
    path reference
    tuple val(sid), path(bai), path(bamFile)
    path fai
    path regions

    output:
    tuple val("${sid}"), path("${sid}.vcf"),         emit: vcf

    script:
    """    
    bcftools mpileup -R ${regions} -f $reference $bamFile -Ou >${sid}.vcf
    #| bcftools call -mv -Ov -o ${sid}.vcf
    """
}
