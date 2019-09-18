"""
ctDNA SNV caller comparison

Pipeline for assessing the ability of different variant callers to call low 
level SNVs for the use in ctDNA variant detection.

Created:       28 August 2019
Last updated:  10 September 2019
Author:        Erik Waskiewicz

Usage: nextflow run ctdna_snv_caller.nf

"""
version = "0.0.1"


/*-------------------------------------------------------------------------*
 *  Config
 *-------------------------------------------------------------------------*/

// get file objects for the fasta and bam files
genome_reference = [
    "fasta": file(params.ref_fasta), 
    "index": file(params.ref_fasta + ".fai"), 
    "dict":  file(params.ref_fasta.replace(".fasta", ".dict"))
]
bam_file = file("${params.indir}/${params.sample_name}.bam")


// print info to log and/or console
println ""
log.info """
Nextflow ctDNA SNV caller comparison  -  v$version
=============================================================================
inputs:
sample      $params.sample_name
ref_fasta   $genome_reference.fasta
bam         $bam_file """

println ""


/*-------------------------------------------------------------------------*
 *  Prepare input data
 *-------------------------------------------------------------------------*/

process preprocess_bam_rg {
    /* 
     * Adds reads groups to the input BAM file so it is compatable with 
     * downstream tools
     */
    publishDir "${params.outdir}/processed_bams", mode: "copy"

    input:
        file bam from bam_file

    output:
        file "${bam.name.replace('.bam', '_rg.bam')}" into bam_rg
        file "${bam.name.replace('.bam', '_rg.bai')}" into bam_rg_index

    script:
        """
        picard AddOrReplaceReadGroups \
            I=$bam \
            O=${bam.name.replace('.bam', '_rg.bam')} \
            RGID=RG1 \
            RGLB=LIBRARY1 \
            RGPL=ILLUMINA \
            RGPU=FLOWCELL1.LANE1 \
            RGSM=$params.sample_name \
            SORT_ORDER=coordinate \
            CREATE_INDEX=true
        """
}


process preprocess_bam_rmdup {
    /* 
     * 
     */
    publishDir "${params.outdir}/processed_bams", mode: "copy"

    input:
        file bam from bam_rg

    output:
        file "${bam.name.replace('.bam', '_rmdup.bam')}" into bam_rmdup
        file "${bam.name.replace('.bam', '_rmdup.bai')}" into bam_rmdup_index

    script:
        """
        picard MarkDuplicates \
            INPUT=$bam \
            OUTPUT=${bam.name.replace('.bam', '_rmdup.bam')} \
            METRICS_FILE=metrics.txt \
            CREATE_INDEX=true
        """
}


/*-------------------------------------------------------------------------*
 *  Variant calling 
 *-------------------------------------------------------------------------*/

process var_call_mutect {
    /* 
     * Call variants with GATK mutect2
     */
    container "broadinstitute/gatk:4.1.3.0"
    publishDir "${params.outdir}/vcfs", mode: "copy"

    input:
        // must input all files so that they can be seen in docker container
        file ref_file from genome_reference.fasta
        file ref_index from genome_reference.index
        file ref_dict from genome_reference.dict
        file bam from bam_rmdup
        file bam_index from bam_rmdup_index

    output:
        file "${params.sample_name}_mutect.vcf" into vcf_mutect

    script:
        """
        gatk Mutect2 \
            --input $bam \
            --output ${params.sample_name}_mutect.vcf \
            --reference $ref_file
        """
}


process var_call_sinvict {
    /* 
     * TODO - placeholder for other variant callers
     * add processes to call variants with other variant callers
     */
    container "erikwaskiewicz/ctdna-sinvict:latest"
    publishDir "${params.outdir}/vcfs", mode: "copy"

    input:
        file bam from bam_rmdup

    output:
        file "${params.sample_name}_sinvict.txt" into vcf_sinvict

    script:
        """
        bash /home/run_sinvict.sh $bam ${params.sample_name}_sinvict.txt
        """
}


/*-------------------------------------------------------------------------*
 *  Combine results and compare
 *-------------------------------------------------------------------------*/

process combine_vcfs {
    /* 
     * TODO
     * once all variant calling steps have been run, merge vcfs here 
     * and compare variant calls
     */
    publishDir "${params.outdir}/combined", mode: "copy"

    input:
        file f1 from vcf_mutect
        file f2 from vcf_sinvict

    output:
        file "${params.sample_name}_all.txt" into combined_vcfs

    script:
        """
        echo $f1 >> ${params.sample_name}_all.txt
        echo $f2 >> ${params.sample_name}_all.txt
        """
}


/*-------------------------------------------------------------------------*
 *  Pipeline completion
 *-------------------------------------------------------------------------*/

workflow.onComplete {

    log.info """
=============================================================================
Completed:  $workflow.complete
Status:     ${ workflow.success ? 'SUCCESS' : 'FAIL' } """

}
