"""
ctDNA SNV caller comparison

Pipeline for assessing the ability of different variant callers to call
low level SNVs for the use in ctDNA variant detection.

Usage: nextflow run ctdna_snv_caller.nf

"""
version = "0.0.2"


/*-------------------------------------------------------------------*
 *  Config
 *-------------------------------------------------------------------*/

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
=======================================================================
inputs:
sample      $params.sample_name
ref_fasta   $genome_reference.fasta
bam         $bam_file """

println ""


/*-------------------------------------------------------------------*
 *  Prepare input data
 *-------------------------------------------------------------------*/

process preprocess_bam_rg {
    /* 
     * Adds reads groups to the input BAM file so it is compatable 
     * with downstream tools. Also sorts and indexes the BAM file.
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
     * Marks and removes any PCR duplicates in the BAM file and creates
     * an index.
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


/*-------------------------------------------------------------------*
 *  Variant calling 
 *-------------------------------------------------------------------*/

process var_call_mutect {
    /* 
     * Call variants with GATK mutect2
     */
    container "broadinstitute/gatk:4.1.3.0"
    publishDir "${params.outdir}/vcfs", mode: "copy"

    input:
        // must input all files so they can be seen in docker container
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
     * Call variants with the SiNVICT variant caller
     * Runs a custom script within Docker dontainer that:
     *   - runs bam-readcount to make pileup file
     *   - runs SiNVICT (outputs 6 text files, one for each filter)
     *   - converts output into VCF file
     */
    container "erikwaskiewicz/ctdna-sinvict:latest"
    publishDir "${params.outdir}/vcfs", mode: "copy"

    input:
        file ref_file from genome_reference.fasta
        file bam from bam_rmdup

    output:
        file "${params.sample_name}_sinvict.vcf" into vcf_sinvict

    script:
        """
        python /home/run_sinvict.py \
            $bam \
            $ref_file \
            ${params.sample_name}_sinvict.vcf
        """
}


process var_call_varscan {
    /* 
     * Call variants with VarScan2
     */
    container "erikwaskiewicz/ctdna-varscan:latest"
    publishDir "${params.outdir}/vcfs", mode: "copy"

    input:
        file ref_file from genome_reference.fasta
        file bam from bam_rmdup

    output:
        file "${params.sample_name}_varscan.vcf" into vcf_varscan

    script:
        """
        # make text file with sample name, as required by varscan
        echo $params.sample_name > sample_list.txt
        
        samtools mpileup -B -f $ref_file $bam | \
        varscan mpileup2snp \
            --output-vcf 1 \
            --vcf-sample-list sample_list.txt \
            > ${params.sample_name}_varscan.vcf
        """
}


/*-------------------------------------------------------------------*
 *  Combine results and compare
 *-------------------------------------------------------------------*/

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
        file f3 from vcf_varscan

    output:
        file "${params.sample_name}_all.txt" into combined_vcfs

    script:
        """
        echo $f1 >> ${params.sample_name}_all.txt
        echo $f2 >> ${params.sample_name}_all.txt
        echo $f3 >> ${params.sample_name}_all.txt
        """
}


/*-------------------------------------------------------------------*
 *  Pipeline completion
 *-------------------------------------------------------------------*/

workflow.onComplete {

    log.info """
=======================================================================
Completed:  $workflow.complete
Status:     ${ workflow.success ? 'SUCCESS' : 'FAIL' } """

}
