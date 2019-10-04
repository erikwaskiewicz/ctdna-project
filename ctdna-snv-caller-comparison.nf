"""
ctDNA SNV caller comparison

Pipeline for assessing the ability of different variant callers to call
low level SNVs for the use in ctDNA variant detection.

"""

help_message = """
===============================================================================

ctdna-snv-caller-comparison.nf  -  v${params.version}

===============================================================================

Pipeline for assessing the ability of different variant callers to call
low level SNVs for the use in ctDNA variant detection.

Usage: nextflow ctdna-snv-caller-comparison.nf [ parameters ]

Required parameters:
  --bam_file        Absolute path to BAM file
  --fasta_ref       Absolute path to FASTA file

Output:
  --outdir          Path where the results to be saved [Default: './data/output']

"""

/*-------------------------------------------------------------------*
 *  Config
 *-------------------------------------------------------------------*/

// show help text if help flag used
params.help = false
if ( params.help ) { log.info help_message; exit 0 }

// check for input files
if ( !params.input_bam )    { log.error "ERROR\tinput bam is required"; exit 1 }
if ( !params.genome_fasta ) { log.error "ERROR\tgenome fasta is required"; exit 1 }

// check if fasta indexes and dict exists
if ( file(params.genome_fasta + ".fai").exists() ) {
    params.genome_fasta_index = params.genome_fasta + ".fai"
} else { params.genome_fasta_index = false }

if ( file(params.genome_fasta.replace(".fasta", ".dict")).exists() ) {
    params.genome_fasta_dict = params.genome_fasta.replace(".fasta", ".dict")
} else { params.genome_fasta_dict = false }

// make file objects
bam_file               = file("$params.input_bam", checkIfExists: true)
genome_fasta_file      = file("$params.genome_fasta", checkIfExists: true)
params.sample_name     = bam_file.name.replace(".bam", "")

// Create a summary for the logfile
def summary = [:]
summary["Version"]     = "$params.version"
summary["Git info"]    = "$workflow.repository - $workflow.revision [$workflow.commitId]"
summary["Command"]     = "$workflow.commandLine"
summary["Start time"]  = "$workflow.start"
summary["Sample"]      = "$params.sample_name"
summary["Project"]     = "$workflow.projectDir"
summary["BAM file"]    = "$params.input_bam"
summary["fasta file"]  = "$params.genome_fasta"
summary["fasta index"] = "$params.genome_fasta_index"
summary["fasta dict"]  = "$params.genome_fasta_dict"

// print info to log and/or console
log.info """
===============================================================================
Inputs:\n${ summary.collect  { k,v -> "  ${k.padRight(12)}: $v" }.join("\n") }
===============================================================================
"""


/*-------------------------------------------------------------------*
 *  Prepare input data
 *-------------------------------------------------------------------*/

// make fasta index if not supplied
if ( !params.genome_fasta_index ) {

    process preprocess_fasta_index {
        publishDir "$baseDir/data/test_data/input/", mode: "copy"

        input:
            file(fasta) from genome_fasta_file

        output:
            file("${fasta.name + '.fai'}") into genome_fasta_index

        script:
        """
        samtools faidx $fasta > ${fasta.name + '.fai'}
        """
    }

} else { genome_fasta_index = file("$params.genome_fasta_index") }


// make fasta dict if not supplied
if ( !params.genome_fasta_dict ) {

    process preprocess_fasta_dict {
        publishDir "$baseDir/data/test_data/input/", mode: "copy"

        input:
            file(fasta) from genome_fasta_file

        output:
            file("${fasta.name.replace(".fasta", ".dict")}") into genome_fasta_dict

        script:
        """
        picard CreateSequenceDictionary \
            R=$fasta \
            O=${fasta.name.replace(".fasta", ".dict")}
        """
    }

} else { genome_fasta_dict = file("$params.genome_fasta_dict") }


process preprocess_bam_rg {
    /* 
     * Adds reads groups to the input BAM file so it is compatable 
     * with downstream tools. Also sorts and indexes the BAM file.
     */
    publishDir "${params.outdir}/processed_bams", mode: "copy"

    input:
        file(bam) from bam_file

    output:
        file("${bam.name.replace('.bam', '_rg.bam')}") into bam_rg
        file("${bam.name.replace('.bam', '_rg.bai')}") into bam_rg_index

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
        file(bam) from bam_rg

    output:
        file("${bam.name.replace('.bam', '_rmdup.bam')}") into bam_rmdup
        file("${bam.name.replace('.bam', '_rmdup.bai')}") into bam_rmdup_index

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
        file(ref_file) from genome_fasta_file
        file(ref_index) from genome_fasta_index
        file(ref_dict) from genome_fasta_dict
        file(bam) from bam_rmdup
        file(bam_index) from bam_rmdup_index

    output:
        file("${params.sample_name}_mutect.vcf") into vcf_mutect

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
        file(ref_file) from genome_fasta_file
        file(bam) from bam_rmdup

    output:
        file("${params.sample_name}_sinvict.vcf") into vcf_sinvict

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
        file(ref_file) from genome_fasta_file
        file(bam) from bam_rmdup

    output:
        file("${params.sample_name}_varscan.vcf") into vcf_varscan

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
        file(f1) from vcf_mutect
        file(f2) from vcf_sinvict
        file(f3) from vcf_varscan

    output:
        file("${params.sample_name}_all.txt") into combined_vcfs

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
        ===============================================================================
        Completed:  $workflow.complete
        Status:     ${ workflow.success ? 'SUCCESS' : 'FAIL' }
        """.stripIndent()

}
