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
if ( !params.roi_bed )      { log.error "ERROR\tROI BED file is required"; exit 1 }

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
roi_bed_file           = file("$params.roi_bed", checkIfExists: true)
params.sample_name     = bam_file.name.replace(".bam", "")

// Create a summary for the logfile
def summary = [:]
summary["Version"]     = "$params.version"
summary["Git info"]    = "$workflow.repository - $workflow.revision [$workflow.commitId]"
summary["Command"]     = "$workflow.commandLine"
summary["Start time"]  = "$workflow.start"
summary["Sample"]      = "$params.sample_name"
summary["Project"]     = "$workflow.projectDir"
summary["ROI BED file"]= "$params.roi_bed"
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
    publishDir "${params.outdir}/vcfs", mode: "copy"

    input:
        // must input all files so they can be seen in docker container
        file(ref_file) from genome_fasta_file
        file(ref_index) from genome_fasta_index
        file(ref_dict) from genome_fasta_dict
        file(bam) from bam_rmdup
        file(bam_index) from bam_rmdup_index
        file(bed_file) from roi_bed_file

    output:
        file("${params.sample_name}_mutect.vcf") into vcf_mutect

    script:
        """
        gatk Mutect2 \
            --input $bam \
            --output ${params.sample_name}_mutect.vcf \
            --reference $ref_file \
            --intervals $bed_file
        """
}


process process_mutect {
    /* 
     * Convert Mutect VCF into table
     */
    input:
        file(vcf) from vcf_mutect

    output:
        file("${params.sample_name}_mutect.txt") into processed_mutect

    script:
        """
        # convert to table
        gatk VariantsToTable \
            -V $vcf \
            -O variants.txt \
            -F CHROM -F POS -F REF -F ALT -F FILTER \
            -GF AF -GF AD -GF DP
        
        # merge chr,pos,ref and alt, change headers and save to file
        awk ' { print \$1":"\$2\$3">"\$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8 } ' variants.txt | \
        sed "1 s/^.*\$/variant\tmutect_filter\tmutect_AF\tmutect_AD\tmutect_DP/" \
            > ${params.sample_name}_mutect.txt
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
    container "${params.singularity_cache}/ctdna-sinvict-latest.simg"
    publishDir "${params.outdir}/vcfs", mode: "copy"

    input:
        file(ref_file) from genome_fasta_file
        file(bam) from bam_rmdup
        file(bam_index) from bam_rmdup_index
        file(bed_file) from roi_bed_file

    output:
        file("${params.sample_name}_sinvict.vcf") into vcf_sinvict

    script:
        """
        python ${workflow.projectDir}/env/docker_images/ctdna-sinvict/run_sinvict.py \
            $bam \
            $ref_file \
            ${params.sample_name}_sinvict.vcf \
            $bed_file
        """
}


process process_sinvict {
    /* 
     * Convert SiNVICT VCF into table
     */
    input:
        file(vcf) from vcf_sinvict

    output:
        file("${params.sample_name}_sinvict.txt") into processed_sinvict

    script:
        """
        # convert to table
        gatk VariantsToTable \
            -V $vcf \
            -O variants.txt \
            -F CHROM -F POS -F REF -F ALT -F FILTER \
            -GF SVVAF -GF SVR -GF DP
        
        # merge chr,pos,ref and alt, change headers and save to file
        awk ' { print \$1":"\$2\$3">"\$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8 } ' variants.txt | \
        sed "1 s/^.*\$/variant\tsinvict_filter\tsinvict_AF\tsinvict_AD\tsinvict_DP/" \
            > ${params.sample_name}_sinvict.txt
        """
}


process var_call_varscan {
    /* 
     * Call variants with VarScan2
     */
    container "${params.singularity_cache}/ctdna-varscan-latest.simg"
    publishDir "${params.outdir}/vcfs", mode: "copy"

    input:
        file(ref_file) from genome_fasta_file
        file(bam) from bam_rmdup
        file(bed_file) from roi_bed_file

    output:
        file("${params.sample_name}_varscan.vcf") into vcf_varscan

    script:
        """
        # make text file with sample name, as required by varscan
        echo $params.sample_name > sample_list.txt
        
        samtools mpileup \
            --fasta-ref $ref_file \
            --no-BAQ \
            --positions $bed_file \
            $bam | \
        varscan mpileup2snp \
            --output-vcf 1 \
            --vcf-sample-list sample_list.txt \
            > ${params.sample_name}_varscan.vcf
        """
}


process process_varscan {
    /* 
     * Convert VarScan VCF into table
     */
    input:
        file(vcf) from vcf_varscan

    output:
        file("${params.sample_name}_varscan.txt") into processed_varscan

    script:
        """
        # convert to table
        gatk VariantsToTable \
            -V $vcf \
            -O variants.txt \
            -F CHROM -F POS -F REF -F ALT -F FILTER \
            -GF FREQ -GF AD -GF DP
        
        # merge chr,pos,ref and alt, change headers and save to file
        awk ' { print \$1":"\$2\$3">"\$4"\t"\$5"\t"\$6"\t"\$7"\t"\$8 } ' variants.txt | \
        sed "1 s/^.*\$/variant\tvarscan_filter\tvarscan_AF\tvarscan_AD\tvarscan_DP/" \
            > ${params.sample_name}_varscan.txt
        """
}


/*-------------------------------------------------------------------*
 *  Combine results and compare
 *-------------------------------------------------------------------*/

process combine_vcfs {
    /* 
     * Once all variant calling steps have been run, merge variant 
     * calls into single dataframe and save as CSV
     */
    publishDir "${params.outdir}/combined", mode: "copy"

    input:
        file(mutect_variants) from processed_mutect
        file(sinvict_variants) from processed_sinvict
        file(varscan_variants) from processed_varscan

    output:
        file("${params.sample_name}_combined_results.txt") into combined_vcfs

    script:
        """
        #!/usr/bin/env python
        import pandas as pd

        # load all files into dataframes with variant ID as index
        df1 = pd.read_table("$mutect_variants", index_col=0)
        df2 = pd.read_table("$sinvict_variants", index_col=0)
        df3 = pd.read_table("$varscan_variants", index_col=0)

        # combine dataframes based on index
        main = pd.concat([df1, df2, df3], axis=1)

        # write to file
        with open("${params.sample_name}_combined_results.txt", 'w') as f:
            main.to_csv(f)
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
