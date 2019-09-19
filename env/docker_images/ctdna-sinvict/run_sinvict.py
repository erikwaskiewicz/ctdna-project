import sys
import subprocess

input_bam = sys.argv[1]
ref_fasta = sys.argv[2]
output_file = sys.argv[3]

# parse sample name from output file
sample_name = '_'.join(output_file.split('_')[:-1])

# make directory to save bam-readcount output, sinvict requires that it is only file in folder
subprocess.call('mkdir /home/analysis_in', shell=True)

# run bam-readcount, suppress warnings to 1 (otherwise will be one per read)
subprocess.call(
    f'''
    bam-readcount \
        --reference-fasta {ref_fasta} \
        --min-mapping-quality 20 \
        --min-base-quality 20 \
        --max-warnings 1 \
        {input_bam} \
        > /home/analysis_in/{sample_name}.txt
    ''',
    shell=True
)

# run sinvict
subprocess.call(
    '/home/sinvict/sinvict -t /home/analysis_in -o /home',
    shell=True
)

# process sinvict output
# TODO merge all output, add column for which file variants came from

# convert sinvict output to vcf
subprocess.call(
    f'''
    bcftools convert \
        --tsv2vcf /home/calls_level1.sinvict \
        --columns CHROM,POS,-,-,-,AA,-,-,-,-,-,- \
        --samples {sample_name} \
        --fasta-ref {ref_fasta} \
        --output {output_file}
    ''',
    shell=True
)
