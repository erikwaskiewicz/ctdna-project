'''
run_sinvict.py

Converts SiNVICT output into VCF
To be run within ctdna-sinvict Docker container:
https://cloud.docker.com/repository/docker/erikwaskiewicz/ctdna-sinvict
'''

import sys
import subprocess
import pandas as pd
import csv


# =====================================================================
# load in variables and run command line tools
# =====================================================================

input_bam = sys.argv[1]
ref_fasta = sys.argv[2]
output_file = sys.argv[3]
roi_bed = sys.argv[4]

# parse sample name from output file
sample_name = '_'.join(output_file.split('_')[:-1])

# change into Nextflow work directory
subprocess.call('cd $NXF_WORK', shell=True)

# make directory to save bam-readcount output, sinvict requires that it is only file in folder
subprocess.call('mkdir analysis_in', shell=True)

# run bam-readcount, suppress warnings to 1 (otherwise will be one per read)
# if fasta ref is not included, sinvict will call lots of Ns
subprocess.call(
    f'''
    bam-readcount \
        --reference-fasta {ref_fasta} \
        --site-list {roi_bed} \
        --min-mapping-quality 20 \
        --min-base-quality 20 \
        --max-warnings 1 \
        {input_bam} \
        > analysis_in/{sample_name}.txt
    ''',
    shell=True
)

# run sinvict
subprocess.call(
    '/var/app/sinvict/sinvict -t analysis_in -o .',
    shell=True
)


# =====================================================================
# process sinvict output
# sinvict outputs 6 tsv files, one for each filter, we need to merge and convert to VCF
# =====================================================================

# make list of sinvict output files and the filter that was applied to each
filter_list = (
    ('calls_level1.sinvict', '1'),
    ('calls_level2.sinvict', '2'),
    ('calls_level3.sinvict', '3'),
    ('calls_level4.sinvict', '4'),
    ('calls_level5.sinvict', '5'),
    ('calls_level6.sinvict', '6')
)

# make empty dataframe to store each iteration
all_calls_dict = {}

# loop through each output file and add a column for filter type
for filepath, filter_tag in filter_list:
    with open(filepath, 'r') as f:
        csvreader = csv.reader(f, delimiter='\t')
        for line in csvreader:
            chrom, position, sample, ref, dp, alt, svr, percent, mapped_plus, mapped_minus, arp, svt = line
            unique_id = f'{chrom}_{position}_{ref}_{alt}'
            if unique_id in all_calls_dict:
                # add filter to filter column
                all_calls_dict[unique_id]['tag'] += f',{filter_tag}'
            else:
                # reformat del/ins in alt field
                if '-' in alt:
                    ref = alt.strip('-')
                    alt = '<DEL>'
                elif '+' in alt:
                    # TODO - this just strips the illegal characters - need to make sure genotypes are correct
                    alt = alt.strip('+')

                # make new object
                all_calls_dict[unique_id] = {
                    'chr': str(chrom),
                    'position': str(position),
                    'sample': sample,
                    'ref': ref,
                    'depth': dp,
                    'alt': alt,
                    'svr': svr,
                    'percent': round(float(percent) / 100, 3),
                    'mapped_plus': str(mapped_plus).strip('+:'),
                    'mapped_minus': str(mapped_minus).strip('-:'),
                    'arp': arp,
                    'svt': svt,
                    'tag': filter_tag
                }


# =====================================================================
# convert output to VCF
# =====================================================================

# process sinvict data before importing into VCF
all_calls = pd.DataFrame.from_dict(all_calls_dict, orient='index')
all_calls = all_calls.sort_values(['chr', 'position'])

formatted_calls = pd.DataFrame()
formatted_calls['#CHROM'] = all_calls['chr']
formatted_calls['POS'] = all_calls['position']
formatted_calls['ID'] = '.'
formatted_calls['REF'] = all_calls['ref']
formatted_calls['ALT'] = all_calls['alt']
formatted_calls['QUAL'] = '.'
formatted_calls['FILTER'] = '.'
formatted_calls['INFO'] = '.'

# make format field
formatted_calls['FORMAT'] = 'SVF:SVVAF:DP:SVT:SVR:MPS:MMS:ARP'
formatted_calls[sample_name] = all_calls[
    ['tag', 'percent', 'depth', 'svt', 'svr', 'mapped_plus', 'mapped_minus', 'arp']
].apply(lambda x: ':'.join(x.map(str)), axis=1)

# make VCF header
vcf_header = [
    '##fileformat=VCFv4.2',
    '##FILTER=<ID=PASS,Description="All filters passed">',
    f'##reference={ref_fasta}',
    '##ALT=<ID=*,Description="Represents allele(s) other than observed.">',
    '##FORMAT=<ID=SVF,Number=.,Type=String,Description="sinvict filter - 1-poisson_model;2-min_read_depth;3-strand_bias;4-average_pos;5-signal_to_noise;6-homopolymer_regions">',
    '##FORMAT=<ID=SVVAF,Number=.,Type=String,Description="sinvict vaf">',
    '##FORMAT=<ID=DP,Number=.,Type=String,Description="depth">',
    '##FORMAT=<ID=SVT,Number=.,Type=String,Description="sinvict type (germline/somatic)">',
    '##FORMAT=<ID=SVR,Number=.,Type=String,Description="sinvict reads supporting">',
    '##FORMAT=<ID=MPS,Number=.,Type=String,Description="mapped_plus_strand">',
    '##FORMAT=<ID=MMR,Number=.,Type=String,Description="mapped_plus_strand">',
    '##FORMAT=<ID=ARP,Number=.,Type=String,Description="average read pos">'
]

# write VCF header to file
with open(output_file, 'w') as f:
    f.writelines(line + '\n' for line in vcf_header)

# write data to file
formatted_calls.to_csv(output_file, index=False, sep='\t', mode='a')
