import sys
import subprocess
import pandas as pd

# =====================================================================
# load in variables and run command line tools
# =====================================================================

input_bam = sys.argv[1]
ref_fasta = sys.argv[2]
output_file = sys.argv[3]

# parse sample name from output file
sample_name = '_'.join(output_file.split('_')[:-1])

# make directory to save bam-readcount output, sinvict requires that it is only file in folder
subprocess.call('mkdir /home/analysis_in', shell=True)

# run bam-readcount, suppress warnings to 1 (otherwise will be one per read)
# if fasta ref is not included, sinvict will call lots of Ns
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


# =====================================================================
# process sinvict output
# sinvict outputs 6 tsv files, one for each filter, we need to merge and convert to VCF
# =====================================================================

# make list of sinvict output files and the filter that was applied to each
filter_list = (
    ('calls_level1.sinvict', '1-poisson_model'),
    ('calls_level2.sinvict', '2-min_read_depth'),
    ('calls_level3.sinvict', '3-strand_bias'),
    ('calls_level4.sinvict', '4-average_pos'),
    ('calls_level5.sinvict', '5-signal_to_noise'),
    ('calls_level6.sinvict', '6-homopolymer_regions')
)

# label columns of sinvict output, 
# capitals are ready to go into final vcf, lower case need processing first
sinvict_out_cols = [
    'chr', 'position', 'sample', 'REF', 'DP', 'ALT', 'SVR', 'percent', 
    'mapped_plus', 'mapped_minus', 'ARP', 'SVT'
]

# make empty dataframe to store each iteration
all_calls = pd.DataFrame(columns=sinvict_out_cols.append('SVF'))

# loop through each output file and add a column for filter type
for filepath, filter_tag in filter_list:
    df = pd.read_csv(
        f'/home/{filepath}', 
        header=None,
        index_col=False,
        sep='\t', 
        names=sinvict_out_cols
    )
    df['SVF'] = filter_tag

    # add to main dataframe
    all_calls = pd.concat([all_calls, df], ignore_index=True)


# =====================================================================
# convert output to VCF
# =====================================================================

# make VCF header
vcf_header = [
    '##fileformat=VCFv4.2',
    '##FILTER=<ID=PASS,Description="All filters passed">',
    f'##reference={ref_fasta}',
    '##ALT=<ID=*,Description="Represents allele(s) other than observed.">',
    '##FORMAT=<ID=SVF,Number=.,Type=String,Description="sinvict filter">',
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

# process sinvict data before importing into VCF
all_calls = all_calls.sort_values(by='position')

all_calls['#CHROM'] = all_calls['chr'].map(str)
all_calls['POS'] = all_calls['position'].map(str)

all_calls['MPS'] = all_calls['mapped_plus'].str.strip('+:')
all_calls['MMS'] = all_calls['mapped_minus'].str.strip('-:')
all_calls['SVVAF'] = all_calls['percent'].div(10).round(3)

all_calls['ID'] = '.'
all_calls['QUAL'] = '.'
all_calls['FILTER'] = '.'
all_calls['INFO'] = '.'

# make format field
all_calls['FORMAT'] = 'SVF:SVVAF:DP:SVT:SVR:MPS:MMS:ARP'

all_calls[sample_name] = all_calls[
    ['SVF', 'SVVAF', 'DP', 'SVT', 'SVR', 'MPS', 'MMS', 'ARP']
].apply(lambda x: ':'.join(x.map(str)), axis=1)

# write VCF body
all_calls = all_calls[
    ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample_name]
]
all_calls.to_csv(output_file, index=False, sep='\t', mode='a')
