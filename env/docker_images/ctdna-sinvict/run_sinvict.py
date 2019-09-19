import sys
import subprocess
import pandas as pd

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
# merge all output, add column for which file variants came from
filter_list = (
    ('calls_level1.sinvict', '1-poisson_model'),
    ('calls_level2.sinvict', '2-min_read_depth'),
    ('calls_level3.sinvict', '3-strand_bias'),
    ('calls_level4.sinvict', '4-average_pos'),
    ('calls_level5.sinvict', '5-signal_to_noise'),
    ('calls_level6.sinvict', '6-homopolymer_regions')
)

# make empty dataframe to store each iteration
all_calls = pd.DataFrame(
    columns=['#CHROM', 'POS', 'sample', 'REF', 'DP', 'ALT', 'allele_count', 'vaf', 'a', 'b', 'c', 'type', 'sinvict_filter']
)

# loop through each output file and add a column for filter type
for item in filter_list:
    filepath, filter_tag = item

    df = pd.read_csv(
        f'/home/{filepath}', 
        header=None,
        index_col=False,
        sep='\t', 
        names=['#CHROM', 'POS', 'sample', 'REF', 'DP', 'ALT', 'allele_count', 'vaf', 'a', 'b', 'c', 'type']
    )

    df['sinvict_filter'] = filter_tag

    all_calls = pd.concat([all_calls, df], ignore_index=True)


# TODO convert merged data into VCF
# TODO ref and alt are wrong way round
all_calls = all_calls.sort_values(by='POS')
all_calls['a'] = all_calls['a'].str.strip('+:')
all_calls['b'] = all_calls['b'].str.strip('-:')
all_calls['#CHROM'] = all_calls['#CHROM'].map(str)
all_calls['POS'] = all_calls['POS'].map(str)
all_calls['ID'] = '.'
all_calls['QUAL'] = '.'
all_calls['FILTER'] = '.'
all_calls['INFO'] = '.'

# make format field
all_calls['FORMAT'] = 'SVF:SVVAF:DP:SVT:SVR:MPS:MMS:ARP'
all_calls[sample_name] = all_calls[[
    'sinvict_filter', 'vaf', 'DP', 'type', 'allele_count', 'a', 'b', 'c'
]].apply(lambda x: ':'.join(x.map(str)), axis=1)

# write to file
 # TODO add VCF header and format field descriptions
with open(output_file, 'w') as f:
    f.writelines(
f'''##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##reference={ref_fasta}
##ALT=<ID=*,Description="Represents allele(s) other than observed.">
##FORMAT=<ID=SVF,Type=String,Description="sinvict filter">
''')

#all_calls = all_calls.astype(str)
all_calls = all_calls[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample_name]]
all_calls.to_csv(output_file, index=False, sep='\t', mode='a')


'''
format field descriptions
SVF   - sinvict filter
SVVAF - sinvict vaf
SVT   - sinvict type
SVR   - sinvict reads supporting
MPS   - mapped_plus_strand
MMS   - mapped_minus_strand
ARP   - average read pos
'''

'''
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##reference=file://../test_data/Homo_sapiens_chr22_assembly19.fasta
##contig=<ID=22,length=51304566>
##ALT=<ID=*,Description="Represents allele(s) other than observed.">
##FORMAT=<ID=SVF,Type=String,Description="sinvict filter">
'''
