## ctdna-base Docker image

Base Docker image for the ctDNA Nextflow pipeline, to be used for all processes unless a specific Docker image has been specified.
See https://github.com/erikwaskiewicz/ctdna-project/ for full details.

Based on the Miniconda3 base image (v4.7.10), with these packages installed through conda:
- picard v2.20.6
- samtools v1.9

Also has awk and ps installed through apt-get, as these are required by Nextflow to generate pipeline execution reports and aren't included in the base image.
