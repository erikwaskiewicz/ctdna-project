#!/bin/bash
#SBATCH --partition=highmem


fastq=$1
filename=$2

bwa mem /home/x.e.w.01/ref/human_g1k_v37.fasta $1 > "$2".bam
