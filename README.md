# ctdna-nextflow

A comparison of variant callers for a ctDNA pipeline

![TravisCI](https://api.travis-ci.com/erikwaskiewicz/ctdna-project.svg?token=6Kdapt3JHU3GMffhB2Mx&branch=master)

## Introduction

A comparison of variant callers for the purpose of calling ctDNA variants. Takes an input BAM file and runs many variant callers in parallel, then combines the outputs in the end. Each variant caller runs within its own Docker container .

This was written as part of my MSc Clinical Science (Clinical Bioinformatics) project at the University of Manchester.

## Requirements and setup

For detailed requirements and setup instructions, see the [Wiki](https://github.com/erikwaskiewicz/ctdna-project/wiki/Setup)

## Running

To run with the included test data:

`nextflow ctdna-snv-caller-comparison.nf -profile test,docker_conf`

## Tests

After running the test data above, run the PyTest script: `pytest --verbose`
