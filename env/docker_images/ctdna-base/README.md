## ctdna-base Docker image

![Docker Cloud Automated build](https://img.shields.io/docker/cloud/automated/erikwaskiewicz/ctdna-base?style=flat-square)
![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/erikwaskiewicz/ctdna-base?style=flat-square)

Base Docker image for the ctDNA Nextflow pipeline, to be used for all processes unless a specific Docker image has been specified.

Based on the Miniconda3 base image (v4.7.10), with these packages installed through conda:
- picard v2.20.6
- samtools v1.9

Also has awk and ps installed through apt-get, as these are required by Nextflow to generate pipeline execution reports and aren't included in the base image.

GitHub: https://github.com/erikwaskiewicz/ctdna-project/

DockerHub: https://cloud.docker.com/repository/docker/erikwaskiewicz/ctdna-base
