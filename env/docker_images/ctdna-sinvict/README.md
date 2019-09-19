## ctdna-sinvict Docker image

![Docker Cloud Automated build](https://img.shields.io/docker/cloud/automated/erikwaskiewicz/ctdna-sinvict?style=flat-square)
![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/erikwaskiewicz/ctdna-sinvict?style=flat-square)

Docker image for running the SiNVICT variant caller as part of the ctDNA Nextflow pipeline (https://github.com/sfu-compbio/sinvict).

Based on the Miniconda3 base image (v4.7.10), with these packages installed through conda:
- bam-readcount v0.8 (requirement for SiNVICT)
- bcftools v1.9 (to produce VCF output files)

SiNVICT isn't available through Conda and is installed through git and make.

The following requirements are installed through apt-get:

- build-essential (required to build SiNVICT)
- git (required to build SiNVICT)
- gawk (Required for Nextflow reports)
- procps (Required for Nextflow reports)

GitHub: https://github.com/erikwaskiewicz/ctdna-project/

DockerHub: https://cloud.docker.com/repository/docker/erikwaskiewicz/ctdna-sinvict
