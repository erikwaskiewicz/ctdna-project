## ctdna-sinvict Docker image

![Docker Cloud Automated build](https://img.shields.io/docker/cloud/automated/erikwaskiewicz/ctdna-sinvict?style=flat-square)
![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/erikwaskiewicz/ctdna-sinvict?style=flat-square)

Docker image for running the SiNVICT variant caller as part of the ctDNA variant caller comparison pipeline ([https://github.com/sfu-compbio/sinvict](https://github.com/sfu-compbio/sinvict)).

Based on the Miniconda3 base image (v4.7.10), with these packages installed through conda:

- bam-readcount v0.8 (requirement for SiNVICT)
- pandas v0.25.1 (to produce VCF output files)
- numpy v1.17.2 (to produce VCF output files)
- Python v3.7 (to produce VCF output files)

**Note**: pandas and numpy are installed through conda-forge to reduce image size due to large license file ([see GitHub issue here](https://github.com/conda-forge/numpy-feedstock/issues/84))

SiNVICT isn't available through Conda and is installed through git and make.

The following requirements are installed through apt-get:

- build-essential (required to build SiNVICT)
- git (required to build SiNVICT)
- gawk (Required for Nextflow reports)
- procps (Required for Nextflow reports)

Requires the custom script, `run_sinvict.py`, to run SiNVICT:

- Runs bam-readcount to generate a pileup file required by SiNVICT
- Runs SiNVICT (uses FASTA reference otherwise SiNVICT will call loads of variants)
- Compiles SiNVICT output into VCF format

GitHub: [https://github.com/erikwaskiewicz/ctdna-project/](https://github.com/erikwaskiewicz/ctdna-project/)

DockerHub: [https://cloud.docker.com/repository/docker/erikwaskiewicz/ctdna-sinvict](https://cloud.docker.com/repository/docker/erikwaskiewicz/ctdna-sinvict)
