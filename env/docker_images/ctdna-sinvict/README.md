## ctdna-sinvict Docker image

![Docker Cloud Automated build](https://img.shields.io/docker/cloud/automated/erikwaskiewicz/ctdna-sinvict?style=flat-square)
![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/erikwaskiewicz/ctdna-sinvict?style=flat-square)

Docker image for running the SiNVICT variant caller as part of the ctDNA Nextflow pipeline (https://github.com/sfu-compbio/sinvict).

Includes bam-readcount for pre-processing BAMs.

The following requirements are installed through apt-get:

- build-essential
- git
- software-properties-common
- cmake
- libncurses-dev
- zlib1g-dev
- patch

GitHub: https://github.com/erikwaskiewicz/ctdna-project/

DockerHub: https://cloud.docker.com/repository/docker/erikwaskiewicz/ctdna-sinvict
