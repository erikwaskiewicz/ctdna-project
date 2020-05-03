## ctdna-freebayes Docker image

![Docker Cloud Automated build](https://img.shields.io/docker/cloud/automated/erikwaskiewicz/ctdna-freebayes?style=flat-square)
![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/erikwaskiewicz/ctdna-freebayes?style=flat-square)

Docker image for running the freebayes variant caller as part of the ctDNA variant caller comparison pipeline ([https://github.com/ekg/freebayes](https://github.com/ekg/freebayes)).

Based on the Miniconda3 base image (v4.7.10), with these packages installed through conda:
- freebayes v1.3.1-0

Also has awk and ps installed through apt-get, as these are required by Nextflow to generate pipeline execution reports and aren't included in the base image.

GitHub: https://github.com/erikwaskiewicz/ctdna-project/

DockerHub: https://cloud.docker.com/repository/docker/erikwaskiewicz/ctdna-freebayes
