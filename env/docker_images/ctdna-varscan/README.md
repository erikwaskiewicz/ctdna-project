## ctdna-varscan Docker image

![Docker Cloud Automated build](https://img.shields.io/docker/cloud/automated/erikwaskiewicz/ctdna-varscan?style=flat-square)
![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/erikwaskiewicz/ctdna-varscan?style=flat-square)

Docker image for running the VarScan2 variant caller as part of the ctDNA variant caller comparison pipeline ([http://varsan.sourceforge.net/](http://varsan.sourceforge.net/)).

Based on the Miniconda3 base image (v4.7.10), with these packages installed through conda:
- varscan v2.4.3
- samtools v1.9

Also has awk and ps installed through apt-get, as these are required by Nextflow to generate pipeline execution reports and aren't included in the base image.

GitHub: https://github.com/erikwaskiewicz/ctdna-project/

DockerHub: https://cloud.docker.com/repository/docker/erikwaskiewicz/ctdna-varscan
