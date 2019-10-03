# =====================================================================
#              erikwaskiewicz/ctdna-base Dockerfile
#
# For use in ctDNA pipeline project - base image to be used by 
# pipeline if a specific image isn't being used
# <https://github.com/erikwaskiewicz/ctdna-project>
#
# =====================================================================

# Use miniconda parent image - debian based
FROM continuumio/miniconda3:4.7.10
MAINTAINER Erik Waskiewicz

# Install awk and ps for use with Nextflow reports
RUN apt-get update \
 && apt-get install --yes --no-install-recommends \
    gawk \
    procps \
# remove apt cache after running to reduce image size
 && rm -rf /var/lib/apt/lists/*

# Set the working directory & copy across current directory
WORKDIR /home
COPY env/docker_images/ctdna-base /home

# Update conda and make the conda environment
RUN conda update -n base -c defaults conda \
 && conda install -c bioconda samtools=1.9 \
 && conda install -c bioconda picard=2.20.6
