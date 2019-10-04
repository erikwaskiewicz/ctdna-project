# =====================================================================
#              erikwaskiewicz/ctdna-sinvict Dockerfile
#
# For use in ctDNA pipeline project - image for running the SiNVICT
# variant caller
# <https://github.com/erikwaskiewicz/ctdna-project>
#
# =====================================================================

# Use miniconda parent image - debian based
FROM continuumio/miniconda3:4.7.10
MAINTAINER Erik Waskiewicz

# Install awk and ps for use with Nextflow reports, 
# Install build essentials and git for SiNVICT
RUN apt-get update \
 && apt-get install --yes --no-install-recommends \
    build-essential \
    gawk \
    git \
    procps \
# remove apt cache after running to reduce image size
 && rm -rf /var/lib/apt/lists/*

# Set the working directory & copy across current directory
WORKDIR /home
COPY . /home

# Update conda, make the conda environment and install sinvict
RUN conda update -n base -c defaults conda \
 && conda install \
    bioconda::bam-readcount=0.8 \
    conda-forge::pandas=0.25.1 \
    conda-forge::numpy=1.17.2 \
    python=3.7 \
 && git clone https://github.com/sfu-compbio/sinvict.git \
 && make -C /home/sinvict
