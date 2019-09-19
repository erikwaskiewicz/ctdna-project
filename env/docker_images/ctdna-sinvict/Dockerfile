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
# build essentials and git for installing sinvict
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
 && conda install -c bioconda bam-readcount=0.8 \
 && conda install -c bioconda bcftools=1.9 \
 && git clone https://github.com/sfu-compbio/sinvict.git \
 && make -C /home/sinvict
