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
LABEL maintainer="Erik Waskiewicz"

# Set the working directory & copy across current directory
WORKDIR /var/app
COPY . /var/app
ENV PYTHONDONTWRITEBYTECODE=true

# Install awk and ps for use with Nextflow reports, 
# Install build essentials and git for SiNVICT
RUN apt-get update \
 && apt-get install --yes --no-install-recommends \
        build-essential \
        gawk \
        git \
        procps \
# Remove apt cache after running to reduce image size
 && rm -rf /var/lib/apt/lists/*

# Make the conda environment
RUN conda install --yes --freeze-installed \
        bioconda::bam-readcount=0.8 \
        conda-forge::pandas=0.25.1 \
        conda-forge::numpy=1.17.2 \
        python=3.7 \
# Clean up conda env
 && conda clean -afy \
 && find /opt/conda/ -follow -type f -name '*.a' -delete \
 && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
 && find /opt/conda/ -follow -type f -name '*.js.map' -delete \
# Install SiNVICT
 && git clone https://github.com/sfu-compbio/sinvict.git \
 && make -C /var/app/sinvict
